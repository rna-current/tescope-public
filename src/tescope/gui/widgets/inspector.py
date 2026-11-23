"""Locus Inspector: one locus in full context.

Display-only: inspecting, navigating, or highlighting never changes
the analysis, its denominators, or any export. Coordinates shown to
the user are 1-based inclusive (the machine-readable locus ID stays
0-based half-open and is shown verbatim, labeled).
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from ...core.coordinates.mapping import CopyAlignment, consensus_interval_to_genomic
from ...core.tracks.base import TrackKind


def _esc(s) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


class LocusInspector(QDialog):
    highlight_requested = Signal(str)  # locus_id

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Locus Inspector")
        self.resize(560, 640)
        self.setModal(False)
        self.result = None
        self._ids: list[str] = []
        self._pos = 0

        lay = QVBoxLayout(self)
        nav = QHBoxLayout()
        self.prev_btn = QPushButton("◀ Previous")
        self.next_btn = QPushButton("Next ▶")
        self.pos_label = QLabel("")
        self.pos_label.setAlignment(Qt.AlignCenter)
        nav.addWidget(self.prev_btn)
        nav.addWidget(self.pos_label, 1)
        nav.addWidget(self.next_btn)
        lay.addLayout(nav)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        body_host = QWidget()
        body_lay = QVBoxLayout(body_host)
        self.body = QLabel()
        self.body.setWordWrap(True)
        self.body.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.body.setAlignment(Qt.AlignTop)
        body_lay.addWidget(self.body)
        body_lay.addStretch(1)
        scroll.setWidget(body_host)
        lay.addWidget(scroll, 1)

        bottom = QHBoxLayout()
        bottom.addWidget(QLabel("Neighbor flank (bp):"))
        self.flank_spin = QSpinBox()
        self.flank_spin.setRange(0, 1_000_000)
        self.flank_spin.setValue(5000)
        self.flank_spin.setSingleStep(1000)
        self.flank_spin.setToolTip(
            "Genomic flank for listing user-supplied neighboring "
            "features (display only).")
        bottom.addWidget(self.flank_spin)
        bottom.addStretch(1)
        self.highlight_btn = QPushButton("Highlight in heatmap")
        self.highlight_btn.setToolTip(
            "Mark this locus's row in the per-copy heatmap "
            "(display-only marker; no value or order changes).")
        bottom.addWidget(self.highlight_btn)
        close_btn = QPushButton("Close")
        bottom.addWidget(close_btn)
        lay.addLayout(bottom)

        self.prev_btn.clicked.connect(lambda: self._step(-1))
        self.next_btn.clicked.connect(lambda: self._step(+1))
        close_btn.clicked.connect(self.close)
        self.flank_spin.valueChanged.connect(lambda _v: self._render())
        self.highlight_btn.clicked.connect(
            lambda: self.highlight_requested.emit(self._current_id()))

    # -- navigation ---------------------------------------------------------

    def show_locus(self, result, ordered_ids: list[str], current: str) -> None:
        self.result = result
        self._ids = ordered_ids
        self._pos = ordered_ids.index(current) if current in ordered_ids else 0
        self._render()

    def _current_id(self) -> str:
        return self._ids[self._pos] if self._ids else ""

    def _step(self, delta: int) -> None:
        if not self._ids:
            return
        self._pos = max(0, min(len(self._ids) - 1, self._pos + delta))
        self._render()

    # -- rendering ----------------------------------------------------------

    def _row(self):
        from ...core.loci.table import locus_ids

        ids = locus_ids(self.result.copies)
        pos = np.where(ids.to_numpy() == self._current_id())[0]
        if not pos.size:
            return None, None
        i = int(pos[0])
        return i, self.result.copies.iloc[i]

    def _render(self) -> None:  # noqa: C901 - one linear report
        if self.result is None or not self._ids:
            return
        i, row = self._row()
        if row is None:
            self.body.setText("Locus is no longer in the eligible "
                              "population.")
            return
        lid = self._current_id()
        self.pos_label.setText(
            f"{self._pos + 1} of {len(self._ids)} (current table order)")
        self.prev_btn.setEnabled(self._pos > 0)
        self.next_btn.setEnabled(self._pos < len(self._ids) - 1)
        r = self.result
        h: list[str] = []

        def sec(title):
            h.append(f"<h3 style='margin-bottom:2px'>{_esc(title)}</h3>")

        def kv(k, v):
            h.append(f"<b>{_esc(k)}:</b> {v}<br>")

        gs, ge = int(row.genome_start), int(row.genome_end)
        sec("Identity")
        kv("Locus ID", f"<code>{_esc(lid)}</code> "
                       "<i>(machine key, 0-based half-open)</i>")
        kv("Assembly", _esc(r.config.assembly))
        kv("Genomic position (1-based)",
           f"{_esc(row.chrom)}:{gs + 1:,}–{ge:,} ({_esc(row.strand)} strand)")
        kv("Genomic length", f"{ge - gs:,} bp")
        prov = r.provenance or {}
        kv("Annotation source",
           _esc(prov.get("source", "?"))
           + (f", retrieved {_esc(prov.get('retrieved_at', ''))}"
              if prov.get("retrieved_at") else ""))

        sec("Classification (from the annotation record)")
        kv("Repeat", f"{_esc(row.repeat_name)} "
                     f"({_esc(row.repeat_class)}/{_esc(row.repeat_family)})")
        kv("Completeness class", _esc(row.completeness))
        cs, ce = int(row.cons_start), int(row.cons_end)
        kv("Consensus footprint (1-based)",
           f"{cs + 1:,}–{ce:,} of {r.cons_length:,} bp "
           f"({float(row.consensus_coverage):.1%} of model)")
        rmf = int((r.copies['rm_id'] == row.rm_id).sum()) \
            if 'rm_id' in r.copies else 1
        kv("RepeatMasker record",
           f"rm_id {int(row.rm_id)}; {rmf} analyzed record(s) share this "
           "ID (fragment linkage — records are never merged)")
        h.append("<i>One row = one RepeatMasker record; this is not "
                 "asserted to be a complete biological insertion.</i><br>")

        sec("Signal at this locus")
        table = r.locus_table()
        trow = table[table["locus_id"] == lid].iloc[0]
        for t in r.tracks:
            if t.kind is TrackKind.CONTINUOUS and t.matrix is not None:
                mv = trow.get(f"mean_{t.name}")
                ob = trow.get(f"obs_bins_{t.name}")
                if mv is not None and np.isfinite(mv):
                    kv(t.name, f"mean {mv:.4g} over {int(ob)} observed "
                               "bin(s) (input units, as supplied)")
                else:
                    kv(t.name, "<i>no observed data (missing — not a "
                               "measured zero)</i>")
            elif t.kind is TrackKind.INTERVAL:
                pv = trow.get(f"peaks_{t.name}")
                if pv is None or (isinstance(pv, float) and np.isnan(pv)):
                    kv(t.name, "<i>unavailable (chromosome absent from "
                               "track)</i>")
                else:
                    kv(t.name, f"{int(pv)} overlapping interval(s)")
        if r.locus_estimates is not None:
            est = r.locus_estimates
            tot = trow.get("est_total")
            if tot is not None and np.isfinite(tot):
                parts = [f"total {tot:g} {est.units}"]
                uq = trow.get("est_unique")
                if uq is not None and np.isfinite(uq):
                    parts.append(f"unique component {uq:g}")
                else:
                    parts.append("unique component unavailable (not zero)")
                sc = trow.get("est_score")
                if sc is not None and np.isfinite(sc):
                    parts.append(f"score {sc:g} (native semantics)")
                kv(f"Upstream estimate ({est.assay}"
                   + (f", {est.tool}" if est.tool else "") + ")",
                   "; ".join(parts)
                   + " — locus-level estimate, no positional resolution")
            else:
                kv("Upstream estimate",
                   "<i>no row for this locus (unavailable — not zero)</i>")

        if r.features is not None and r.features.n_total:
            sec("Consensus features on this locus (projected)")
            h.append("<i>Projection is affine within the RepeatMasker hit "
                     "endpoints — approximate, not gap-aware or "
                     "base-accurate.</i><br>")
            aln = CopyAlignment(gs, ge, cs, ce, str(row.strand))
            for f in r.features.all_features():
                res = consensus_interval_to_genomic(aln, f.start, f.end)
                tag = " <i>(user)</i>" if f.origin == "user-defined" else ""
                if res is None:
                    kv(f.name + tag,
                       "<i>outside this copy's consensus footprint — no "
                       "verified correspondence; shown on the consensus "
                       "view only</i>")
                else:
                    g0, g1, partial = res
                    kv(f.name + tag,
                       f"projected to {_esc(row.chrom)}:{g0 + 1:,}–{g1:,}"
                       + (" <b>(partial — feature extends beyond this "
                          "copy)</b>" if partial else ""))

        if r.genome_features:
            flank = int(self.flank_spin.value())
            sec(f"User-supplied genomic features within ±{flank:,} bp")
            any_hit = False
            for gset in r.genome_features:
                for f in gset.overlapping(str(row.chrom), gs, ge,
                                          flank=flank):
                    any_hit = True
                    rel = ("overlaps locus" if f.start < ge and f.end > gs
                           else (f"{gs - f.end:,} bp upstream of locus "
                                 "start" if f.end <= gs
                                 else f"{f.start - ge:,} bp downstream"))
                    kv(f.label or f.feature_id,
                       f"{_esc(f.chrom)}:{f.start + 1:,}–{f.end:,} "
                       f"({_esc(f.strand)}) [{_esc(f.type)}] — {rel}")
            if not any_hit:
                h.append("<i>none in this window (absence from the "
                         "supplied annotation file — not proof of "
                         "biological absence)</i><br>")

        self.body.setText("".join(h))
