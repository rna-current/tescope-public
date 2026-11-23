"""Interactive consensus browser: structural tests (offscreen Qt).

Verifies display-only navigation (zoom/pan change x limits ONLY),
synchronized axes, zoom-dependent sequence rendering, capability gating
of sequence controls, and scientific invariance of the whole sequence /
interaction layer. No pixel comparisons — structure and numbers.
"""

from __future__ import annotations

import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402

from tescope.config.model import (  # noqa: E402
    AnalysisConfig,
    RepeatSelection,
    SequenceConfig,
    TrackConfig,
)
from tescope.core.analysis import run_analysis  # noqa: E402
from tescope.core.sequence.cache import store_cached  # noqa: E402
from tescope.core.sequence.coords import RangeError  # noqa: E402
from tescope.core.sequence.model import (  # noqa: E402
    STATUS_VERIFIED,
    ConsensusSequence,
    SequenceProvenance,
    SequenceResolution,
)

CONS_LEN = 1000  # golden TST1 model length
SEQ1000 = ("ACGT" * 250)


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture(scope="module")
def seq_fasta(tmp_path_factory):
    p = tmp_path_factory.mktemp("seq") / "tst1.fa"
    p.write_text(f">TST1#LINE/TST synthetic test consensus\n{SEQ1000}\n")
    return p


def _cfg(golden_dir, seq_path=None, **kw):
    c = AnalysisConfig(
        assembly="testasm",
        repeat=RepeatSelection(name="TST1"),
        tracks=[TrackConfig(file=golden_dir["bigwig"], name="sig_a")],
        **kw,
    )
    c.aggregation.bin_size = 100
    if seq_path is not None:
        c.sequence = SequenceConfig(enabled=True, source="fasta",
                                    path=str(seq_path))
    return c


@pytest.fixture()
def view(qapp, golden_dir, golden_store, seq_fasta):
    from tescope.gui.widgets.results import ResultsView

    res = run_analysis(_cfg(golden_dir, seq_fasta), store=golden_store)
    rv = ResultsView()
    rv.show_result(res)
    yield rv
    rv.navigator.detach()
    rv.deleteLater()


def _axis_windows(rv):
    return {
        tuple(round(v, 6) for v in ax.get_xlim())
        for ax in rv.profile_canvas.canvas.figure.axes
    }


def _seq_ax(rv):
    return next(a for a in rv.profile_canvas.canvas.figure.axes
                if hasattr(a, "_tescope_sequence"))


class TestNavigation:
    def test_initial_full_view(self, view):
        assert view.navigator.visible() == (0.0, float(CONS_LEN))

    def test_goto_range_synchronizes_every_panel(self, view):
        view.goto_edit.setText("101-140")
        view._go_to()
        assert view.navigator.visible() == (100.0, 140.0)
        assert _axis_windows(view) == {(100.0, 140.0)}  # ALL panels together

    def test_goto_single_position_centers_window(self, view):
        view.navigator.go_to("500")
        lo, hi = view.navigator.visible()
        assert lo < 499.5 < hi and (hi - lo) < 200

    def test_goto_bounds_validated(self, view):
        for bad in ("0", "1001", "abc", "900-1200"):
            with pytest.raises(RangeError):
                view.navigator.go_to(bad)

    def test_wheel_zoom_in_shrinks_window_around_cursor(self, view):
        nav = view.navigator
        nav.full_view()
        ax = view.profile_canvas.canvas.figure.axes[0]
        nav._on_scroll(SimpleNamespace(inaxes=ax, xdata=250.0, button="up"))
        lo, hi = nav.visible()
        assert (hi - lo) < CONS_LEN and lo < 250.0 < hi
        # cursor's relative position is preserved (centered zoom)
        assert abs((250.0 - lo) / (hi - lo) - 0.25) < 1e-6

    def test_wheel_zoom_out_clamped_to_consensus(self, view):
        nav = view.navigator
        ax = view.profile_canvas.canvas.figure.axes[0]
        nav.set_view(100, 900)
        for _ in range(10):
            nav._on_scroll(SimpleNamespace(inaxes=ax, xdata=500.0,
                                           button="down"))
        assert nav.visible() == (0.0, float(CONS_LEN))  # never past the ends

    def test_minimum_span_enforced(self, view):
        view.navigator.set_view(500.0, 500.5)
        lo, hi = view.navigator.visible()
        assert hi - lo >= 10.0

    def test_full_view_resets(self, view):
        view.navigator.set_view(300, 350)
        view.full_btn.click()
        assert view.navigator.visible() == (0.0, float(CONS_LEN))


class TestSequenceTrackRendering:
    def test_letters_only_below_threshold(self, view):
        ax = _seq_ax(view)
        view.navigator.full_view()
        view.profile_canvas.canvas.draw()
        assert len(ax._tescope_letter_artists) == 0  # never 1000 letters
        view.navigator.set_view(100, 140)
        view.profile_canvas.canvas.draw()
        assert len(ax._tescope_letter_artists) == 40
        view.navigator.full_view()
        view.profile_canvas.canvas.draw()
        assert len(ax._tescope_letter_artists) == 0

    def test_letters_are_the_true_bases_1based(self, view):
        view.navigator.set_view(100, 140)
        arts = _seq_ax(view)._tescope_letter_artists
        for t in arts:
            x = t.get_position()[0]  # base center b0 + 0.5
            b0 = int(x)  # 0-based index
            assert t.get_text() == SEQ1000[b0]
        # explicit corner: axis [100,140) shows bases 101..140 (1-based)
        xs = sorted(t.get_position()[0] for t in arts)
        assert xs[0] == 100.5 and xs[-1] == 139.5

    def test_sequence_follows_toolbar_style_xlim_change(self, view):
        # direct xlim change (toolbar path) refreshes letters on draw
        ax = _seq_ax(view)
        view.profile_canvas.canvas.figure.axes[0].set_xlim(200, 230)
        view.profile_canvas.canvas.draw()
        assert len(ax._tescope_letter_artists) == 30

    def test_show_hide_sequence_track(self, view):
        assert view.seq_check.isChecked()
        view.seq_check.setChecked(False)
        assert not any(hasattr(a, "_tescope_sequence")
                       for a in view.profile_canvas.canvas.figure.axes)
        view.seq_check.setChecked(True)
        assert any(hasattr(a, "_tescope_sequence")
                   for a in view.profile_canvas.canvas.figure.axes)


class TestCursorReadout:
    def test_readout_reports_position_base_signal_copies(self, view):
        txt = view.navigator._readout(105, 1)
        assert "base 105 (1-based)" in txt
        assert SEQ1000[104] in txt.split("|")[1]
        assert "sig_a" in txt and "copies:" in txt

    def test_motion_event_updates_label(self, view):
        ax = view.profile_canvas.canvas.figure.axes[0]
        view.navigator._on_motion(
            SimpleNamespace(inaxes=ax, xdata=104.3, ydata=1.0))
        assert "base 105" in view.cursor_label.text()

    def test_click_places_marker_on_every_panel(self, view):
        ax = view.profile_canvas.canvas.figure.axes[0]
        n_axes = len(view.profile_canvas.canvas.figure.axes)
        view.navigator._on_click(
            SimpleNamespace(inaxes=ax, xdata=300.0, button=1))
        assert len(view.navigator._markers) == n_axes


class TestSearchAndExport:
    def test_search_zooms_to_first_hit(self, view):
        view.search_edit.setText("ACGTACGTACGT")
        view._search_sequence()
        assert view.hits_combo.count() > 0
        lo, hi = view.navigator.visible()
        h = view._search_hits[0]
        assert lo <= h.base_start - 1 and h.base_end <= hi

    def test_search_no_hits_message(self, view):
        view.search_edit.setText("GGGGGGGGGG")
        view._search_sequence()
        assert "No literal match" in view.cursor_label.text()

    def test_copy_visible_sequence_fasta(self, view, qapp):
        view.navigator.set_view(100.0, 140.0)
        view._copy_sequence()
        text = qapp.clipboard().text()
        assert text.startswith(">")
        assert "positions=101-140(1-based)" in text
        assert "".join(text.strip().split("\n")[1:]) == SEQ1000[100:140]


class TestCapabilityGating:
    def test_no_sequence_controls_disabled_with_tooltips(
            self, qapp, golden_dir, golden_store):
        from tescope.gui.widgets.results import ResultsView

        res = run_analysis(_cfg(golden_dir), store=golden_store)
        rv = ResultsView()
        rv.show_result(res)
        for w in (rv.seq_check, rv.copy_seq_btn, rv.fasta_btn, rv.search_edit):
            assert not w.isEnabled()
            assert "source" in w.toolTip()
        # navigation itself stays available (orthogonal to sequence)
        assert rv.full_btn.isEnabled() and rv.goto_edit.isEnabled()
        rv.navigator.go_to("101-140")
        assert rv.navigator.visible() == (100.0, 140.0)

    def test_incompatible_sequence_controls_disabled(
            self, qapp, golden_dir, golden_store, tmp_path):
        from tescope.gui.widgets.results import ResultsView

        bad = tmp_path / "bad.fa"
        bad.write_text(">TST1\n" + "ACGT" * 300 + "\n")  # 1200 != 1000
        res = run_analysis(_cfg(golden_dir, bad), store=golden_store)
        assert res.sequence.status == "incompatible"
        rv = ResultsView()
        rv.show_result(res)
        assert not rv.fasta_btn.isEnabled()
        assert "does not match" in rv.fasta_btn.toolTip()
        # no sequence axis is drawn for an incompatible model
        assert not any(hasattr(a, "_tescope_sequence")
                       for a in rv.profile_canvas.canvas.figure.axes)

    def test_cached_dfam_sequence_enables_offline(
            self, qapp, golden_dir, golden_store):
        from tescope.gui.widgets.results import ResultsView

        # simulate an earlier online lookup now used offline (cache-first)
        store_cached(
            "TST1", CONS_LEN,
            SequenceResolution(
                status=STATUS_VERIFIED, detail="from cache",
                sequence=ConsensusSequence(SEQ1000, SequenceProvenance(
                    provider="Dfam", accession="DFTEST")),
                provenance=SequenceProvenance(provider="Dfam")),
            cache_dir=golden_store.cache_dir)
        cfg = _cfg(golden_dir)
        cfg.sequence = SequenceConfig(enabled=True, source="dfam")
        res = run_analysis(cfg, store=golden_store)
        assert res.sequence is not None and res.sequence.verified
        rv = ResultsView()
        rv.show_result(res)
        assert rv.fasta_btn.isEnabled() and rv.search_edit.isEnabled()

    def test_missing_library_file_explained(
            self, qapp, golden_dir, golden_store, tmp_path):
        from tescope.gui.widgets.results import ResultsView

        cfg = _cfg(golden_dir, tmp_path / "gone.fa")
        res = run_analysis(cfg, store=golden_store)
        assert res.sequence.status == "unavailable"
        assert any("not found" in w for w in res.warnings)
        rv = ResultsView()
        rv.show_result(res)
        assert not rv.copy_seq_btn.isEnabled()
        assert "not found" in rv.copy_seq_btn.toolTip()


class TestSequencePanel:
    def test_config_roundtrip(self, qapp):
        from tescope.gui.widgets.panels import SequencePanel

        p = SequencePanel()
        cfg = SequenceConfig(enabled=True, source="fasta", path="/x/l.fa")
        p.apply_sequence_config(cfg)
        back = p.sequence_config()
        assert back.enabled and back.source == "fasta" and back.path == "/x/l.fa"
        p.set_source_key("dfam")
        assert not p.browse_btn.isEnabled()
        p.set_source_key("library")
        assert p.browse_btn.isEnabled()


class TestScientificInvariance:
    """Zoom/pan/sequence/search must not alter ANY canonical result."""

    # analysis.yaml and qc.tsv are provenance carriers: they are SUPPOSED
    # to record the sequence source. Every scientific table is compared
    # byte-for-byte; QC is additionally checked to differ only in
    # sequence/warning rows.
    _PROVENANCE_FILES = {"analysis.yaml", "qc.tsv"}

    def _canonical(self, res, tmp_path, tag):
        from tescope.export.tables import export_tables

        out = tmp_path / tag
        files = {p.name: p.read_bytes() for p in export_tables(res, out)}
        science = {k: v for k, v in files.items()
                   if k not in self._PROVENANCE_FILES}
        assert science, "expected scientific tables"
        return science, files.get("qc.tsv", b"")

    @staticmethod
    def _qc_science_rows(qc_bytes):
        return [line for line in qc_bytes.decode().splitlines()
                if not line.startswith(("sequence\t", "warnings\t"))]

    def test_sequence_and_interaction_leave_science_untouched(
            self, qapp, golden_dir, golden_store, seq_fasta, tmp_path):
        from tescope.gui.widgets.results import ResultsView

        base = run_analysis(_cfg(golden_dir), store=golden_store)
        with_seq = run_analysis(_cfg(golden_dir, seq_fasta),
                                store=golden_store)
        a, qc_a = self._canonical(base, tmp_path, "base")
        b, qc_b = self._canonical(with_seq, tmp_path, "seq")
        assert a == b  # byte-identical canonical scientific tables
        # QC gains sequence provenance rows and nothing else changes
        assert self._qc_science_rows(qc_a) == self._qc_science_rows(qc_b)

        # exercise the full interactive surface, then re-export
        rv = ResultsView()
        rv.show_result(with_seq)
        rv.navigator.go_to("101-140")
        rv.navigator.set_view(0, 1000)
        rv.navigator.zoom_to_bases(400, 450)
        rv.search_edit.setText("ACGT")
        rv._search_sequence()
        rv._copy_sequence()
        rv.seq_check.setChecked(False)
        rv.seq_check.setChecked(True)
        rv.profile_canvas.canvas.draw()
        c, _ = self._canonical(rv.result, tmp_path, "after")
        assert a == c
        # numeric spot check on top of the byte comparison
        pd.testing.assert_frame_equal(
            base.continuous_tracks()[0].profile,
            rv.result.continuous_tracks()[0].profile)
        np.testing.assert_array_equal(base.element_covered,
                                      rv.result.element_covered)
