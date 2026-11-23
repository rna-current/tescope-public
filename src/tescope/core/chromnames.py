"""Chromosome naming reconciliation (chr1 vs 1, chrM vs MT).

Policy: never silently remap. :func:`reconcile` computes an explicit
mapping from *track* chromosome names to *reference* names using only
unambiguous, well-known aliases, and reports everything it did (and could
not do) in a :class:`ChromReport` that surfaces in QC output.
"""

from __future__ import annotations

from dataclasses import dataclass, field


def _alias_candidates(name: str) -> list[str]:
    """Unambiguous alias spellings for a chromosome name (besides itself)."""
    out = []
    if name.startswith("chr"):
        stripped = name[3:]
        out.append(stripped)
        if stripped == "M":
            out.append("MT")
    else:
        out.append("chr" + name)
        if name == "MT":
            out.extend(["chrM", "M"])
        if name == "M":
            out.append("chrM")
    return out


@dataclass
class ChromReport:
    """Result of reconciling one track against the reference naming."""

    matched: dict[str, str] = field(default_factory=dict)  # track name -> reference name
    aliased: dict[str, str] = field(default_factory=dict)  # subset of matched via alias
    unmatched_track: list[str] = field(default_factory=list)
    unmatched_reference: list[str] = field(default_factory=list)

    @property
    def n_aliased(self) -> int:
        return len(self.aliased)

    def translate(self, ref_name: str) -> str | None:
        """Reference chromosome name -> the track's name for it (or None)."""
        return self._ref_to_track.get(ref_name)

    def __post_init__(self):
        self._ref_to_track: dict[str, str] = {}

    def _finalize(self):
        self._ref_to_track = {ref: trk for trk, ref in self.matched.items()}

    def summary(self) -> str:
        parts = [f"{len(self.matched)} chromosomes matched"]
        if self.aliased:
            ex = next(iter(self.aliased.items()))
            parts.append(f"{len(self.aliased)} via alias (e.g. {ex[0]} -> {ex[1]})")
        if self.unmatched_track:
            parts.append(
                f"{len(self.unmatched_track)} track chromosomes not in reference "
                f"(e.g. {self.unmatched_track[0]})"
            )
        return "; ".join(parts)


class ChromMismatchError(RuntimeError):
    pass


def reconcile(
    track_chroms: list[str] | dict,
    reference_chroms: list[str] | dict,
    *,
    allow_alias: bool = True,
) -> ChromReport:
    """Map track chromosome names onto reference names.

    Exact matches always win. With ``allow_alias`` (default), the standard
    chr-prefix and mitochondrial aliases are applied when they resolve to
    exactly one reference name. Anything else is reported unmatched, never
    guessed.
    """
    track_names = list(track_chroms)
    ref_names = set(reference_chroms)
    report = ChromReport()
    for t in track_names:
        if t in ref_names:
            report.matched[t] = t
            continue
        if allow_alias:
            hits = [c for c in _alias_candidates(t) if c in ref_names]
            if len(hits) == 1:
                report.matched[t] = hits[0]
                report.aliased[t] = hits[0]
                continue
        report.unmatched_track.append(t)
    used = set(report.matched.values())
    report.unmatched_reference = sorted(ref_names - used)
    report._finalize()
    return report


def require_overlap(report: ChromReport, track_name: str) -> None:
    """Raise a clear error when a track shares no chromosomes with the reference."""
    if not report.matched:
        raise ChromMismatchError(
            f"Track '{track_name}' uses chromosomes not found in the selected "
            f"reference (e.g. {', '.join(report.unmatched_track[:3])}). "
            "Check that the track was produced on the selected assembly."
        )
