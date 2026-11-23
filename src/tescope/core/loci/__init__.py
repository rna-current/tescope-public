"""Locus-level views over an analysis result.

A "locus" here is one analyzed TE copy: one row of the filtered
annotation frame (one RepeatMasker record; ``rm_id`` links fragments
that RepeatMasker attributes to the same insertion, and that linkage
is displayed, never silently merged). Locus views are *derived* from
the data the analysis already computed — the per-copy matrices and
projected intervals — and never trigger re-projection, re-estimation,
or reassignment.
"""

from .estimates import EstimateMatchReport, LocusEstimates, load_estimates
from .selection import SelectionError, export_selection, import_selection
from .subset import SubsetResult, subset_result
from .table import build_locus_table, locus_id_for_row, locus_ids

__all__ = [
    "EstimateMatchReport",
    "LocusEstimates",
    "SelectionError",
    "SubsetResult",
    "build_locus_table",
    "export_selection",
    "import_selection",
    "load_estimates",
    "locus_id_for_row",
    "locus_ids",
    "subset_result",
]
