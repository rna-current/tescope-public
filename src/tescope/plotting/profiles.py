"""Figure construction for analysis results.

Uses the matplotlib object API (``Figure``) without pyplot, so the same
code renders headlessly (Agg) in the CLI and embeds in the Qt GUI canvas.
Display transformations (normalization, shared axes) are applied only when
selected in the config and are recorded there.
"""

from __future__ import annotations

import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure

from ..core.aggregation.aggregate import normalize_values
from ..core.analysis import AnalysisResult, TrackResult
from .plotstyle import PlotStyleConfig
from .style import STYLE
from .transforms import copy_quantiles, smooth


def _style_of(result: AnalysisResult) -> PlotStyleConfig:
    """The analysis' plot style (a default for results without one)."""
    return getattr(result.config, "style", None) or PlotStyleConfig()


def _bin_centers(result: AnalysisResult) -> np.ndarray:
    return (result.bin_edges[:-1] + result.bin_edges[1:]) / 2.0


def series_labels(result: AnalysisResult) -> list[str]:
    """Stable series identifiers in display order (for style editors)."""
    return [label for label, *_ in _continuous_series(result, _style_of(result))]


def _continuous_series(result: AnalysisResult, style: PlotStyleConfig):
    """Yield (label, color, line_style, line_width, alpha, values, sem,
    replicate_curves) to plot.

    Grouped tracks are combined into one series per group; ungrouped tracks
    are one series each. Colors resolve deterministically: explicit series
    style > per-track display color hint > palette by display order.
    """
    cfg = result.config
    stat = cfg.aggregation.statistic
    idx = 0
    grouped_names = set()

    def resolved(label, hint):
        color = style.series_style(label).color or hint or style.series_color(label, idx)
        return (
            color,
            style.series_line_style(label, idx),
            style.series_line_width(label),
            style.series_alpha(label),
        )

    for gname, prof in result.group_profiles.items():
        members = [
            t for t in result.continuous_tracks() if t.config.group == gname
        ]
        grouped_names.update(t.name for t in members)
        reps = [
            (t.name, t.profile[stat].to_numpy(dtype=float)) for t in members
        ] if cfg.plot.show_replicates else []
        hint = next((t.config.color for t in members if t.config.color), "")
        yield (
            gname,
            *resolved(gname, hint),
            prof[stat].to_numpy(dtype=float),
            prof["sem"].to_numpy(dtype=float),
            reps,
            None,  # grouped series: no single per-copy matrix (IQR band n/a)
        )
        idx += 1
    for t in result.continuous_tracks():
        if t.name in grouped_names or t.profile is None or not t.config.visible:
            continue
        yield (
            t.name,
            *resolved(t.name, t.config.color),
            t.profile[stat].to_numpy(dtype=float),
            t.profile["sem"].to_numpy(dtype=float),
            [],
            t.matrix,
        )
        idx += 1


def build_figure(result: AnalysisResult, mode: str = "stacked") -> Figure:
    """Build the main results figure ('stacked' or 'overlay').

    When consensus-feature annotation is present and enabled for display,
    a compact feature lane is prepended above the signal panels, aligned
    to the same consensus axis. It is purely visual: signal panels and
    their values are produced by exactly the same code with or without it
    (the no-annotation path is byte-identical to previous releases).
    """
    if mode == "fragments":
        return build_fragment_comparison(result)
    if mode not in ("stacked", "overlay"):
        raise ValueError("mode must be 'stacked', 'overlay' or 'fragments'")
    cfg = result.config
    style = _style_of(result)
    feats = None
    if (
        getattr(result, "features", None) is not None
        and result.features.n_total > 0
        and cfg.annotations.show_in_plots
    ):
        feats = result.features
    series = list(_continuous_series(result, style))
    interval_tracks = [t for t in result.interval_tracks() if t.config.visible]
    point_tracks = [t for t in result.point_tracks() if t.config.visible]

    # sequence track: drawn ONLY for a verified (model-compatible)
    # sequence with display enabled — availability never changes the
    # signal panels themselves
    seq_res = getattr(result, "sequence", None)
    show_seq = bool(
        seq_res is not None and seq_res.verified
        and getattr(cfg, "sequence", None) is not None
        and cfg.sequence.show_track
    )

    n_cont_panels = 1 if (mode == "overlay" and series) else len(series)
    n_extra = 1 if _unique_multi_pair(result) is not None else 0
    n_panels = (n_cont_panels + len(interval_tracks) + len(point_tracks)
                + n_extra + 1)  # +1 coverage
    height = 1.0 + 1.25 * n_panels

    # rc context: typography/line defaults resolve from the style for the
    # duration of figure construction only (pure rendering concern)
    rc = {**STYLE, **style.rc_overrides()}
    with mpl.rc_context(rc):
        fig_w = style.fig_width_in()
        if feats is None:
            lane_rows, feat_ratio = None, 0.0
        else:
            lane_rows = _feature_lane_rows(feats)
            n_rows = sum(len(rows) for rows, _, _ in lane_rows) or 1
            feat_ratio = 0.22 + 0.22 * n_rows
            height += 1.25 * feat_ratio
        seq_ratio = 0.22 if show_seq else 0.0
        height += 1.25 * seq_ratio
        ratios = (([feat_ratio] if feats is not None else [])
                  + [1.0] * n_panels
                  + ([seq_ratio] if show_seq else []))
        fig = Figure(figsize=(fig_w, style.fig_height_in(height)),
                     constrained_layout=True)
        fig.set_dpi(STYLE["figure.dpi"])
        axes = fig.subplots(
            len(ratios), 1, sharex=True,
            gridspec_kw={"height_ratios": ratios} if len(ratios) > 1 else None,
        )
        axes = list(np.atleast_1d(axes))
        seq_ax = axes.pop() if show_seq else None
        if feats is not None:
            _draw_feature_panel(axes[0], lane_rows, result, style)
            axes = axes[1:]
        _draw_panels(result, cfg, style, mode, series, interval_tracks,
                     point_tracks, axes)
        if seq_ax is not None:
            from .sequence_track import draw_sequence_band, refresh_letters

            draw_sequence_band(seq_ax, seq_res.sequence.sequence,
                               result.cons_length)
            # the x label moves to the bottom-most (sequence) axis
            seq_ax.set_xlabel(axes[-1].get_xlabel())
            axes[-1].set_xlabel("")
            # base letters for the configured view window (zoom-dependent:
            # wide windows keep the color band only) — export path; the
            # GUI navigator re-runs this on every zoom/pan
            refresh_letters(seq_ax)
    return fig


def _unique_multi_pair(result):
    """(unique_profile, multi_profile) means when the analysis contains an
    explicit unique+multi track pair from the same BAM file; else None."""
    by_file: dict = {}
    for t in result.continuous_tracks():
        mc = t.config.mapping_class
        if mc in ("unique", "multi") and t.profile is not None:
            by_file.setdefault(t.config.file, {})[mc] = t
    for entry in by_file.values():
        if "unique" in entry and "multi" in entry:
            stat = result.config.aggregation.statistic
            return (entry["unique"].profile[stat].to_numpy(dtype=float),
                    entry["multi"].profile[stat].to_numpy(dtype=float))
    return None


def _draw_panels(result, cfg, style, mode, series, interval_tracks,
                 point_tracks, axes):
    x = _bin_centers(result)
    panel = 0

    norm = cfg.plot.normalization
    ylabel_stat = cfg.aggregation.statistic + {"zscore": " (z)", "minmax": " (0-1)", "none": ""}[norm]

    def draw_series(ax, label, color, ls, lw, alpha, values, sem, reps,
                    matrix=None):
        v = normalize_values(values, norm)
        # display-only smoothing: raw values are never modified; the
        # canonical profile export is untouched (docs/signal-representation.md)
        if style.smoothing_method != "none":
            v_smooth = smooth(v, style.smoothing_method, style.smoothing_window,
                              style.smoothing_sigma, style.smoothing_polyorder)
            if style.smoothing_show_raw:
                ax.plot(x, v, color=color, linestyle=ls,
                        linewidth=max(lw * 0.5, 0.4), alpha=min(alpha, 0.35))
            v_line = v_smooth
        else:
            v_line = v
        if style.representation == "filled":
            fill_c = style.series_style(label).fill_color or color
            ax.fill_between(x, 0, v_line, color=fill_c,
                            alpha=style.fill_alpha, linewidth=0, zorder=1)
        ax.plot(x, v_line, color=color, linestyle=ls, linewidth=lw,
                alpha=alpha, label=label, zorder=3)
        # bands: SEM = uncertainty of the mean; IQR = variation across copies
        if style.band == "sem":
            if cfg.plot.show_sem and norm == "none" and np.isfinite(sem).any():
                ax.fill_between(x, v_line - sem, v_line + sem, color=color,
                                alpha=style.band_alpha, linewidth=0, zorder=2)
        elif style.band == "iqr" and matrix is not None and norm == "none":
            qs = (10.0, 25.0, 75.0, 90.0) if style.outer_band else (25.0, 75.0)
            qv = copy_quantiles(matrix, qs)
            ax.fill_between(x, qv[25.0], qv[75.0], color=color,
                            alpha=style.band_alpha, linewidth=0, zorder=2)
            if style.outer_band:
                ax.fill_between(x, qv[10.0], qv[90.0], color=color,
                                alpha=style.band_alpha * 0.45, linewidth=0,
                                zorder=1.5)
        # (band == "iqr" with a grouped series degrades to no band: a group
        # has no single per-copy matrix; documented in the style dialog/docs)
        for rname, rvals in reps:
            ax.plot(x, normalize_values(rvals, norm), color=color, alpha=0.4,
                    linewidth=max(lw * 0.6, 0.5), label=f"{rname}")

    signal_axes = []
    if series:
        if mode == "overlay":
            ax = axes[panel]
            for label, color, ls, lw, alpha, values, sem, reps, matrix in series:
                draw_series(ax, label, color, ls, lw, alpha, values, sem, reps,
                            matrix)
            ax.set_ylabel(style.y_label or ylabel_stat)
            if style.show_legend:
                if style.legend_position == "outside right":
                    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
                else:
                    ax.legend(loc=style.legend_position)
            signal_axes.append(ax)
            panel += 1
        else:
            for label, color, ls, lw, alpha, values, sem, reps, matrix in series:
                ax = axes[panel]
                draw_series(ax, label, color, ls, lw, alpha, values, sem, reps,
                            matrix)
                ax.set_ylabel(label, rotation=0, ha="right", va="center")
                signal_axes.append(ax)
                panel += 1
            if cfg.plot.y_axis == "shared":
                lims = [axes[i].get_ylim() for i in range(panel)]
                lo = min(lim[0] for lim in lims)
                hi = max(lim[1] for lim in lims)
                for i in range(panel):
                    axes[i].set_ylim(lo, hi)

    for t in interval_tracks:
        ax = axes[panel]
        _draw_interval_panel(ax, t, result, x)
        panel += 1

    for t in point_tracks:
        ax = axes[panel]
        occ = t.point_occurrence
        if occ is not None:
            widths = result.bin_edges[1:] - result.bin_edges[:-1]
            ax.bar(x, occ["point_count"], width=widths * 0.9, color="#666666")
        ax.set_ylabel(f"{t.name}\n(points)", rotation=0, ha="right", va="center")
        panel += 1

    # mapping-support panel: only when the SAME BAM was analyzed as an
    # explicit unique + multi pair (never fabricated for other inputs)
    pair = _unique_multi_pair(result)
    if pair is not None:
        u_prof, m_prof = pair
        ax = axes[panel]
        with np.errstate(invalid="ignore", divide="ignore"):
            tot = u_prof + m_prof
            frac = np.where(tot > 0, u_prof / tot, np.nan)
        ax.plot(x, frac, drawstyle="steps-mid", color="#444444", linewidth=1.0)
        ax.fill_between(x, 0, np.where(np.isfinite(frac), frac, 0),
                        step="mid", color="#444444", alpha=0.15, linewidth=0)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("fraction unique\n(of classified\ncoverage)",
                      rotation=0, ha="right", va="center")
        panel += 1

    # element coverage panel
    ax = axes[panel]
    n_elements = result.element_covered.sum(axis=0)
    ax.fill_between(x, n_elements, step="mid", color="#bbbbbb")
    ax.plot(x, n_elements, drawstyle="steps-mid", color="#555555", linewidth=0.8)
    ax.set_ylabel("TE copies\nper bin", rotation=0, ha="right", va="center")
    ax.set_ylim(bottom=0)

    # feature-region background shading (display only; coordinates come
    # verbatim from the FeatureSet — never invented)
    if style.annotation_shading and getattr(result, "features", None) is not None:
        lane_colors = {"database": style.database_lane_color or "#C7DAE3",
                       "user-defined": style.user_lane_color or "#EBDDC3"}
        for f in result.features.all_features():
            c = lane_colors.get(f.origin, "#C7DAE3")
            for ax_ in signal_axes:
                ax_.axvspan(f.start, f.end, color=c,
                            alpha=style.annotation_shading_alpha, zorder=0)

    # display-only view windows and grid (data are never filtered here)
    if style.y_min is not None or style.y_max is not None:
        for ax_ in signal_axes:
            ax_.set_ylim(bottom=style.y_min, top=style.y_max)
    if style.show_grid:
        for ax_ in axes:
            ax_.grid(True, alpha=style.grid_alpha, linewidth=0.5)
            ax_.set_axisbelow(True)

    axes[-1].set_xlabel(
        style.x_label
        or f"{result.repeat_name} consensus position (bp, 5'→3')"
    )
    x_lo = style.x_min if style.x_min is not None else 0
    x_hi = style.x_max if style.x_max is not None else result.cons_length
    axes[-1].set_xlim(x_lo, x_hi)
    fig = axes[-1].figure
    fig.suptitle(
        style.title
        or f"{result.repeat_name} ({result.config.assembly}) — "
           f"{result.n_selected}/{result.n_annotated} copies",
        fontsize=style.font("title"),
    )
    if style.show_legend and style.legend_position == "outside right":
        # signal to export: use tight bbox so the outside legend is kept
        fig._tescope_tight_bbox = True


def _feature_lane_rows(fs):
    """Greedy row layout per origin lane; returns [(rows, lane_label), ...].

    rows is a list of rows; each row is a list of ConsensusFeature that do
    not overlap within the row (small gap enforced for label readability).
    Database and user-defined features always occupy separate lanes so an
    origin is never conveyed by color alone.
    """
    gap = max(1, fs.cons_length // 100)
    lanes = []
    for feat_list, label, kind in (
        (fs.database_features,
         (fs.provenance.source if fs.provenance else "database"), "database"),
        (fs.user_features, "user-defined", "user"),
    ):
        if not feat_list:
            continue
        rows: list[list] = []
        for f in sorted(feat_list, key=lambda z: (z.start, z.end)):
            for row in rows:
                if f.start >= row[-1].end + gap:
                    row.append(f)
                    break
            else:
                rows.append([f])
        lanes.append((rows, label, kind))
    return lanes


def _draw_feature_panel(ax, lane_rows, result: AnalysisResult,
                        style: PlotStyleConfig | None = None) -> None:
    """Consensus-feature lane: verbatim source/user labels, no data drawn.

    Styling never weakens origin distinctions: the user-defined lane keeps
    its hatch and both lanes keep their lane labels whatever colors are
    chosen; feature coordinates come straight from the FeatureSet.
    """
    style = style or PlotStyleConfig()
    # styles keyed by lane KIND (not position): the user lane always keeps
    # its hatch so origin stays visible even if only one lane is present
    # or lane colors are customized
    lane_defs = {
        "database": {"facecolor": style.database_lane_color or "#C7DAE3",
                     "edgecolor": "#4A7A8C", "hatch": None},
        "user": {"facecolor": style.user_lane_color or "#EBDDC3",
                 "edgecolor": "#9A7B4F", "hatch": "//"},
    }
    y = 0.0
    yticks, ylabels = [], []
    row_h, row_gap, lane_gap = 1.0, 0.25, 0.7
    for rows, label, kind in lane_rows:
        st = lane_defs[kind]
        lane_top = y
        for row in rows:
            for f in row:
                ax.broken_barh(
                    [(f.start, f.end - f.start)], (y - row_h, row_h),
                    facecolor=st["facecolor"], edgecolor=st["edgecolor"],
                    hatch=st["hatch"], alpha=style.annotation_alpha,
                    linewidth=style.annotation_border_width,
                )
                ax.annotate(
                    f.name, ((f.start + f.end) / 2, y - row_h / 2),
                    ha="center", va="center", fontsize=style.font("annotation"),
                    clip_on=True, annotation_clip=True,
                )
            y -= row_h + row_gap
        yticks.append((lane_top + (y + row_gap)) / 2)
        ylabels.append(label)
        y -= lane_gap
    ax.set_ylim(y + lane_gap, 0.35)
    ax.set_yticks(yticks, ylabels, fontsize=max(style.font("annotation"), 7))
    ax.tick_params(axis="y", length=0)
    for side in ("left",):
        ax.spines[side].set_visible(False)
    ax.set_ylabel("")


def _draw_interval_panel(ax, t: TrackResult, result: AnalysisResult, x: np.ndarray) -> None:
    occ = t.occurrence
    if occ is not None:
        ax.plot(
            x,
            occ["fraction_copies_with_peak"],
            drawstyle="steps-mid",
            color="#D55E00",
        )
        ax.set_ylim(bottom=0)
    ax.set_ylabel(
        f"{t.name}\n(fraction of copies)", rotation=0, ha="right", va="center"
    )
    # projected interval segments as a rug band at the top
    if t.projected is not None and len(t.projected):
        top = ax.get_ylim()[1] or 1.0
        band_y = top * 1.05
        for row in t.projected.itertuples():
            ax.plot(
                [row.cons_start, row.cons_end],
                [band_y, band_y],
                color="#8B4513",
                linewidth=3,
                solid_capstyle="butt",
                alpha=0.7,
            )
        ax.set_ylim(top=band_y * 1.1)


FRAGMENT_SUBSETS = (
    ("all copies", None),
    ("full-length", "full_length"),
    ("5'-truncated", "five_prime_truncated"),
)


def build_fragment_comparison(result: AnalysisResult) -> Figure:
    """All vs full-length vs 5'-truncated copies on one consensus axis.

    Display-level regrouping of the EXISTING per-copy projections using the
    EXISTING completeness classification attached to the analyzed copies:
    a subset's profile equals what a completeness-filtered analysis would
    produce (projection is per-copy independent). No new categories, no new
    statistics — subset mean/median of already-computed per-copy values,
    with per-subset copy-support panels.
    """
    cfg = result.config
    style = _style_of(result)
    stat = cfg.aggregation.statistic
    if "completeness" not in result.copies.columns:
        raise ValueError("fragment comparison requires completeness classes")
    comp = result.copies["completeness"].astype(str).to_numpy()
    x = _bin_centers(result)
    series = [srs for srs in _continuous_series(result, style)]

    def subset_profile(matrix, members, mask):
        import warnings as _w

        with _w.catch_warnings():
            _w.simplefilter("ignore", RuntimeWarning)
            if matrix is not None:
                vals = matrix.values[mask].astype("float64")
                return (np.nanmedian(vals, axis=0) if stat == "median"
                        else np.nanmean(vals, axis=0))
            # grouped series: mean across member tracks' subset profiles
            profs = []
            for t in members:
                vals = t.matrix.values[mask].astype("float64")
                profs.append(np.nanmedian(vals, axis=0) if stat == "median"
                             else np.nanmean(vals, axis=0))
            return np.nanmean(np.stack(profs), axis=0)

    rc = {**STYLE, **style.rc_overrides()}
    n_sub = len(FRAGMENT_SUBSETS)
    with mpl.rc_context(rc):
        fig_w = style.fig_width_in()
        height = 1.0 + 1.6 * n_sub
        fig = Figure(figsize=(fig_w, style.fig_height_in(height)),
                     constrained_layout=True)
        fig.set_dpi(STYLE["figure.dpi"])
        ratios = []
        for _ in range(n_sub):
            ratios += [1.0, 0.35]
        axes = fig.subplots(2 * n_sub, 1, sharex=True,
                            gridspec_kw={"height_ratios": ratios})
        axes = list(np.atleast_1d(axes))
        members_of = {
            g: [t for t in result.continuous_tracks() if t.config.group == g]
            for g in result.group_profiles
        }
        for i, (label, cls) in enumerate(FRAGMENT_SUBSETS):
            mask = np.ones(len(comp), bool) if cls is None else comp == cls
            ax, ax_sup = axes[2 * i], axes[2 * i + 1]
            for sname, color, ls, lw, alpha, _v, _s, _r, matrix in series:
                if matrix is None and not members_of.get(sname):
                    continue
                prof = subset_profile(matrix, members_of.get(sname, []), mask)
                if style.representation == "filled":
                    ax.fill_between(x, 0, prof, color=color,
                                    alpha=style.fill_alpha, linewidth=0)
                ax.plot(x, prof, color=color, linestyle=ls, linewidth=lw,
                        alpha=alpha, label=sname)
            ax.set_ylabel(f"{label}\n(n={int(mask.sum()):,})", rotation=0,
                          ha="right", va="center")
            if i == 0 and style.show_legend and len(series) > 1:
                ax.legend(fontsize=style.font("legend"))
            support = result.element_covered[mask].sum(axis=0)
            ax_sup.fill_between(x, support, step="mid", color="#bbbbbb")
            ax_sup.plot(x, support, drawstyle="steps-mid", color="#555555",
                        linewidth=0.7)
            ax_sup.set_ylabel("copies\nper bin", rotation=0, ha="right",
                              va="center", fontsize=style.font("annotation"))
            ax_sup.set_ylim(bottom=0)
        axes[-1].set_xlabel(
            style.x_label
            or f"{result.repeat_name} consensus position (bp, 5'→3')")
        axes[-1].set_xlim(0, result.cons_length)
        fig.suptitle(
            style.title
            or f"{result.repeat_name} ({cfg.assembly}) — copies by "
               "completeness class",
            fontsize=style.font("title"))
    return fig


def build_subset_comparison(result: AnalysisResult, subset,
                            compare: bool = True) -> Figure:
    """Explicitly applied locus subset vs the eligible population.

    ``subset`` is a :class:`tescope.core.loci.subset.SubsetResult` —
    aggregates recomputed from the ORIGINAL per-copy matrices by the
    same :func:`aggregate_profile` as the main analysis. The eligible
    population's aggregate is drawn unchanged for context (``compare``)
    or omitted (subset-only view); both denominators are always in the
    labels. SEM bands are uncertainty of the mean across the loci of
    each population — copy-to-copy variation, NOT replicate variation,
    and not a confidence interval.
    """
    style = _style_of(result)
    stat = result.config.aggregation.statistic
    x = _bin_centers(result)
    tracks = [t for t in result.continuous_tracks()
              if t.name in subset.profiles]
    if not tracks:
        raise ValueError("subset comparison requires a continuous track")
    rc = {**STYLE, **style.rc_overrides()}
    with mpl.rc_context(rc):
        fig_w = style.fig_width_in()
        fig = Figure(figsize=(fig_w, style.fig_height_in(1.0 + 1.9 * len(tracks))),
                     constrained_layout=True)
        fig.set_dpi(STYLE["figure.dpi"])
        ratios: list[float] = []
        for _ in tracks:
            ratios += [1.0, 0.35]
        axes = list(np.atleast_1d(fig.subplots(
            2 * len(tracks), 1, sharex=True,
            gridspec_kw={"height_ratios": ratios})))
        elig_label = f"eligible (N={subset.n_eligible:,})"
        sub_label = f"subset (n={subset.n_subset:,})"
        for i, t in enumerate(tracks):
            ax, ax_sup = axes[2 * i], axes[2 * i + 1]
            color = style.series_color(t.name, i)
            sub_prof = subset.profiles[t.name]
            if compare and t.profile is not None:
                ep = t.profile
                ax.plot(x, ep[stat].to_numpy(float), color="#8a8a8a",
                        linewidth=style.resolved_line_width(),
                        label=elig_label)
                if style.band == "sem":
                    ax.fill_between(
                        x,
                        ep[stat].to_numpy(float) - ep["sem"].to_numpy(float),
                        ep[stat].to_numpy(float) + ep["sem"].to_numpy(float),
                        color="#8a8a8a", alpha=0.18, linewidth=0)
            sp = sub_prof[stat].to_numpy(float)
            if style.representation == "filled" and not compare:
                ax.fill_between(x, 0, sp, color=color,
                                alpha=style.fill_alpha, linewidth=0)
            ax.plot(x, sp, color=color,
                    linewidth=style.resolved_line_width(), label=sub_label)
            if style.band == "sem":
                ax.fill_between(
                    x,
                    sp - sub_prof["sem"].to_numpy(float),
                    sp + sub_prof["sem"].to_numpy(float),
                    color=color, alpha=0.25, linewidth=0)
            ax.set_ylabel(t.name, fontsize=style.font("label"))
            if i == 0 and style.show_legend:
                ax.legend(fontsize=style.font("legend"))
            ax_sup.plot(x, sub_prof["n_elements"].to_numpy(),
                        drawstyle="steps-mid", color=color, linewidth=0.9,
                        label=sub_label)
            if compare and t.profile is not None:
                ax_sup.plot(x, t.profile["n_elements"].to_numpy(),
                            drawstyle="steps-mid", color="#8a8a8a",
                            linewidth=0.7, label=elig_label)
            ax_sup.set_ylabel("copies\nper bin", rotation=0, ha="right",
                              va="center", fontsize=style.font("annotation"))
            ax_sup.set_ylim(bottom=0)
        axes[-1].set_xlabel(
            style.x_label
            or f"{result.repeat_name} consensus position (bp, 5'→3')")
        axes[-1].set_xlim(0, result.cons_length)
        fig.suptitle(
            style.title
            or (f"{result.repeat_name} — explicit locus subset "
                + ("vs eligible population" if compare else "(subset only)")),
            fontsize=style.font("title"))
    return fig


def heatmap_order(result: AnalysisResult, track: TrackResult,
                  sort: str) -> np.ndarray:
    """Public row-order accessor (row i of the drawn heatmap =
    ``result.copies`` positional row ``heatmap_order(...)[i]``)."""
    return _heatmap_order(result, track, sort)


def _heatmap_order(result: AnalysisResult, track: TrackResult,
                   sort: str) -> np.ndarray:
    """Row ORDER for the heatmap (values are never changed)."""
    if sort == "genome_length":
        key = (result.copies["genome_end"] - result.copies["genome_start"]).to_numpy()
    elif sort == "consensus_coverage" and "consensus_coverage" in result.copies:
        key = result.copies["consensus_coverage"].to_numpy()
    elif sort == "total_signal":
        import warnings as _w

        with _w.catch_warnings():
            _w.simplefilter("ignore", RuntimeWarning)
            key = np.nansum(track.matrix.values.astype("float64"), axis=1)
    else:
        key = result.copies["cons_start"].to_numpy()
    return np.argsort(key, kind="stable")


def build_heatmap(result: AnalysisResult, track: TrackResult,
                  highlight_row: int | None = None) -> Figure:
    """Copies x bins heatmap for one continuous track (copies sorted by cons_start).

    Colormap and vmin/vmax are DISPLAY settings only — the underlying
    matrix is never rescaled or modified (display color limits are not
    normalization). ``highlight_row`` (a positional row of
    ``result.copies``) draws a display-only marker on that copy's
    heatmap row; it changes no value and no row order.
    """
    if track.matrix is None:
        raise ValueError(f"track {track.name} has no continuous matrix")
    style = _style_of(result)
    m = track.matrix
    order = _heatmap_order(result, track, style.heatmap_sort)
    data = m.values[order]
    with mpl.rc_context({**STYLE, **style.rc_overrides()}):
        w = style.fig_width_in()
        if style.heatmap_show_profile and track.profile is not None:
            fig = Figure(figsize=(w, style.fig_height_in(6.0)),
                         constrained_layout=True)
            ax_prof, ax = fig.subplots(
                2, 1, sharex=True, gridspec_kw={"height_ratios": [0.3, 1.0]})
            stat = result.config.aggregation.statistic
            x = _bin_centers(result)
            prof = track.profile[stat].to_numpy(dtype=float)
            color = style.series_color(track.name, 0)
            if style.representation == "filled":
                ax_prof.fill_between(x, 0, prof, color=color,
                                     alpha=style.fill_alpha, linewidth=0)
            ax_prof.plot(x, prof, color=color,
                         linewidth=style.resolved_line_width())
            ax_prof.set_ylabel(stat, fontsize=style.font("annotation"))
        else:
            fig = Figure(figsize=(w, style.fig_height_in(4.8)),
                         constrained_layout=True)
            ax = fig.subplots(1, 1)
        im = ax.imshow(
            data,
            aspect="auto",
            interpolation="nearest",
            extent=(0, result.cons_length, len(data), 0),
            cmap=style.heatmap_cmap,
            vmin=style.heatmap_vmin,
            vmax=style.heatmap_vmax,
        )
        if highlight_row is not None:
            at = np.where(order == highlight_row)[0]
            if at.size:
                y = float(at[0])
                ax.axhspan(y, y + 1, facecolor="none",
                           edgecolor="#B22222", linewidth=1.2, zorder=5)
        if style.heatmap_colorbar:
            fig.colorbar(im, ax=ax, label=track.name, shrink=0.8)
        ax.set_xlabel(style.x_label
                      or f"{result.repeat_name} consensus position (bp)")
        ax.set_ylabel(f"TE copies (n={len(data)})")
        ax.set_title(style.title
                     or f"{track.name} — per-copy signal "
                        f"(rows sorted by {style.heatmap_sort})")
    return fig
