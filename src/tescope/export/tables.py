"""Machine-readable result export (TSV/CSV).

Every analysis writes tables sufficient to reconstruct figures
independently (spec: no image-only outputs):

```
<out>/profile_<track>.tsv          per-bin continuous profiles
<out>/group_<group>.tsv            combined replicate profiles
<out>/occurrence_<track>.tsv       interval occurrence per bin
<out>/projected_intervals_<track>.tsv   projected peak coordinates
<out>/points_<track>.tsv           projected point events
<out>/point_occurrence_<track>.tsv
<out>/copies.tsv                   the selected TE copies (post-filter)
<out>/qc.tsv                       QC summary
<out>/analysis.yaml                the exact config (reproducibility)
```
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from ..config.io import save_config
from ..core.analysis import AnalysisError, AnalysisResult
from ..core.qc import qc_table

_SEP = {"tsv": "\t", "csv": ","}


def _slug(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    return s or "track"


def _annotate_profile(df: pd.DataFrame, result: AnalysisResult, signal_name: str,
                      group: str, aggregation: str) -> pd.DataFrame:
    out = df.copy()
    out.insert(0, "repeat_name", result.repeat_name)
    out.insert(1, "signal_name", signal_name)
    out.insert(2, "group", group)
    out.insert(3, "aggregation", aggregation)
    return out


def export_tables(result: AnalysisResult, out_dir: str | Path,
                  formats: list[str] | None = None) -> list[Path]:
    """Write all result tables; returns the list of files written."""
    formats = formats or result.config.output.table_formats or ["tsv"]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    def write(df: pd.DataFrame, stem: str):
        for fmt in formats:
            p = out / f"{stem}.{fmt}"
            df.to_csv(p, sep=_SEP[fmt], index=False)
            written.append(p)

    stat = result.config.aggregation.statistic
    for t in result.continuous_tracks():
        if t.profile is not None:
            write(
                _annotate_profile(t.profile, result, t.name, t.config.group, stat),
                f"profile_{_slug(t.name)}",
            )
    for gname, prof in result.group_profiles.items():
        write(
            _annotate_profile(prof, result, f"group:{gname}", gname, stat),
            f"group_{_slug(gname)}",
        )
    for t in result.interval_tracks():
        if t.occurrence is not None:
            write(
                _annotate_profile(t.occurrence, result, t.name, t.config.group, "occurrence"),
                f"occurrence_{_slug(t.name)}",
            )
        if t.projected is not None:
            proj = t.projected.copy()
            proj.insert(0, "repeat_name", result.repeat_name)
            proj.insert(1, "signal_name", t.name)
            write(proj, f"projected_intervals_{_slug(t.name)}")
    for t in result.point_tracks():
        if t.points is not None:
            pts = t.points.copy()
            pts.insert(0, "repeat_name", result.repeat_name)
            pts.insert(1, "signal_name", t.name)
            write(pts, f"points_{_slug(t.name)}")
        if t.point_occurrence is not None:
            write(
                _annotate_profile(t.point_occurrence, result, t.name, t.config.group, "occurrence"),
                f"point_occurrence_{_slug(t.name)}",
            )

    # -- derived DISPLAY exports (clearly labeled; canonical profile_*.tsv
    # -- files above are never replaced or modified) -----------------------
    style = getattr(result.config, "style", None)
    if style is not None and style.smoothing_method != "none":
        from tescope import __version__ as _ver

        from ..plotting.transforms import smooth as _smooth

        for t in result.continuous_tracks():
            if t.profile is None:
                continue
            sm = t.profile.copy()
            for col in ("mean", "median"):
                sm[col] = _smooth(
                    sm[col].to_numpy(dtype=float), style.smoothing_method,
                    style.smoothing_window, style.smoothing_sigma,
                    style.smoothing_polyorder)
            sm.insert(0, "derived_from", f"profile_{_slug(t.name)}.tsv")
            sm.insert(1, "smoothing_method", style.smoothing_method)
            sm.insert(2, "smoothing_params",
                      f"window={style.smoothing_window};"
                      f"sigma={style.smoothing_sigma};"
                      f"polyorder={style.smoothing_polyorder}")
            sm.insert(3, "tescope_version", _ver)
            write(sm, f"profile_smoothed_{_slug(t.name)}")
    if style is not None and style.band == "iqr":
        import numpy as _np

        from ..plotting.transforms import copy_quantiles as _cq

        for t in result.continuous_tracks():
            if t.matrix is None or t.profile is None:
                continue
            qv = _cq(t.matrix, (10.0, 25.0, 50.0, 75.0, 90.0))
            qt = t.profile[["consensus_bin", "consensus_start",
                            "consensus_end", "n_observations"]].copy()
            for q in (10, 25, 50, 75, 90):
                qt[f"p{q}"] = qv[float(q)]
            qt.insert(0, "summary", "per-copy distribution (variation across TE copies)")
            write(qt, f"profile_quantiles_{_slug(t.name)}")
    # mapping support: explicit unique+multi pair from the same BAM
    from ..plotting.profiles import _unique_multi_pair

    pair = _unique_multi_pair(result)
    if pair is not None:
        import numpy as _np

        u, m = pair
        base = result.continuous_tracks()[0].profile
        with _np.errstate(invalid="ignore", divide="ignore"):
            tot = u + m
            frac = _np.where(tot > 0, u / tot, _np.nan)
        ms = pd.DataFrame({
            "consensus_bin": base["consensus_bin"],
            "consensus_start": base["consensus_start"],
            "consensus_end": base["consensus_end"],
            "unique_coverage": u,
            "multi_coverage": m,
            "fraction_unique": frac,
        })
        write(ms, "mapping_support")

    if result.features is not None and result.features.n_total:
        feats = result.features
        rows = [
            {
                "repeat_name": result.repeat_name,
                "feature_name": f.name,
                "consensus_start": f.start,
                "consensus_end": f.end,
                "origin": f.origin,
                "type": f.type,
                "note": f.note,
            }
            for f in feats.all_features()
        ]
        write(pd.DataFrame(rows), "features")
    if result.features is not None:
        # reproducibility snapshot: the exact annotation records + provenance
        import json as _json

        fp = out / "features_used.json"
        fp.write_text(_json.dumps(result.features.to_dict(), indent=2),
                      encoding="utf-8")
        written.append(fp)

    copies = result.copies.copy()
    copies.insert(0, "repeat_name_", copies.pop("repeat_name"))
    copies = copies.rename(columns={"repeat_name_": "repeat_name"})
    write(copies, "copies")

    # per-locus table: stable locus_id + copy_idx (joins projected_* /
    # points_* exports) + per-track per-locus summaries (+ est_* when
    # an estimate table is attached)
    loci = result.locus_table().copy()
    loci.insert(1, "copy_idx", loci.index.to_numpy())
    write(loci, "loci")

    for gset in result.genome_features:
        grows = [
            {
                "assembly": gset.assembly,
                "chrom": f.chrom,
                "start": f.start,
                "end": f.end,
                "strand": f.strand,
                "feature_id": f.feature_id,
                "label": f.label,
                "type": f.type,
                "parent_id": f.parent_id,
                "source_file": gset.path,
                "coordinate_source": gset.format,
            }
            for f in gset.features
        ]
        write(pd.DataFrame(grows), f"genome_features_{_slug(Path(gset.path).stem)}")

    if result.config.loci.subset_applied and result.config.loci.selected:
        # explicit subset analysis: recomputed from the ORIGINAL per-locus
        # values (same derived-export contract as profile_smoothed_*)
        from ..core.loci.subset import SubsetError, subset_result

        try:
            sub = subset_result(result, list(result.config.loci.selected))
        except SubsetError as e:
            raise AnalysisError(f"Subset export failed: {e}") from e
        for tname, prof in sub.profiles.items():
            df = prof.copy()
            df.insert(0, "population", sub.label)
            df.insert(0, "signal_name", tname)
            df.insert(0, "repeat_name", result.repeat_name)
            write(df, f"subset_profile_{_slug(tname)}")
        for tname, occ in sub.occurrences.items():
            df = occ.copy()
            df.insert(0, "population", sub.label)
            df.insert(0, "signal_name", tname)
            df.insert(0, "repeat_name", result.repeat_name)
            write(df, f"subset_occurrence_{_slug(tname)}")
        if sub.invalid_ids:
            rp = out / "subset_invalid_ids.txt"
            rp.write_text("\n".join(sub.invalid_ids) + "\n", encoding="utf-8")
            written.append(rp)

    write(qc_table(result), "qc")

    cfg_path = out / "analysis.yaml"
    save_config(result.config, cfg_path)
    written.append(cfg_path)
    return written
