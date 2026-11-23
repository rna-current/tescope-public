"""Feature resolution: config -> FeatureSet, cache-first, failure-tolerant.

Kept strictly separate from parsing/projection/aggregation: annotation
failure degrades to "no annotation" with a status message and never fails
the scientific analysis unless the config marks annotation as required.
"""

from __future__ import annotations

from .cache import load_cached, store_cached
from .dfam import DfamClient, DfamError, lookup
from .model import (
    Compatibility,
    ConsensusFeature,
    FeatureProvenance,
    FeatureSet,
)
from .user import UserFeatureError, feature_from_record, parse_user_features_file


class AnnotationRequiredError(RuntimeError):
    """Raised only when the config demands annotation and none is available."""


def resolve_features(
    ann_cfg,
    repeat_name: str,
    cons_length: int,
    cache_dir=None,
    client: DfamClient | None = None,
    allow_network: bool = True,
) -> FeatureSet:
    """Build the FeatureSet for an analysis.

    ``ann_cfg`` is :class:`tescope.config.model.AnnotationConfig`.
    Order: user features (validated against the consensus length) always;
    database features from cache first, then network (when allowed).
    """
    fs = FeatureSet(repeat_name=repeat_name, cons_length=cons_length)
    statuses: list[str] = []

    # -- user-defined features --------------------------------------------
    user_feats: list[ConsensusFeature] = []
    for i, rec in enumerate(ann_cfg.user_features or []):
        user_feats.append(feature_from_record(dict(rec), where=f"user_features[{i}]"))
    if ann_cfg.user_features_file:
        user_feats.extend(parse_user_features_file(ann_cfg.user_features_file))
    for f in user_feats:
        f.validate(cons_length)  # bounds against this consensus
    fs.user_features = user_feats
    if user_feats:
        statuses.append(f"{len(user_feats)} user-defined feature(s)")

    # -- database features -------------------------------------------------
    if ann_cfg.provider == "dfam":
        snapshot = load_cached(repeat_name, cons_length, "dfam", cache_dir)
        from_cache = snapshot is not None
        if snapshot is None and allow_network:
            try:
                snapshot = lookup(repeat_name, cons_length, client=client)
                store_cached(repeat_name, cons_length, snapshot, "dfam", cache_dir)
            except DfamError as e:
                statuses.append(f"Dfam unavailable ({e})")
                snapshot = None
        elif snapshot is None:
            statuses.append("Dfam: not cached (network disabled)")

        if snapshot is not None:
            prov = FeatureProvenance.from_dict(snapshot.get("provenance") or {})
            fs.provenance = prov
            if prov.compatibility == Compatibility.VERIFIED.value:
                feats = [
                    ConsensusFeature.from_dict(x)
                    for x in snapshot.get("features", [])
                ]
                for f in feats:
                    f.validate(cons_length)
                fs.database_features = feats
                origin = "cached" if from_cache else "fetched"
                statuses.append(
                    f"Dfam {prov.accession} verified ({origin}), "
                    f"{len(feats)} feature(s)"
                )
            elif prov.compatibility == Compatibility.INCOMPATIBLE.value:
                statuses.append("Dfam: incompatible consensus — features withheld")
            else:
                statuses.append("Dfam: no matching family (unverified)")
    elif ann_cfg.provider not in ("", "none"):
        raise ValueError(f"unknown annotation provider {ann_cfg.provider!r}")

    fs.status = "; ".join(statuses) if statuses else "no annotation"
    fs.validate()

    if ann_cfg.required and fs.n_total == 0:
        raise AnnotationRequiredError(
            f"Annotation was marked required but none is available for "
            f"'{repeat_name}' ({fs.status})."
        )
    return fs


__all__ = [
    "AnnotationRequiredError",
    "DfamClient",
    "DfamError",
    "UserFeatureError",
    "resolve_features",
]
