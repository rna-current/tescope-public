"""Dfam feature provider.

Uses the official Dfam REST API (https://dfam.org/api). Dfam data is
released under CC0 1.0 (public-domain dedication; stated in Dfam release
notes), so retrieving, caching and snapshotting records is permitted.

Mapping RepeatMasker name -> Dfam record:

- query ``/families?name=<rm_name>``;
- accept a result only when its own ``name`` equals the RepeatMasker name
  or one of its aliases equals it (UCSC RepeatMasker runs use Repbase
  names, which Dfam records carry as ``aliases``);
- name equality alone NEVER implies coordinate compatibility: each
  candidate is checked by exact consensus-length match against the
  annotation-derived model length
  (:func:`tescope.core.features.model.verify_length_compatibility`).

Extracted feature records (verbatim, source terminology preserved):

- ``features``: typed records such as ``protein_match`` ("BlastX match to
  TE protein DB") with 1-based ``model_start_pos``/``model_end_pos`` and a
  source ``label``;
- ``coding_seqs``: curated CDS records where Dfam provides them.

Nothing is inferred: families without source features yield an empty,
honestly-labeled result.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

from .model import (
    SOURCE_DATABASE,
    Compatibility,
    ConsensusFeature,
    FeatureProvenance,
    verify_length_compatibility,
)

DFAM_API_BASE = "https://dfam.org/api"
_TIMEOUT = 30


class DfamError(RuntimeError):
    """Network/API failure. Callers must degrade gracefully (no analysis failure)."""


class DfamClient:
    """Thin JSON client for the official Dfam API (tests replace _get_json)."""

    def __init__(self, base_url: str = DFAM_API_BASE, timeout: int = _TIMEOUT):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def _get_json(self, path: str, params: dict | None = None) -> dict:
        url = f"{self.base_url}{path}"
        if params:
            url += "?" + urllib.parse.urlencode(params)
        req = urllib.request.Request(url, headers={"User-Agent": "TEscope"})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, OSError, ValueError) as e:
            raise DfamError(f"Dfam API request failed ({url}): {e}") from e

    def _get_text(self, path: str, params: dict | None = None) -> str:
        url = f"{self.base_url}{path}"
        if params:
            url += "?" + urllib.parse.urlencode(params)
        req = urllib.request.Request(url, headers={"User-Agent": "TEscope"})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                return resp.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError, OSError, ValueError) as e:
            raise DfamError(f"Dfam API request failed ({url}): {e}") from e

    def version(self) -> dict:
        return self._get_json("/version")

    def search_families(self, name: str, limit: int = 20) -> list[dict]:
        data = self._get_json("/families", {"name": name, "limit": limit})
        return list(data.get("results", []))

    def family(self, accession: str) -> dict:
        return self._get_json(f"/families/{urllib.parse.quote(accession)}")

    def family_sequence_fasta(self, accession: str) -> str:
        """Consensus sequence of one family (official FASTA endpoint)."""
        return self._get_text(
            f"/families/{urllib.parse.quote(accession)}/sequence",
            {"format": "fasta"},
        )


def _matches_rm_name(record: dict, rm_name: str) -> str | None:
    """Return how the record matches the RepeatMasker name, or None."""
    if record.get("name") == rm_name:
        return "name"
    for a in record.get("aliases") or []:
        if a.get("alias") == rm_name:
            return f"alias:{a.get('database', '?')}"
    return None


def _extract_features(record: dict) -> list[ConsensusFeature]:
    """Verbatim conversion of Dfam feature records to consensus features.

    Dfam model coordinates are 1-based inclusive; converted to 0-based
    half-open. Labels/types/attributes are passed through unchanged.
    """
    out: list[ConsensusFeature] = []
    for f in record.get("features") or []:
        try:
            start = int(f["model_start_pos"]) - 1
            end = int(f["model_end_pos"])
        except (KeyError, TypeError, ValueError):
            continue  # record without usable coordinates: skip, never guess
        attrs = {
            str(a.get("attribute")): str(a.get("value"))
            for a in (f.get("attributes") or [])
            if a.get("attribute") is not None
        }
        out.append(
            ConsensusFeature(
                name=str(f.get("label") or f.get("type") or "feature"),
                start=start,
                end=end,
                origin=SOURCE_DATABASE,
                type=str(f.get("type") or ""),
                note=str(f.get("description") or ""),
                source_attributes=attrs,
            )
        )
    for c in record.get("coding_seqs") or []:
        try:
            start = int(c["cds_start"]) - 1
            end = int(c["cds_end"])
        except (KeyError, TypeError, ValueError):
            continue
        out.append(
            ConsensusFeature(
                name=str(c.get("product") or "coding sequence"),
                start=start,
                end=end,
                origin=SOURCE_DATABASE,
                type="coding_seq",
                note=str(c.get("description") or ""),
                source_attributes={
                    k: str(v)
                    for k, v in c.items()
                    if k not in ("cds_start", "cds_end", "product", "description")
                    and v is not None
                },
            )
        )
    out.sort(key=lambda x: (x.start, x.end))
    return out


def lookup(
    rm_name: str,
    expected_cons_length: int,
    client: DfamClient | None = None,
) -> dict:
    """Look up Dfam annotation for a RepeatMasker repeat name.

    Returns a plain-dict snapshot (cache- and JSON-friendly):

    ``{"features": [...], "provenance": {...}, "candidates": [...]}``

    Only a VERIFIED candidate contributes features. When name-matching
    candidates exist but none is length-compatible, the result carries an
    INCOMPATIBLE provenance explaining every rejection — features are
    never displayed at possibly wrong coordinates.
    """
    client = client or DfamClient()
    version = client.version()
    db_version = f"Dfam {version.get('dfam_version', '?')}"
    release_date = str(version.get("release_date", ""))
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")

    results = client.search_families(rm_name)
    # The search listing omits `aliases`; a result whose name differs may
    # still alias the RepeatMasker name (UCSC uses Repbase names). Check
    # full records for the first few non-name matches.
    candidates = []
    for r in results:
        if _matches_rm_name(r, rm_name):
            candidates.append(r)
    if not candidates:
        for r in results[:10]:
            try:
                full = client.family(r["accession"])
            except DfamError:
                continue
            if _matches_rm_name(full, rm_name):
                full.setdefault("length", r.get("length"))
                candidates.append(full)

    base = FeatureProvenance(
        source="Dfam",
        database_version=db_version,
        database_release_date=release_date,
        retrieved_at=now,
        api_url=client.base_url,
        source_coordinate_convention="1-based inclusive (model_start_pos/model_end_pos)",
    )

    if not candidates:
        base.compatibility = Compatibility.UNVERIFIED.value
        base.compatibility_detail = (
            f"no Dfam family with name or alias '{rm_name}' ({db_version})"
        )
        return {"features": [], "provenance": base.to_dict(), "candidates": []}

    rejections = []
    for r in candidates:
        status, detail = verify_length_compatibility(
            expected_cons_length, r.get("length")
        )
        if status is Compatibility.VERIFIED:
            record = client.family(r["accession"])
            feats = _extract_features(record)
            prov = base
            prov.accession = str(record.get("accession", r["accession"]))
            prov.record_name = str(record.get("name", ""))
            prov.record_version = str(record.get("version", ""))
            prov.source_consensus_length = int(r.get("length"))
            prov.matched_by = _matches_rm_name(r, rm_name) or "name"
            prov.compatibility = Compatibility.VERIFIED.value
            prov.compatibility_detail = detail
            return {
                "features": [f.to_dict() for f in feats],
                "provenance": prov.to_dict(),
                "candidates": [
                    {"accession": c.get("accession"), "name": c.get("name"),
                     "length": c.get("length")}
                    for c in candidates
                ],
            }
        rejections.append(f"{r.get('accession')} ({r.get('name')}): {detail}")

    base.compatibility = Compatibility.INCOMPATIBLE.value
    base.compatibility_detail = "; ".join(rejections)
    return {
        "features": [],
        "provenance": base.to_dict(),
        "candidates": [
            {"accession": c.get("accession"), "name": c.get("name"),
             "length": c.get("length")}
            for c in candidates
        ],
    }
