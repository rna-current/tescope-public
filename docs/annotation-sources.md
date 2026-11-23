# Consensus-annotation sources: feasibility audit & decision record

Audited 2026-08-31 against live official interfaces (not website
appearance). Conservative policy: integrate only sources that are
scientifically defensible, reproducible, legally clear, and testable.

## Dfam — IMPLEMENTED

| item | finding |
|---|---|
| Official URL | https://dfam.org (API: https://dfam.org/api; releases: https://dfam.org/releases) |
| Access mechanism | official REST API (JSON), verified live: `/version`, `/families?name=`, `/families/{accession}`; full downloadable releases (FamDB) also exist |
| License | **CC0 1.0** public-domain dedication (stated in Dfam release notes) — caching/snapshotting permitted |
| Identifiers | `DF…` accessions, family names, aliases incl. **Repbase names** (the names used by UCSC RepeatMasker runs) |
| Consensus availability | full `consensus_sequence` + `length` per record |
| Feature information | `features[]`: typed records (observed: `protein_match`, "BlastX match to TE protein DB") with `model_start_pos`/`model_end_pos`, source `label`, attributes (evalue, orientation, protein coords); `coding_seqs[]`: curated CDS where present |
| Coordinate convention | 1-based inclusive model coordinates (converted; original preserved in provenance) |
| Version info | `/version`: dfam_version + release date (recorded per lookup) |
| RM-family mapping | name equality or Repbase-alias equality (search listing omits aliases → full records of top candidates are checked) |
| Coordinate verification | exact consensus-length match vs TEscope's annotation-derived model length; mismatch ⇒ features withheld |
| Caveat found | Dfam 4.0 splits some classic families (L1HS → `_5end`/`_3end`; likewise mouse L1Md*) — such families are honestly `incompatible` with monolithic UCSC models. Verified live positive case: `IAPEz-int` (DF000004146, len 6481 == mm39/mm10 model) with 3 protein_match features |
| Decision | **implement** (first and only automatic provider in this milestone) |

## Repbase (GIRI) — REJECTED for automatic integration

Subscription/licensed; no open machine-readable API; redistribution and
programmatic use not clearly permitted. TEscope touches Repbase only as
*names* (Dfam's alias metadata), never content.

## Protein/domain databases (UniProt, Pfam/InterPro) — REJECTED for now

They describe protein space. Mapping protein-domain coordinates onto a TE
DNA consensus requires an alignment/inference step — exactly the kind of
biological interpolation this subsystem forbids. Dfam's `protein_match`
features already deliver consensus-space homology evidence with
provenance, making this integration unnecessary.

## UCSC Repeat Browser — REJECTED

Distributes consensus "hub" genomes for browser use; no stable
per-feature API; itself derived from Dfam/Repbase; integration would mean
fragile scraping.

## L1Base and similar curated element sets — DEFERRED

Curated *genomic element* annotations (not consensus-feature records),
no machine API; small scope. Revisit only if a manuscript figure needs it.

---

# Implementation decision record (2026-08-31 milestone)

## IMPLEMENTED
- Provider abstraction with one real provider (Dfam) + user-defined
  features; local cache; per-analysis `features_used.json` snapshot.
- Consensus compatibility states verified/incompatible/unverified with
  stored reasons; incompatible/unverified ⇒ features withheld.
- Full provenance (database version, accession, match route, retrieval
  time, source convention, per-feature attributes).
- Config section (`annotations:`), GUI panel (enable, availability
  status, threaded Dfam check, user-feature table + TSV/YAML import,
  show-in-plots), CLI (`tescope annotation show`, `profile --annotate
  --user-features`), figure lane (origin-separated, overlap-stacked),
  QC (status, provenance, per-feature copy coverage), exports.
- Feature-coverage QC (overlap / full-cover counts over analyzed copies;
  exact definitions in docs/annotations.md).

## DEFERRED (plot-style milestone 2026-09-02 additions)
- Per-figure style overrides (one style per analysis is enough for the
  manuscript workflow); style themes beyond the three presets; custom
  dash patterns; TIFF compression options; per-panel grid control.

## DEFERRED
- Annotation-aware copy filtering ("copies overlapping feature X"):
  semantics are clear, but it would couple annotation retrieval into the
  filtering stage that this milestone deliberately keeps annotation-free;
  revisit as its own reviewed change if manuscript work needs it.
- Additional providers (see above), sequence-checksum compatibility
  (would need the RM run's consensus sequences, which UCSC `.out` does
  not provide), Dfam FamDB offline bulk import.

## REJECTED
- Any Claude/AI-generated feature table (hard rule; nothing of the kind
  exists in the repository — test fixtures use SYNTH_* names).
- Automatic biological interpretation of peaks; ORF/domain prediction;
  cross-family equivalence; multi-database conflict resolution.
