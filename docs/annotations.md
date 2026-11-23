# Consensus feature annotations

TEscope can display known features of a TE consensus (for example coding
or protein-homology regions) in a lane aligned to the same consensus axis
as the signal profiles.

**The cardinal rule: TEscope does not invent biological feature
annotations.** A feature is displayed only when it comes verbatim from an
authoritative database record (with full provenance) or was explicitly
entered by you. If no reliable annotation exists, TEscope displays no
annotation. Database evidence keeps its source terminology — a Dfam
`protein_match` ("BlastX match to TE protein DB") is never promoted to
"confirmed ORF" or a gene name.

## Automatic source: Dfam

[Dfam](https://dfam.org) (CC0 1.0 public-domain data) is the supported
automatic provider, via its official REST API.

Workflow when annotation is enabled:

```
RepeatMasker repeat name
   → Dfam family search (matched by name or Repbase alias —
     UCSC RepeatMasker annotations use Repbase names)
   → consensus COMPATIBILITY CHECK (below)
   → verified record's features + coding sequences retrieved verbatim
   → cached locally with provenance
   → displayed in the "Dfam" lane
```

Feature types actually obtainable from Dfam records: `protein_match`
evidence (with e-value, orientation and protein coordinates preserved as
attributes) and curated `coding_seq` records, where the database provides
them. Many families have none — that is shown honestly as "no database
features".

## Consensus compatibility check

A name match never proves coordinate compatibility. Before any database
feature is placed on the analysis axis, TEscope compares the source
consensus length with the annotation-derived model length (the same
model length used everywhere else in the analysis):

| state | meaning | features shown? |
|---|---|---|
| `verified` | exact consensus-length match | yes |
| `incompatible` | lengths differ — coordinates cannot be trusted | **no** (withheld, with a per-candidate explanation) |
| `unverified` | no matching family / missing data | no |

Real example of the guard working: recent Dfam releases split some
classic families (e.g. `L1HS` into `L1HS_5end`/`L1HS_3end`), so a
monolithic RepeatMasker `L1HS` model has no length-compatible Dfam record
— TEscope reports *incompatible* and displays nothing, rather than
placing features at wrong positions.

## Provenance

Every database annotation retains: source database and release version
(e.g. "Dfam 4.0", release date), accession and record name/version, how
it was matched (name vs alias), retrieval timestamp, API URL, the
source's coordinate convention, the source consensus length, the
compatibility state and the reason for it, and each feature's original
type/label/attributes. This appears in the QC table and in
`features_used.json`.

## Caching, offline use, reproducible snapshots

Lookups are cached under `<cache>/features/dfam/` keyed by repeat name
and model length. After the first lookup everything works offline;
network or database failure never fails an analysis — it degrades to "no
annotation" with a status message (unless you set `required: true`).
Refresh with `tescope annotation show <asm> <te> --refresh` or the GUI's
check button.

Every analysis output directory contains `features_used.json` — the exact
records and provenance used — plus `features.tsv`, so a figure can be
reproduced later without any live API call even if the database changes.

## User-defined features

You may annotate the consensus yourself. User features are always labeled
**user-defined** (separate lane, distinct hatched style, never presented
as database annotations) and are your responsibility.

Coordinates are **0-based, half-open** on the TEscope consensus axis
(identical to every other TEscope coordinate). Formats:

TSV (header required; `label`, `note` optional):

```text
name	start	end	label	note
my_region	100	900	regulatory	from our 2024 mapping
```

YAML/JSON (list of mappings with the same fields), or inline in the
analysis config:

```yaml
annotations:
  enabled: true
  provider: dfam          # or "none" for user features only
  user_features:
    - {name: my_region, start: 100, end: 900, label: regulatory}
  user_features_file: extra_features.tsv   # optional import
  show_in_plots: true
```

Start/end/bounds are validated against the consensus length; malformed
records are rejected with the offending row named. Overlapping features
are allowed (they stack into rows in the lane).

## Feature coverage QC

For each feature the QC table reports, over the **analyzed (post-filter)
copies**:

- *copies overlapping*: copies whose consensus interval intersects the
  feature interval (`cons_start < feature.end and cons_end > feature.start`);
- *copies fully covering*: copies with
  `cons_start <= feature.start and cons_end >= feature.end`.

This is descriptive only — apparent signal depletion inside a feature may
simply mean few fragments reach it. It never alters the profile
calculation.

## Scientific invariance

Annotations are a metadata/visualization layer. Enabling them changes no
projected coordinate, no binned value, no aggregate, no `n_elements`/
`n_observations`, no normalization and no non-annotation QC — this is
pinned by a regression test comparing full result tables with annotations
on and off.

## What TEscope deliberately does NOT do

- invent, guess or interpolate feature names or coordinates;
- assume family-typical structures (not every LINE record has ORF
  annotations; not every LTR family has gag/pol/env records);
- display features from a length-mismatched consensus;
- promote homology evidence into gene/ORF claims;
- reconcile conflicting databases (only Dfam is integrated);
- interpret peaks biologically ("this peak is in X" is your call, using
  the displayed evidence).
