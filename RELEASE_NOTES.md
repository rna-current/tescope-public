# TEscope 0.2.1 — release notes

Stable release; **use this version instead of 0.2.0.** Two defects of
0.2.0 are corrected:

1. **GUI reproducibility.** A project opened with *File ▸ Open project*
   lost the settings that have no GUI control — most importantly the
   aggregation bin size — and ran with the defaults (10-bp bins), so a
   GUI run did not reproduce the command-line run of the same project
   file. Opened projects now keep these settings (regression test
   added). Command-line runs were never affected.
2. **Linux bundle.** The 0.2.0 bundle was built on a newer build image
   than documented and required glibc 2.43 (it does not start on
   Rocky/RHEL 10 or Ubuntu 24.04). The bundle is now built on a pinned
   image (glibc 2.39) and the build fails if that requirement is
   exceeded. A bundle built by this workflow from this exact commit was
   started on Rocky Linux 10.1 and driven through the real-data example
   (open project, run, Loci tab with applied subset, figure export) by
   synthetic keyboard and mouse input before the version was tagged.

Archived at
https://doi.org/10.5281/zenodo.24458256 (source, offline acceptance
pack, User Guide).

## Highlights of the 0.2 line

- **TE metagene analysis** on RepeatMasker consensus coordinates:
  fragment- and strand-aware projection, explicit per-bin copy-support
  denominators, completeness-class comparison, per-copy heatmaps.
- **Loci tab + Locus Inspector**: every analyzed copy is addressable
  (stable coordinate IDs; missing evidence is labeled "unavailable",
  never zero); explicit subset analysis re-aggregates the original
  per-copy values with both denominators always shown.
- **User TE feature files** (template TSV in consensus or genome
  coordinates, BED, GFF3 under their specification conventions;
  amino-acid coordinates rejected) and **upstream per-locus estimate
  import** with declared, enforced semantics (never positional signal;
  totals never summed with unique components). The estimate import is
  validated on schema and synthetic data only.
- Reproducible single-file configurations (schema 2); provenance
  sidecars on every figure export.

## Validation scope

| evidence | what it is | what it is not |
|---|---|---|
| Test suite (unit + synthetic golden end-to-end) | run on Linux, Windows and macOS runners (Python 3.10 and 3.12) for the release commit | — |
| Linux x86_64 bundle | CI build on a pinned glibc-2.39 image with an automated glibc-requirement check; start under Xvfb with the real xcb platform; a pre-release build of the same commit was run on Rocky Linux 10.1 (scripted input, not a human session) | other distributions untested; needs glibc ≥ 2.39 |
| Windows x86_64 installer / portable zip | CI build; the installer is executed silently and the installed app is started (automated, offscreen) | **no interactive validation on a physical Windows PC** |
| macOS arm64 app | CI build and offscreen start; ad-hoc signature | **not notarized; no interactive validation on a physical Mac; Intel Macs unsupported** |
| Wheel / sdist | clean-environment install; relocated offline acceptance pack (all canonical tables must match) | — |

Builds are **unsigned**: expect SmartScreen / Gatekeeper warnings
(macOS: right-click ▸ Open once).

## Known limitations

- Multimapping ambiguity is inherent to upstream alignment; TEscope
  exposes NH/MAPQ observation classes, it cannot recover true origins.
- Consensus coordinates are annotation-model coordinates; the
  name+length compatibility gate is not sequence-identity proof.
- Consensus positions are interpolated linearly within each
  RepeatMasker hit.
- Scientific interpretation requires the per-bin support tracks; see
  `docs/limitations.md`.

## Verify a download

    sha256sum -c TEscope-0.2.1-<platform>.<ext>.sha256

Then run the synthetic acceptance pack
(`tescope-acceptance-pack-0.2.1.tar.gz`): `python verify.py` inside the
unpacked directory, with the `tescope` CLI installed. All canonical
tables must match.

Developer: Hongseok Ha · License: MIT · Citation: `CITATION.cff`
