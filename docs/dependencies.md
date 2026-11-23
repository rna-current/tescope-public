# Dependency policy and license audit

Mature, actively maintained libraries only. Development versions are
locked in `uv.lock`; the published package uses conservative lower bounds.
The table below is the **audited state of the locked environment**
(2026-08-27, extracted from installed package metadata).

## Runtime (core)

| Dependency | Locked | Purpose | License | Redistribution concern |
|---|---|---|---|---|
| numpy | 2.5.2 | array math, binning | BSD-3 (+0BSD/MIT/Zlib/CC0 components) | none |
| pandas | 3.0.5 | annotation/result tables | BSD-3 | none |
| PyYAML | 6.0.3 | project files | MIT | none |
| click | 8.5.0 | CLI | BSD-3 | none |
| matplotlib | 3.11.1 | figures (Agg headless + QtAgg embedded) | PSF-based (BSD-compatible) | none |
| pillow | 12.3.0 | matplotlib image backend | MIT-CMU | none |
| platformdirs | 4.11.4 | per-OS cache directory | MIT | none |
| pyBigWig | 0.3.25 | bigWig reading (extra backend, Linux-only) | MIT | none |
| pybigtools | 0.3.0 | bigWig reading (all OSes incl. Windows) | MIT | none |

Two bigWig backends exist deliberately: pybigtools (Rust "bigtools") ships wheels for every platform and is the portable backend; pyBigWig has manylinux wheels only, so it is installed as an additional backend on Linux. One small interface
wraps both, tests assert numeric parity, and the platform choice cannot
change results.

## Optional extras

| Extra | Dependency | Locked | Purpose | License | Bundle notice requirement |
|---|---|---|---|---|---|
| `gui` | PySide6 (+Essentials/Addons, shiboken6) | 6.11.2 | desktop GUI | **LGPL-3.0-only** (of the LGPL/GPL tri-license) | YES — license texts must ship in binary bundles; dynamic linking only |
| `bam` | pysam | 0.24.0 | BAM/CRAM coverage | MIT (htslib MIT/BSD) | none; POSIX only |

## Build/dev only (never shipped as dependencies)

pytest, pytest-cov, coverage, ruff (MIT/Apache-2.0) — test/lint only.
**PyInstaller 6.22.2 is GPL-2.0 with the Bootloader Exception**: it is a
build tool; the exception explicitly permits distributing the frozen
application under TEscope's MIT license. Its bootloader is embedded in the
desktop bundles and is covered by the exception (noted in
`THIRD-PARTY-NOTICES.txt`).

## Icons and theme assets

All GUI icons (application mark, track-kind badges) are original minimal
SVGs authored for TEscope and MIT-licensed with the project
(`src/tescope/gui/theme/icons.py`); no external icon library is bundled.
The theme uses system fonts only.

## Binary-bundle obligations (implemented)

Desktop bundles include a generated `THIRD-PARTY-NOTICES.txt`
(`packaging/collect_licenses.py`) containing the license texts of every
frozen package, which satisfies the PySide6/Qt LGPL-3.0 notice obligation
under dynamic linking. Do not switch PySide6/Qt to static linking without
a fresh license review.

## Audit conclusion

All runtime licenses are permissive and compatible with TEscope's MIT
license; the single copyleft component (PySide6, LGPL-3.0) is an optional
extra used via dynamic linking with notices shipped. **Recommendation:
keep MIT for TEscope itself.** Reference annotation data is downloaded by
each user directly from its source and never redistributed
(docs/reference-data.md) — data licensing is independent of code
licensing. Re-run `packaging/collect_licenses.py` and refresh this table
whenever locked versions change before a release.
