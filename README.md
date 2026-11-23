<p align="center">
  <img src="assets/branding/logo/tescope_logo_main@2x.png"
       alt="TEscope" width="300">
</p>

# TEscope

**Fragment-aware visualization and aggregation of genomic signals on
transposable element (TE) consensus coordinates.**

TEscope lets researchers load their own molecular signal files (RNA-seq
coverage, ChIP-seq/CUT&Tag signal, ATAC-seq, m6A enrichment, peak calls,
modification sites, ...) and visualize / aggregate those signals across the
copies of a transposable element family, positioned on the element's
**repeat consensus coordinates** as reported by RepeatMasker.

```text
genomic TE copies
        ↓
RepeatMasker annotation
        ↓
map each TE copy to its repeat-consensus coordinates
        ↓
optional fragment/completeness filtering
        ↓
project user-provided genomic signals
        ↓
aggregate across TE copies
        ↓
visualize profiles / intervals / events
```

## Start here

**Use TEscope (the application)**
[Installation](#installation) · [GUI quick start](#quick-start-gui) ·
[CLI quick start](#quick-start-cli) ·
[runnable example](examples/README.md) ·
[documentation](docs/index.md) · [User Guide](docs/user-guide/README.md)

**Reproduce the published real-data example**
[examples/h3k9me3_l1hs](examples/h3k9me3_l1hs/README.md) · [REPRODUCE.md](REPRODUCE.md)

## Why consensus coordinates?

Most genomic TE copies are fragmented or truncated. Naively rescaling every
copy to 0–100 % of its own length destroys positional information: a 500 bp
3'-end fragment of LINE-1 would be stretched across the whole plot. TEscope
instead uses the repeat-consensus alignment coordinates that RepeatMasker
already reports, so each fragment occupies its **true position** along the
element model:

```text
Consensus:  5' ======================================= 3'
Copy A:           ===============================
Copy B:                          ================
Copy C:           ========
```

## Key features

- Canonical RepeatMasker parser (`.out` files and UCSC `rmsk` tables) that
  **preserves repeat-consensus coordinates**
- Assembly registry with automatic annotation download + local caching
  (hg38, hg19, mm39, mm10) and custom RepeatMasker import for any assembly
- Fragment filtering: genomic length, consensus coverage, positional
  completeness (full-length / 5'-truncated / 3'-truncated / internal)
- Continuous tracks (bigWig, bedGraph), interval tracks (BED, narrowPeak,
  broadPeak), and point tracks (site lists); optional BAM coverage
- Strand-aware projection: every copy is oriented 5'→3' along the consensus
- Aggregation with explicit per-bin denominators (`n_elements`,
  `n_observations`) — no silent averaging over missing data
- Multi-track, group/replicate support
- Cross-platform GUI (PySide6) and fully headless CLI that share the same
  analysis core
- Reproducible YAML project files; publication-quality PNG/SVG/PDF export
  plus machine-readable TSV/CSV results
- **Local-first**: your experimental data never leave your machine. Network
  access is only used to fetch public reference annotation once.

## Installation

Download a build for your platform from the
[releases page](https://github.com/rna-current/tescope-public/releases)
(Windows installer or portable zip, macOS arm64 app, Linux x86_64
bundle; all **unsigned** — see [RELEASE_NOTES.md](RELEASE_NOTES.md) for
the validation scope of each), or install the Python package from the
wheel attached to the same release:

```bash
pip install tescope-0.2.1-py3-none-any.whl            # core + CLI
pip install "tescope-0.2.1-py3-none-any.whl[gui]"     # with the desktop GUI
pip install "tescope-0.2.1-py3-none-any.whl[bam]"     # BAM coverage (Linux/macOS)
```

Every release asset has a `.sha256` companion
(`sha256sum -c <file>.sha256`).

or, for development, with [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/rna-current/tescope-public
cd tescope-public
uv sync --extra gui --extra bam
```

## Quick start (CLI)

```bash
# one-time: fetch the RepeatMasker annotation for an assembly
tescope reference fetch hg38

# list matching TE families
tescope te list hg38 --search L1HS

# run a profile analysis
tescope profile \
    --genome hg38 \
    --te L1HS \
    --min-consensus-coverage 0.5 \
    --signal rnaseq.bw \
    --peaks m6a_peaks.narrowPeak \
    --out results/

# or run a saved / hand-written project file (create one with the GUI's
# "Save project", with --save-config above, or from examples/ — the name
# analysis.yaml is a placeholder, not a bundled file)
tescope run analysis.yaml
```

## Quick start (GUI)

```bash
tescope-gui
```

Choose genome → choose TE → drag & drop signal files → set filters → Run →
inspect profile / peaks / QC → export figure, data, and project file.

## Documentation

See [`docs/`](docs/index.md): installation, GUI guide, CLI guide, input
formats, TE filtering, output interpretation, reproducibility, limitations.

## Scientific limitations (short form)

TEscope deliberately does **not** attempt TE expression quantification or
multi-mapping read rescue. Signals in young TE families are affected by
mappability; consensus positions are interpolated linearly within each
RepeatMasker hit (per-base alignments are not available from `.out` files).
Read [`docs/limitations.md`](docs/limitations.md) before publishing
TEscope-derived figures.

## License

MIT (see `LICENSE`). Reference annotations are downloaded from public
sources at run time and cached locally; they are not redistributed with
TEscope. See `docs/reference-data.md` for provenance and licensing notes.

## Citation

Please cite the archived release you used. Version 0.2.1:

> Ha, H. (2026). TEscope (v0.2.1). Zenodo. https://doi.org/10.5281/zenodo.24458256

Machine-readable metadata: [`CITATION.cff`](CITATION.cff). Report the
exact version (`tescope --version`) with any TEscope-derived figure.
