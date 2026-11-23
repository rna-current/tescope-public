# Installation

## Desktop application (Windows / macOS) — recommended for GUI users

Prebuilt installers are planned for the public release (no Python, conda,
samtools or bedtools required). Until then, install from source below.

## From PyPI / source (all platforms)

Requires Python ≥ 3.10.

```bash
pip install tescope              # core + CLI (headless-safe)
pip install "tescope[gui]"       # + desktop GUI (PySide6)
pip install "tescope[bam]"       # + BAM/CRAM coverage (pysam; Linux/macOS only)
```

Isolated-tool installs (recommended on shared machines / HPC):

```bash
pipx install tescope
# or
uv tool install tescope
```

## Development install

```bash
git clone https://github.com/rna-current/tescope-public
cd tescope-public
uv sync --extra gui --extra bam   # creates .venv from the committed uv.lock
uv run pytest                      # run the test suite
uv run tescope --help
```

## Platform notes

| | bigWig | bedGraph/BED/peaks | BAM/CRAM | GUI |
|---|---|---|---|---|
| Linux | pyBigWig + pybigtools | yes | yes (`[bam]`) | optional |
| macOS | pybigtools | yes | yes (`[bam]`) | yes |
| Windows | pybigtools | yes | **no** (pysam has no Windows wheels) | yes |

- No display server is required for the CLI: `tescope run` works on
  headless HPC nodes.
- The Linux **GUI** needs the standard Qt X11 runtime libraries; on minimal
  servers install (Debian/Ubuntu names) `libxcb-cursor0 libxkbcommon-x11-0
  libxcb-icccm4 libxcb-keysyms1 libxcb-image0 libxcb-render-util0
  libegl1` (RPM: `xcb-util-cursor xcb-util-wm xcb-util-keysyms
  xcb-util-image xcb-util-renderutil libxkbcommon-x11`). Desktop
  distributions usually have them already. If they are missing, Qt reports
  "Could not load the Qt platform plugin xcb".
- On HPC, set `TESCOPE_CACHE_DIR` to a shared/scratch location if your
  home quota is small; annotation caches are ~100–800 MB per assembly.

## First run

Reference annotations are downloaded once and cached:

```bash
tescope reference fetch hg38
```

After that, everything runs offline. Your signal files are never uploaded
anywhere.
