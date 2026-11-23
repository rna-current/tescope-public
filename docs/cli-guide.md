# CLI guide

The CLI is fully headless (no display server needed) and calls the same
analysis core as the GUI. `tescope <command> --help` shows all options.

## Reference management

```bash
tescope reference list                 # known + installed assemblies
tescope reference fetch hg38           # download + index (one-time)
tescope reference fetch hg38 --force   # re-download
tescope reference add --assembly myasm --repeatmasker my.fa.out \
    [--chrom-sizes my.chrom.sizes] [--notes "..."] [--force]
tescope reference remove myasm
```

`--cache-dir DIR` (any command) or `TESCOPE_CACHE_DIR` overrides the cache
location.

## Browsing TEs

```bash
tescope te list hg38                          # class/family/name + copy counts
tescope te list hg38 --search L1HS
tescope te list hg38 --te-class LINE --family L1 --limit 100
tescope te info hg38 L1HS                     # consensus length, completeness
```

## Consensus feature annotation

```bash
tescope annotation show mm39 IAPEz-int          # status, provenance, features
tescope annotation show mm39 IAPEz-int --refresh  # drop cache + re-fetch
tescope profile ... --annotate                  # attach Dfam features (verified only)
tescope profile ... --user-features my.tsv      # attach user-defined features
```

Config-driven runs use the `annotations:` section (docs/annotations.md).
Annotation failure never fails an analysis unless `required: true`.

## Consensus sequence

```bash
tescope sequence show hg38 AluY                   # source status + provenance
tescope sequence export hg38 AluY --out aluy.fa --range 100-200
tescope sequence search hg38 AluY GGCCGG          # literal; 1-based hits
tescope sequence show myasm MyTE --source library --path my_repeats.fa
```

Dfam is cache-first (`--no-fetch` = offline only, `--refresh` re-fetch);
user-provided libraries are read in place. Only verified,
model-compatible sequences are exported/searched — see
docs/sequence-browser.md. Config-driven runs use the `sequence:` section.

## Running analyses

### From a project file (recommended)

```bash
tescope run analysis.yaml [--out DIR] [--no-figures] [--style-preset publication]
```

Figure appearance is config-driven via the `style:` section
(docs/plot-style.md); `--style-preset` is a quick preset override.

`analysis.yaml` here is schematic: it is a project file **you** created —
saved by the GUI, written by `tescope profile --save-config`, or authored
by hand. TEscope does not install one. To try the CLI immediately,
generate the self-contained demo (`examples/README.md`):

```bash
python examples/make_example_data.py example-data
tescope reference add --assembly demo \
    --repeatmasker example-data/demo.fa.out \
    --chrom-sizes example-data/demo.chrom.sizes
tescope run example-data/analysis.yaml
```

Track paths inside the YAML resolve relative to the YAML's directory. The
config schema is documented in [reproducibility.md](reproducibility.md);
the GUI's *Save project* writes the same format.

### Directly from options

```bash
tescope profile \
    --genome hg38 --te L1HS \
    --signal a.bw --signal b.bedGraph \
    --peaks peaks.narrowPeak \
    --points sites.bed \
    --min-length 1000 \
    --min-consensus-coverage 0.5 \
    --completeness full_length --completeness three_prime_truncated \
    --bin-size 10 --statistic mean \
    --out results/ \
    --save-config analysis.yaml
```

`--save-config` writes the equivalent YAML so an exploratory run becomes a
reproducible one. `--signal/--peaks/--points` are repeatable. Groups and
replicates require the YAML form (`group:`/`replicate:` per track).

## Exit behavior

User-facing problems (missing reference, no copies after filtering,
chromosome mismatch, malformed track — always naming the file) exit with
code 1 and a one-line message, not a stack trace. Progress goes to stderr,
results summary to stdout.

## HPC notes

- Pure CPU, single process; memory scales with (copies × bins) per
  continuous track — a full AluY set at 300 bins is ~100 MB matrix.
- Reference caches are read-only after installation and can be shared
  between users via `TESCOPE_CACHE_DIR`.
- Everything works offline once references are cached.
