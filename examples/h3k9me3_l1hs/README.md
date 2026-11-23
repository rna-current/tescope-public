# Real-data example: H3K9me3 × L1HS (hg38)

Public data: ENCODE experiment ENCSR000APE (K562) — H3K9me3 fold change
over control (bigWig ENCFF601JGK) and replicated peaks (ENCFF963GZJ),
profiled across the 1,661 hg38 RepeatMasker L1HS records. Units are
ENCODE's upstream-normalized fold change; TEscope reports per-copy means
with per-bin copy support (`n_elements`).

The ENCODE files and the UCSC RepeatMasker annotation are **not
redistributed** here. `fetch_data.sh` downloads the two ENCODE files
from encodeproject.org and checks their MD5 sums; TEscope downloads the
hg38 annotation from UCSC on first use
(`hg38.fa.out.gz`, md5 `820796b72974d077281b3386c4fb8295` at the time
the canonical outputs were produced; recorded in `outputs/qc.tsv`).

| file | purpose |
|---|---|
| `fetch_data.sh` | download + checksum the two ENCODE inputs into `data/` |
| `h3k9me3_l1hs_all.yaml` | configuration: all 1,661 L1HS records |
| `h3k9me3_l1hs_subset.yaml` | same analysis with the 20-locus subset applied (the state shown in the manuscript figure) |
| `select_demo_loci.py` | the documented subset criterion: the 20 full-length copies with the highest mean fold change over observed bins. **Exploratory** — chosen to demonstrate the subset workflow, enriched by construction, not a biological finding |
| `demo_selection.yaml` | the resulting locus list (importable in the Loci tab) |
| `outputs/` | canonical output tables (and the app-rendered profile figures) for the all-copies analysis; these are the tables underlying the manuscript figure |
| `figure/` | capture and composition code for the manuscript figure |

## Reproduce

Run from this directory (paths in the configurations are relative to it):

```bash
./fetch_data.sh
tescope run h3k9me3_l1hs_all.yaml        # writes results/all/
python select_demo_loci.py results/all demo_selection.regenerated.yaml
tescope run h3k9me3_l1hs_subset.yaml     # writes results/subset/
```

The output tables regenerate byte-identically for the same TEscope
version and annotation snapshot; compare with `outputs/`:

```bash
for f in outputs/*.tsv; do
  case "$f" in
    *qc.tsv) continue ;;                                    # see below
    *subset_*) cmp "$f" "results/subset/$(basename "$f")" ;;
    *) cmp "$f" "results/all/$(basename "$f")" ;;
  esac
done
```

`qc.tsv` records the absolute paths of the input files on the machine
that ran the analysis, so its two `file` rows differ between machines
(here they are given relative to this directory); every other row is
reproducible.

Figure panels (needs a display; unedited `QWidget.grab()` captures of
the running application, then side-by-side composition):

```bash
python figure/capture_figure_panels.py figure/out
python figure/compose_figure.py figure/out figure/out   # needs Pillow and a bold sans-serif TTF (FIGURE_FONT=/path/to/font.ttf)
```

Interpretation caveats: see `../../docs/limitations.md` (multimapping is
decided upstream; coordinates are annotation-model coordinates; read
profiles together with their copy-support tracks).
