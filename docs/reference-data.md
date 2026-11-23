# Reference data management, provenance and licensing

## Design: download-and-cache, never redistribute

TEscope does not bundle annotation data. On first use of an assembly the
RepeatMasker annotation is downloaded directly from its public source to a
local cache, the **original file is retained**, and a provenance record is
written. This keeps the installer small and keeps data licensing clean:
each user obtains the data from the source under that source's terms.

## Built-in assemblies

All four built-in assemblies have been fetched, parsed and validated
end-to-end (2026-08-27):

| Assembly | Organism | Source | Hits | Artifact lines skipped |
|---|---|---|---|---|
| hg38 | human (GRCh38) | UCSC goldenPath `hg38.fa.out.gz` | 5,519,906 | 212 |
| hg19 | human (GRCh37) | UCSC goldenPath `hg19.fa.out.gz` | 5,298,122 | 8 |
| mm39 | mouse (GRCm39) | UCSC goldenPath `mm39.fa.out.gz` | 5,320,765 | 6 |
| mm10 | mouse (GRCm38) | UCSC goldenPath `mm10.fa.out.gz` | 5,147,735 | 2 |

Note: repeat naming differs between annotation vintages — mm10's older
RepeatMasker run names mouse L1 subfamilies `L1Md_A`/`L1Md_F`, while mm39
uses `L1MdA_I`/`L1MdTf_I`. Use `tescope te list <assembly> --search` to
find the spelling in your assembly.

The `.fa.out.gz` files are native RepeatMasker output as run by UCSC;
`chrom.sizes` is fetched alongside. UCSC genome annotation downloads are
free for research; the underlying RepeatMasker runs historically used the
Repbase library — which is **why TEscope derives consensus model lengths
from the annotation itself** instead of shipping consensus sequence
libraries. Adding new assemblies is a registry entry
(`tescope/core/reference/registry.py`).

## Custom assemblies

Any RepeatMasker output works, for any organism:

```bash
tescope reference add \
    --assembly customGenome \
    --repeatmasker custom.fa.out \
    --chrom-sizes custom.chrom.sizes \
    --notes "RepeatMasker 4.1.5, Dfam 3.7, our assembly v2"
```

Both native `.out(.gz)` and UCSC `rmsk` table dumps are auto-detected.

## What is recorded per assembly

`<cache>/references/<assembly>/provenance.json`:

- source type (`ucsc-goldenpath` / `user-file`), URL or original path
- retrieval timestamp (UTC), MD5 of the retained source file
- detected format, hit count, malformed-line count
- TEscope version that parsed it, free-text notes

This record is copied into every analysis output (`analysis.yaml`,
`qc.tsv`, `.figure.json`), so a figure can always be traced to the exact
annotation file.

## Cache location

1. `TESCOPE_CACHE_DIR` environment variable, if set
2. otherwise the platform user cache dir (`~/.cache/tescope` on Linux,
   `~/Library/Caches/tescope` on macOS, `%LOCALAPPDATA%\tescope\Cache` on
   Windows)

Layout per assembly: `source/` (retained original), `annotation.npz`
(parsed cache; automatically rebuilt from `source/` if its format
version changes), `hierarchy.npz`, `cons_lengths.tsv`, `chrom.sizes`,
`provenance.json`. `tescope reference remove <assembly>` deletes it.

## Offline use

After `tescope reference fetch`, no network access is used. User signal
data is never uploaded anywhere under any circumstances.
