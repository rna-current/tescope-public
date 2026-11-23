# Troubleshooting

## "Reference annotation for X is not yet downloaded"

Run `tescope reference fetch <assembly>` (CLI) or click *Download
annotation* (GUI). One-time, ~80–120 s per assembly on a fast connection;
cached permanently (~250–290 MB per assembly).

## "Could not load the Qt platform plugin xcb" (Linux GUI)

Your system is missing Qt's X11 runtime libraries. Install
`libxcb-cursor0 libxkbcommon-x11-0 libxcb-icccm4 libxcb-keysyms1
libxcb-image0 libxcb-render-util0` (Debian/Ubuntu) or the equivalent
`xcb-util-*`/`libxkbcommon-x11` RPMs. The CLI does not need any of this.

## "Track ... shares no chromosomes with the selected TE copies"

The track's chromosome names have no overlap with the copies selected on
the chosen assembly. Causes, in order of likelihood: the track was made on
a different assembly; the track covers only chromosomes where no selected
copy resides (e.g. a single-chromosome bigWig combined with a strict
completeness filter); a naming scheme TEscope refuses to guess
(`scaffold_1` etc.). TEscope reconciles `chr1`↔`1` and `chrM`↔`MT`
automatically (disable with `chrom_alias: false`); anything ambiguous is
reported, never guessed.

## "BAM index missing"

Sort and index the file (`samtools sort` + `samtools index`), or enable
automatic indexing (`create_index`). BAM support requires the `[bam]`
extra and is unavailable on Windows.

## "N hits excluded: consensus coordinates exceed the model length"

Not an error. Real RepeatMasker annotations contain rare artifact hits
whose consensus coordinates imply a longer model than the family mode
(e.g. 25 of 1,686 hg38 L1HS hits). Projecting them would distort the last
bin, so they are excluded and counted. See docs/concepts.md,
"Consensus model length".

## "No TE annotations satisfy the current filters"

Your filter combination removed every copy. The GUI filter panel and the
CLI output show per-step counts — relax the step that drops to zero.

## Profile looks empty / near zero

Check `n_observations` in the profile table and the QC tab: your signal
may simply have no data on the selected copies (common for
transcriptome-focused data over intergenic TE copies), or all copies with
signal live on chromosomes absent from the track (see the "copies
skipped" QC line).

## Non-Latin characters in figure labels show as boxes

Publication figures use DejaVu Sans with a fallback list of common
CJK-capable families (Noto Sans KR/CJK, Malgun Gothic, Hiragino, ...).
Boxes mean no listed family is installed — install one (e.g.
`google-noto-sans-kr-fonts` / `fonts-noto-cjk`). If the font IS installed
but boxes persist, matplotlib's font cache predates it: delete
`~/.cache/matplotlib/fontlist*.json` and rerun.

## Large BED/site files are slow to load

Interval and point tracks are held in memory (a 10.6M-line BED loads in
~19 s and ~3 GB RAM). Pre-filter site lists to the chromosomes or regions
of interest, or use bigWig for dense genome-wide signal.

## GUI starts but the reference dropdown shows "(not downloaded)"

That label is informational: select the assembly and click *Download
annotation*. If a download fails mid-way, simply retry — partial files
are never installed (atomic writes).

## Where is my reference cache?

`TESCOPE_CACHE_DIR` if set; otherwise the platform cache dir
(`~/.cache/tescope` on Linux, `~/Library/Caches/tescope` on macOS,
`%LOCALAPPDATA%\tescope\Cache` on Windows). `tescope reference list`
prints it. Removal: `tescope reference remove <assembly>`.
