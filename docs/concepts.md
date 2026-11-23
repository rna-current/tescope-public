# How TEscope maps signals onto TE consensus coordinates

## The problem with 0–100 % scaling

Most genomic copies of a transposable element are fragments. LINE-1 is the
classic case: the genome contains a few thousand near-full-length copies and
hundreds of thousands of short 3'-end fragments (5'-truncation during
retrotransposition). If every copy is independently stretched to the same
plot width, a 500 bp 3' fragment and a 6 kb full-length element occupy the
same axis, and positional signal (a promoter mark at the 5' UTR, an RNA
modification hot spot in ORF2) is smeared or destroyed.

## What TEscope does instead

RepeatMasker reports, for every genomic hit, **where that hit aligns on the
repeat consensus model** (the "position in repeat" `begin`/`end`/`(left)`
columns). TEscope keeps those coordinates through parsing and uses them to
place each copy at its true location along the consensus:

```text
Consensus:  5' ======================================= 3'
Copy A:           ===============================
Copy B:                          ================
Copy C:           ========
```

Signal from copy B contributes only to the consensus interval that copy B
actually represents. Aggregate profiles therefore carry a per-bin
denominator: **n_elements** (how many copies cover the bin) and
**n_observations** (how many contributed data). Both are always exported
and shown — a mean over 3 copies is not presented as if it were a mean
over 1,500.

## Strand handling

TE copies lie on both genomic strands. TEscope orients every copy 5'→3'
along the consensus:

- For a `+` strand hit, increasing genomic position = increasing consensus
  position.
- For a `-` strand hit, the **last** genomic base corresponds to the hit's
  **first** consensus base; per-base signal is reversed before binning, and
  interval/point projections are mirrored accordingly.

RepeatMasker already reports consensus coordinates in the consensus's own
orientation for both strands (`+` rows as `begin end (left)`, `C` rows as
`(left) end begin`); TEscope's parsers normalize both to the same canonical
form. Strand correctness is covered by dedicated tests using a synthetic
dataset in which any orientation error inverts a signal ramp.

## Within-copy interpolation (read this)

RepeatMasker `.out` files record only the aligned interval per hit, not the
per-base alignment. Within a hit, TEscope maps genomic bases to consensus
positions **linearly** (affine): insertions/deletions inside the copy are
distributed proportionally across the hit. For typical hits (where indels
are a small fraction of the alignment) the positional error is small
relative to common bin sizes (≥10 bp), but it is an approximation — see
[limitations](limitations.md).

## Consensus model length

The consensus length of a repeat is derived from the annotation itself:
each hit reports `end + (left)`, i.e. how much consensus exists in the
model used by that RepeatMasker run. TEscope uses the **mode** (most
frequent value) across hits of the repeat name — real annotations contain
rare artifact hits with inflated consensus coordinates (in hg38, ~30 L1HS
hits imply a 6.3–7.8 kb model where 1,192 of 1,686 hits agree on 6,155 bp),
and a maximum-based length would misclassify every true full-length copy.
Hits whose coordinates exceed the model length are excluded at analysis
time with an explicit warning. This approach avoids redistributing
consensus sequence libraries (Repbase/Dfam licensing) and guarantees the
coordinates match the annotation actually being projected.

## Binning and aggregation

The consensus axis `[0, L)` is divided into fixed-size bins (default 10 bp;
the last bin may be shorter). Per copy and bin, the mean of the per-base
values mapped into the bin is computed; per bin across copies, TEscope
reports mean, median, SEM, n_elements and n_observations. Interval (peak)
tracks are never converted into curves: they are projected as consensus
intervals, and occurrence statistics report both **peak counts** per bin
and the **fraction of covering copies with a peak** (different
denominators, both stated).
