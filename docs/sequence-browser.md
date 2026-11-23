# Interactive consensus browser & sequence support

TEscope's result view is an interactive consensus browser: zoom, pan and
navigate along the TE consensus axis, with an optional base-level
consensus sequence track. All of it is **display-only** — no zoom, pan,
sequence source, cursor, search or selection can change any projected,
aggregated or exported scientific value (pinned by regression tests).

## Navigation

- **Mouse wheel** zooms in/out centered near the cursor, clamped to the
  consensus (minimum window 10 bp).
- **Toolbar** (unchanged matplotlib navigation): box zoom, pan by drag,
  Home, and Back/Forward view history.
- **Go to** accepts a base position (`4127`) or range (`4000-4500`);
  **Full view** returns to the whole consensus.
- Every aligned panel — signal, annotation lanes, copy support, mapping
  fraction, sequence — shares one x axis: they always show **exactly the
  same consensus interval** and move together.
- Hovering shows a readout of the position, base (when a verified
  sequence exists), signal values, contributing copies and any feature at
  that position; clicking marks a position across all panels. Only fields
  that actually exist are shown.

## Coordinate convention

Internally TEscope stays 0-based half-open (bins, feature intervals,
array indices). **Everything user-facing about bases — cursor readout,
Go to, search hits, selections, FASTA headers — is 1-based inclusive**,
the convention of sequence tools. Base *b* occupies axis interval
[*b*−1, *b*). The conversion lives in one tested module
(`tescope/core/sequence/coords.py`); no other code converts.

## The sequence must match the coordinate model

The displayed sequence MUST correspond to the same consensus coordinate
model TEscope projects onto. **A matching family name never implies a
matching model**: models change and split between database releases.
Every candidate sequence passes the same non-inferential compatibility
gate the feature-annotation layer uses — exact length match against the
annotation-derived model — plus a check that the actually fetched
sequence has the promised length.

- **verified** — positively matched; all sequence functions active.
- **incompatible** — positively determined not to match (e.g. current
  Dfam represents L1HS as split 5'/3' fragments while the RepeatMasker
  model is full-length 6,155 bp). The sequence is **not displayed** and
  coordinate-exact functions are disabled with the reason shown. TEscope
  never concatenates fragments, infers missing bases, or remaps an
  incompatible model to force a fit.
- **unverified** — compatibility could not be established; treated like
  incompatible for display.
- **unavailable** — no source, no matching family, or offline with
  nothing cached. Normal analysis is never blocked by this.

## Sequence sources

Selected in the GUI Sequence panel or the config `sequence:` section:

| source | meaning |
|---|---|
| `dfam` (Automatic) | official Dfam API, cache-first (offline after the first successful lookup) |
| `library` | user-provided multi-record FASTA (e.g. a RepeatMasker-associated or custom research library) |
| `fasta` | user-provided FASTA (single family) |
| `none` | no sequence |

**Dfam** lookups reuse the existing Dfam integration (name + Repbase
alias matching). Retrieved sequences are cached under
`sequences/dfam/<name>__<model_length>.json` with full provenance
(accession, release, retrieval date, endpoint, length, sha256,
compatibility result). The cache key includes the expected model length —
never the family name alone — so a changed annotation model can never
silently reuse a sequence cached for a different model.

**User-provided libraries** are read in place: standard FASTA (plain or
gzip), including RepeatMasker `name#class/family` headers. Resolution is
careful and never guesses — exact header match first, then the pre-`#`
name; ambiguous or missing names fail with the candidates listed. When
the RepeatMasker annotation came from a custom run, pointing TEscope at
the **exact consensus library used for that run** is the best source and
takes precedence over substituting the newest Dfam model; if name and
exact model length verify, sequence-aware visualization activates with
provenance marked user-provided (never presented as Dfam-curated).

The config records the file path (not the sequence contents), plus
checksum and provenance for reproducibility snapshots. A missing file at
run time is a clear, explained "unavailable" — never a silent fallback.

## Sequence track

Enabled with *Show sequence* (config `sequence.show_track`) when a
verified sequence exists. Rendering is zoom-dependent:

- **full view**: a compact per-base color strip (presence/overview);
- **base-level zoom** (roughly ≤ a few hundred visible bases, based on
  available pixel width): actual A/C/G/T/N letters aligned to their
  positions. TEscope never prints thousands of unreadable letters.

Base colors are muted and are never the only carrier of identity — the
letter itself is visible at readable zoom, so grayscale output stays
interpretable. Exports (SVG/PDF keep letters as vector text) use the
same zoom-dependent rules for the configured view window
(`style.x_min`/`x_max`), so a saved config reproduces the same
representation headlessly.

## Copy, export, search

- **Copy sequence** puts the visible region on the clipboard as FASTA
  with its 1-based range in the header.
- **Export FASTA…** writes the full consensus or the visible region;
  headers carry family, provider/accession, source release, coordinate
  range, source checksum and TEscope version.
- **Search** is literal only: IUPAC letters, case-insensitive, forward
  strand (`N` matches any base; CLI `--revcomp` adds clearly flagged
  reverse-complement hits). Results are coordinate hits — selecting one
  zooms to it. No motif prediction, ORF calling, domain inference or
  translated annotation, ever: if TEscope can't read it verbatim from an
  authoritative source, it doesn't claim it (see
  [annotations](annotations.md)).

All of these require a **verified** sequence; otherwise the controls are
visible but disabled, with a tooltip explaining exactly why (no source /
incompatible model / file missing / offline and uncached) — the same
capability-driven behavior as the mapping controls. Sequence capability
is orthogonal to signal-track capability: a bigWig analysis can carry a
verified sequence track, and mapping classification stays a BAM/CRAM
matter.

## CLI

```bash
tescope sequence show   hg38 AluY                     # status + provenance
tescope sequence show   hg38 L1HS                     # honest: incompatible
tescope sequence export hg38 AluY --out aluy.fa --range 100-200
tescope sequence search hg38 AluY GGCCGG
tescope sequence show   myasm MyTE --source library --path my_repeats.fa
```

`--no-fetch` restricts Dfam to the cache; `--refresh` re-fetches.

## Licensing

- **Dfam** data is released under CC0 1.0; TEscope retrieves and caches
  it via the official API according to those terms.
- **TEscope does not distribute third-party licensed repeat libraries**
  (Repbase is an example of such a licensed resource). TEscope never
  bundles, redistributes, embeds or automatically downloads licensed
  content, and provides no mechanism to bypass a provider's licensing.
- **Users may load locally obtained libraries** through the generic
  user-provided library interface, subject to the terms of their
  respective licenses; the file is read in place, stays on the local
  machine, and is never uploaded anywhere (TEscope is local-first — a
  Dfam lookup sends only the repeat name/model length needed for the
  query). Users remain responsible for complying with their libraries'
  license terms. This is an operational description, not legal advice.

## Provenance

Every used sequence records provider, family, accession (when any),
source release, retrieval date or local path, length, sha256 checksum,
the expected model length, and the compatibility outcome — in the QC
tab, in exported `qc.tsv`, and in every figure's `.figure.json` sidecar
(`sequence_provenance`).
