"""TEscope command-line interface.

Fully headless: no display server required. The CLI calls exactly the same
:func:`tescope.core.analysis.run_analysis` the GUI uses.
"""

from __future__ import annotations

import sys
from pathlib import Path

import click

from tescope import __version__
from tescope.plotting.plotstyle import PRESETS

_USER_ERRORS: tuple[type[Exception], ...] = ()


def _lazy_errors():
    global _USER_ERRORS
    from tescope.config.model import ConfigError
    from tescope.core.analysis import AnalysisError
    from tescope.core.chromnames import ChromMismatchError
    from tescope.core.reference.store import ReferenceError

    _USER_ERRORS = (AnalysisError, ReferenceError, ConfigError, ChromMismatchError)
    return _USER_ERRORS


def _fail(e: Exception) -> None:
    raise click.ClickException(str(e))


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="tescope")
def cli():
    """TEscope: signals on transposable element consensus coordinates."""


# --------------------------------------------------------------------------
# reference management
# --------------------------------------------------------------------------


@cli.group()
def reference():
    """Manage reference assemblies and RepeatMasker annotations."""


@reference.command("list")
@click.option("--cache-dir", type=click.Path(), default=None, help="Cache directory override.")
def reference_list(cache_dir):
    """List known and installed assemblies."""
    from tescope.core.reference.store import ReferenceStore

    store = ReferenceStore(cache_dir)
    installed = set(store.installed_assemblies())
    click.echo(f"cache: {store.cache_dir}")
    click.echo(f"{'assembly':<12} {'status':<12} description")
    for name, spec in store.known_assemblies().items():
        status = "installed" if name in installed else "available"
        click.echo(f"{name:<12} {status:<12} {spec.description} ({spec.organism})")
    for name in sorted(installed - set(store.known_assemblies())):
        click.echo(f"{name:<12} {'installed':<12} custom assembly")


@reference.command("fetch")
@click.argument("assembly")
@click.option("--cache-dir", type=click.Path(), default=None)
@click.option("--force", is_flag=True, help="Re-download even if already installed.")
def reference_fetch(assembly, cache_dir, force):
    """Download and index the RepeatMasker annotation for a known assembly."""
    from tescope.core.reference.store import ReferenceStore

    store = ReferenceStore(cache_dir)

    def progress(stage, done, total):
        if stage == "download" and total:
            pct = 100 * done // total
            click.echo(f"\r  downloading... {done // (1 << 20)} MB ({pct}%)", nl=False, err=True)
        elif stage in ("parse", "index"):
            click.echo(f"\n  {stage}ing...", nl=False, err=True)

    try:
        prov = store.fetch(assembly, progress=progress, force=force)
    except _lazy_errors() as e:
        _fail(e)
    click.echo("", err=True)
    click.echo(f"installed {assembly}: {prov.extra.get('n_hits', '?')} hits from {prov.source_url}")


@reference.command("add")
@click.option("--assembly", required=True, help="Name for the (custom) assembly.")
@click.option(
    "--repeatmasker",
    "rm_path",
    required=True,
    type=click.Path(exists=True),
    help="RepeatMasker .out(.gz) or UCSC rmsk table for the assembly.",
)
@click.option("--chrom-sizes", type=click.Path(exists=True), default=None)
@click.option("--notes", default="", help="Free-text provenance notes.")
@click.option("--cache-dir", type=click.Path(), default=None)
@click.option("--force", is_flag=True, help="Replace an existing installation.")
def reference_add(assembly, rm_path, chrom_sizes, notes, cache_dir, force):
    """Import a user-provided RepeatMasker annotation (custom assemblies)."""
    from tescope.core.reference.store import ReferenceStore

    store = ReferenceStore(cache_dir)
    try:
        prov = store.add_custom(
            assembly, rm_path, chrom_sizes_path=chrom_sizes, notes=notes, force=force
        )
    except _lazy_errors() as e:
        _fail(e)
    click.echo(
        f"installed {assembly}: {prov.extra.get('n_hits', '?')} hits "
        f"({prov.source_format}) from {rm_path}"
    )


@reference.command("remove")
@click.argument("assembly")
@click.option("--cache-dir", type=click.Path(), default=None)
@click.confirmation_option(prompt="Remove this assembly from the local cache?")
def reference_remove(assembly, cache_dir):
    """Remove an installed assembly from the cache."""
    from tescope.core.reference.store import ReferenceStore

    ReferenceStore(cache_dir).remove(assembly)
    click.echo(f"removed {assembly}")


# --------------------------------------------------------------------------
# TE browsing
# --------------------------------------------------------------------------


@cli.group()
def te():
    """Browse transposable element annotations."""


@te.command("list")
@click.argument("assembly")
@click.option("--search", default="", help="Substring match on repeat name.")
@click.option("--te-class", "te_class", default="", help="Filter by repeat class (e.g. LINE).")
@click.option("--family", default="", help="Filter by repeat family (e.g. L1).")
@click.option("--limit", default=50, show_default=True)
@click.option("--cache-dir", type=click.Path(), default=None)
def te_list(assembly, search, te_class, family, limit, cache_dir):
    """List repeat names (class/family hierarchy) with copy counts."""
    from tescope.core.reference.store import ReferenceStore

    store = ReferenceStore(cache_dir)
    try:
        h = store.hierarchy(assembly)
    except _lazy_errors() as e:
        _fail(e)
    if te_class:
        h = h[h["repeat_class"].astype(str).str.lower() == te_class.lower()]
    if family:
        h = h[h["repeat_family"].astype(str).str.lower() == family.lower()]
    if search:
        h = h[h["repeat_name"].astype(str).str.contains(search, case=False, regex=False)]
    total = len(h)
    click.echo(f"{'class':<16} {'family':<16} {'name':<20} {'copies':>8}")
    for r in h.head(limit).itertuples():
        click.echo(
            f"{r.repeat_class:<16} {r.repeat_family:<16} {r.repeat_name:<20} {r.n_copies:>8}"
        )
    if total > limit:
        click.echo(f"... {total - limit} more (use --limit/--search to narrow)")


@te.command("info")
@click.argument("assembly")
@click.argument("repeat_name")
@click.option("--cache-dir", type=click.Path(), default=None)
def te_info(assembly, repeat_name, cache_dir):
    """Copy counts, consensus length and completeness distribution for one repeat."""
    from tescope.core.annotation.schema import add_derived_columns
    from tescope.core.filtering.filters import classify_completeness
    from tescope.core.reference.store import ReferenceStore

    store = ReferenceStore(cache_dir)
    try:
        sub = store.select_repeat(assembly, repeat_name)
    except _lazy_errors() as e:
        _fail(e)
    sub = add_derived_columns(sub, store.consensus_lengths(assembly))
    comp = classify_completeness(sub)
    r0 = sub.iloc[0]
    click.echo(f"{repeat_name} ({r0['repeat_class']}/{r0['repeat_family']}) in {assembly}")
    click.echo(f"  copies:            {len(sub):,}")
    click.echo(f"  consensus length:  {int(sub['cons_model_length'].iloc[0]):,} bp")
    click.echo(f"  median length:     {int(sub['genome_length'].median()):,} bp")
    click.echo(f"  median coverage:   {sub['consensus_coverage'].median():.2f}")
    click.echo("  completeness:")
    for cls, n in comp.value_counts().items():
        click.echo(f"    {cls:<24} {n:>8,}")


# --------------------------------------------------------------------------
# consensus-feature annotation
# --------------------------------------------------------------------------


@cli.group()
def annotation():
    """Consensus-feature annotations (database + user-defined)."""


@annotation.command("show")
@click.argument("assembly")
@click.argument("repeat_name")
@click.option("--fetch/--no-fetch", default=True, show_default=True,
              help="Allow a network lookup when nothing is cached.")
@click.option("--refresh", is_flag=True, help="Discard the cache entry and re-fetch.")
@click.option("--cache-dir", type=click.Path(), default=None)
def annotation_show(assembly, repeat_name, fetch, refresh, cache_dir):
    """Show the Dfam annotation state, provenance and features for a repeat."""
    from tescope.config.model import AnnotationConfig
    from tescope.core.features.cache import cache_path, drop_cached
    from tescope.core.features.resolve import resolve_features
    from tescope.core.reference.store import ReferenceStore

    store = ReferenceStore(cache_dir)
    try:
        store.select_repeat(assembly, repeat_name)
        cons_length = int(store.consensus_lengths(assembly)[repeat_name])
    except _lazy_errors() as e:
        _fail(e)
    except KeyError:
        raise click.ClickException(
            f"no consensus length for '{repeat_name}' in {assembly}"
        ) from None
    if refresh:
        drop_cached(repeat_name, cons_length, "dfam", store.cache_dir)
    fs = resolve_features(
        AnnotationConfig(enabled=True, provider="dfam"),
        repeat_name, cons_length,
        cache_dir=store.cache_dir, allow_network=fetch,
    )
    click.echo(f"{repeat_name} ({assembly}, consensus {cons_length:,} bp)")
    click.echo(f"  status: {fs.status}")
    if fs.provenance:
        pv = fs.provenance
        click.echo(f"  source: {pv.source} {pv.database_version} ({pv.database_release_date})")
        if pv.accession:
            click.echo(f"  record: {pv.accession} '{pv.record_name}' v{pv.record_version} "
                       f"(matched by {pv.matched_by})")
        click.echo(f"  compatibility: {pv.compatibility} — {pv.compatibility_detail}")
        click.echo(f"  retrieved: {pv.retrieved_at}")
        click.echo(f"  cache: {cache_path(repeat_name, cons_length, 'dfam', store.cache_dir)}")
    for f in fs.database_features:
        click.echo(f"  [{f.type}] {f.start}-{f.end}  {f.name}"
                   + (f"  ({f.note})" if f.note else ""))
    if not fs.database_features:
        click.echo("  (no database features to display)")


# --------------------------------------------------------------------------
# consensus sequence
# --------------------------------------------------------------------------


@cli.group()
def sequence():
    """Consensus sequences (Dfam automatic or user-provided libraries)."""


def _sequence_options(f):
    opts = [
        click.option("--source", type=click.Choice(["dfam", "library", "fasta"]),
                     default="dfam", show_default=True,
                     help="Sequence source (library/fasta = user-provided local file)."),
        click.option("--path", "seq_path", type=click.Path(exists=True), default=None,
                     help="User-provided FASTA/library file (library/fasta sources)."),
        click.option("--family", default="", help="Library record name "
                     "(default: the repeat name)."),
        click.option("--library-name", default="",
                     help="Optional provenance label for a user-provided library."),
        click.option("--fetch/--no-fetch", default=True, show_default=True,
                     help="Allow a Dfam network lookup when nothing is cached."),
        click.option("--refresh", is_flag=True,
                     help="Discard the Dfam cache entry and re-fetch."),
        click.option("--cache-dir", type=click.Path(), default=None),
    ]
    for opt in reversed(opts):
        f = opt(f)
    return f


def _resolve_cli_sequence(assembly, repeat_name, source, seq_path, family,
                          library_name, fetch, refresh, cache_dir):
    from tescope.config.model import SequenceConfig
    from tescope.core.reference.store import ReferenceStore
    from tescope.core.sequence.cache import drop_cached
    from tescope.core.sequence.resolve import resolve_sequence

    store = ReferenceStore(cache_dir)
    try:
        store.select_repeat(assembly, repeat_name)
        cons_length = int(store.consensus_lengths(assembly)[repeat_name])
    except _lazy_errors() as e:
        _fail(e)
    except KeyError:
        raise click.ClickException(
            f"no consensus length for '{repeat_name}' in {assembly}"
        ) from None
    if source in ("library", "fasta") and not seq_path:
        raise click.UsageError(f"--source {source} requires --path")
    if refresh:
        drop_cached(repeat_name, cons_length, "dfam", store.cache_dir)
    cfg = SequenceConfig(enabled=True, source=source, path=str(seq_path or ""),
                         family=family, library_name=library_name)
    res = resolve_sequence(cfg, repeat_name, cons_length,
                           cache_dir=store.cache_dir, allow_network=fetch)
    return res, cons_length


def _echo_sequence_status(repeat_name, assembly, cons_length, res) -> None:
    click.echo(f"{repeat_name} ({assembly}, consensus model {cons_length:,} bp)")
    click.echo(f"  status: {res.status} — {res.detail}")
    if res.provenance is not None:
        pv = res.provenance
        if pv.provider:
            click.echo(f"  provider: {pv.provider}"
                       + (f" ({pv.source_version})" if pv.source_version else ""))
        if pv.accession:
            click.echo(f"  record: {pv.accession} '{pv.record_name}' "
                       f"(matched by {pv.matched_by})")
        elif pv.record_name:
            click.echo(f"  record: '{pv.record_name}' (matched by {pv.matched_by})")
        if pv.source_path:
            click.echo(f"  file: {pv.source_path}")
        if pv.sequence_length is not None:
            click.echo(f"  sequence length: {pv.sequence_length:,} bp")
        if pv.checksum:
            click.echo(f"  sha256: {pv.checksum}")
        if pv.retrieved_at:
            click.echo(f"  retrieved: {pv.retrieved_at}")


@sequence.command("show")
@click.argument("assembly")
@click.argument("repeat_name")
@_sequence_options
def sequence_show(assembly, repeat_name, source, seq_path, family, library_name,
                  fetch, refresh, cache_dir):
    """Show the consensus-sequence source status and provenance."""
    res, cons_length = _resolve_cli_sequence(
        assembly, repeat_name, source, seq_path, family, library_name,
        fetch, refresh, cache_dir)
    _echo_sequence_status(repeat_name, assembly, cons_length, res)
    if res.verified:
        s = res.sequence.sequence
        click.echo(f"  first 60 bp: {s[:60]}")
        click.echo("  (base positions are 1-based inclusive)")


@sequence.command("export")
@click.argument("assembly")
@click.argument("repeat_name")
@click.option("--out", "out_path", type=click.Path(), required=True,
              help="Output FASTA file.")
@click.option("--range", "base_range", default="",
              help="1-based inclusive base range, e.g. 4000-4500 (default: full).")
@_sequence_options
def sequence_export(assembly, repeat_name, out_path, base_range, source,
                    seq_path, family, library_name, fetch, refresh, cache_dir):
    """Export the verified consensus sequence (or a region) as FASTA."""
    from tescope.core.sequence.coords import RangeError, parse_position_or_range
    from tescope.core.sequence.ops import fasta_text

    res, cons_length = _resolve_cli_sequence(
        assembly, repeat_name, source, seq_path, family, library_name,
        fetch, refresh, cache_dir)
    if not res.verified:
        _echo_sequence_status(repeat_name, assembly, cons_length, res)
        raise click.ClickException(
            f"no verified consensus sequence to export ({res.status}: {res.detail})"
        )
    b1 = b2 = None
    if base_range:
        try:
            b1, b2 = parse_position_or_range(base_range, res.sequence.length)
        except RangeError as e:
            raise click.ClickException(str(e)) from None
    text = fasta_text(res.sequence, b1, b2, tescope_version=__version__)
    Path(out_path).write_text(text, encoding="utf-8")
    n = (b2 - b1 + 1) if b1 is not None else res.sequence.length
    click.echo(f"wrote {n:,} bp to {out_path} (coordinates 1-based inclusive)")


@sequence.command("search")
@click.argument("assembly")
@click.argument("repeat_name")
@click.argument("query")
@click.option("--revcomp", is_flag=True,
              help="Also report reverse-complement matches (flagged).")
@_sequence_options
def sequence_search(assembly, repeat_name, query, revcomp, source, seq_path,
                    family, library_name, fetch, refresh, cache_dir):
    """Literal sequence search (case-insensitive; N = any base). Hits are
    1-based inclusive consensus coordinates — no motif prediction."""
    from tescope.core.sequence.model import SequenceError
    from tescope.core.sequence.ops import search_literal

    res, cons_length = _resolve_cli_sequence(
        assembly, repeat_name, source, seq_path, family, library_name,
        fetch, refresh, cache_dir)
    if not res.verified:
        _echo_sequence_status(repeat_name, assembly, cons_length, res)
        raise click.ClickException(
            f"no verified consensus sequence to search ({res.status})"
        )
    try:
        hits = search_literal(res.sequence, query, include_revcomp=revcomp)
    except SequenceError as e:
        raise click.ClickException(str(e)) from None
    click.echo(f"{len(hits)} hit(s) for '{query.upper()}' "
               "(1-based inclusive consensus coordinates)")
    for h in hits:
        click.echo(f"  {h.base_start:>8} - {h.base_end:<8} strand {h.strand}")


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------


def _echo_progress(msg: str, step: int, total: int) -> None:
    click.echo(f"[{step}/{total}] {msg}", err=True)


def _execute(config, out_dir: Path, no_figures: bool, cache_dir=None) -> None:
    from tescope.core.analysis import run_analysis
    from tescope.core.reference.store import ReferenceStore
    from tescope.export.figures import render_result_figures
    from tescope.export.tables import export_tables

    store = ReferenceStore(cache_dir)
    try:
        result = run_analysis(config, store, progress=_echo_progress)
    except _lazy_errors() as e:
        _fail(e)
    for w in result.warnings:
        click.secho(f"warning: {w}", fg="yellow", err=True)
    files = export_tables(result, out_dir)
    if not no_figures:
        files += render_result_figures(result, out_dir)
    click.echo(f"\n{result.repeat_name}: {result.n_selected}/{result.n_annotated} copies analyzed")
    for label, n in result.filter_report.to_rows():
        click.echo(f"  {label:<38} {n:>8,}")
    click.echo(f"\nwrote {len(files)} files to {out_dir}/")


@cli.command("run")
@click.argument("config_file", type=click.Path(exists=True))
@click.option("--out", "out_dir", type=click.Path(), default=None,
              help="Output directory (default: from config).")
@click.option("--style-preset",
              type=click.Choice([p for p in PRESETS if p != "custom"]),
              default=None,
              help="Override the plot-style preset (rendering only; "
                   "individual style settings in the config still apply).")
@click.option("--no-figures", is_flag=True, help="Export tables only.")
@click.option("--cache-dir", type=click.Path(), default=None)
def run_cmd(config_file, out_dir, style_preset, no_figures, cache_dir):
    """Run a saved analysis config (YAML/JSON) headlessly."""
    from tescope.config.io import load_config

    try:
        config = load_config(config_file)
    except _lazy_errors() as e:
        _fail(e)
    if style_preset:
        config.style.preset = style_preset
    out = Path(out_dir or config.output.dir)
    _execute(config, out, no_figures, cache_dir)


@cli.command("profile")
@click.option("--genome", required=True, help="Assembly name (e.g. hg38, or a custom one).")
@click.option("--te", "te_name", required=True, help="Repeat name (e.g. L1HS).")
@click.option("--signal", "signals", multiple=True, type=click.Path(exists=True),
              help="Continuous track (bigWig/bedGraph/BAM); repeatable.")
@click.option("--peaks", "peaks", multiple=True, type=click.Path(exists=True),
              help="Interval track (BED/narrowPeak/broadPeak); repeatable.")
@click.option("--points", "points", multiple=True, type=click.Path(exists=True),
              help="Point track (BED-like site list); repeatable.")
@click.option("--min-length", type=int, default=None, help="Minimum genomic length (bp).")
@click.option("--min-consensus-coverage", type=float, default=None,
              help="Minimum consensus coverage (0-1).")
@click.option("--completeness", multiple=True,
              type=click.Choice(["full_length", "five_prime_truncated",
                                 "three_prime_truncated", "internal_fragment"]),
              help="Keep only these completeness classes; repeatable.")
@click.option("--bin-size", type=int, default=10, show_default=True,
              help="Consensus bin size (bp).")
@click.option("--statistic", type=click.Choice(["mean", "median"]), default="mean",
              show_default=True)
@click.option("--out", "out_dir", type=click.Path(), default="tescope-results",
              show_default=True)
@click.option("--save-config", type=click.Path(), default=None,
              help="Also write the equivalent analysis YAML here.")
@click.option("--annotate", is_flag=True,
              help="Attach Dfam consensus features (verified-compatible only).")
@click.option("--user-features", "user_features", type=click.Path(exists=True),
              default=None, help="TSV/YAML file of user-defined consensus features.")
@click.option("--no-figures", is_flag=True)
@click.option("--cache-dir", type=click.Path(), default=None)
def profile_cmd(genome, te_name, signals, peaks, points, min_length,
                min_consensus_coverage, completeness, bin_size, statistic,
                out_dir, save_config, annotate, user_features, no_figures, cache_dir):
    """Run a profile analysis directly from command-line options."""
    from tescope.config.io import save_config as save_cfg
    from tescope.config.model import (
        AggregationConfig,
        AnalysisConfig,
        RepeatSelection,
        TrackConfig,
    )
    from tescope.core.filtering.filters import FilterConfig

    if not (signals or peaks or points):
        raise click.UsageError("provide at least one --signal, --peaks or --points file")
    tracks = (
        [TrackConfig(file=str(p)) for p in signals]
        + [TrackConfig(file=str(p), kind="interval") for p in peaks]
        + [TrackConfig(file=str(p), kind="point") for p in points]
    )
    config = AnalysisConfig(
        assembly=genome,
        repeat=RepeatSelection(name=te_name),
        tracks=tracks,
        filters=FilterConfig(
            min_genome_length=min_length,
            min_consensus_coverage=min_consensus_coverage,
            completeness=list(completeness),
        ),
        aggregation=AggregationConfig(bin_size=bin_size, statistic=statistic),
    )
    config.output.dir = str(out_dir)
    if annotate or user_features:
        config.annotations.enabled = True
        config.annotations.provider = "dfam" if annotate else "none"
        if user_features:
            config.annotations.user_features_file = str(user_features)
    if save_config:
        save_cfg(config, save_config)
        click.echo(f"config written to {save_config}", err=True)
    _execute(config, Path(out_dir), no_figures, cache_dir)


@cli.command("gui")
def gui_cmd():
    """Launch the TEscope desktop GUI (requires the [gui] extra)."""
    try:
        from tescope.gui.app import main as gui_main
    except ImportError:
        raise click.ClickException(
            "The GUI requires PySide6. Install with: pip install 'tescope[gui]'"
        ) from None
    sys.exit(gui_main())


def main():  # pragma: no cover - console entry
    cli()


if __name__ == "__main__":  # pragma: no cover
    main()
