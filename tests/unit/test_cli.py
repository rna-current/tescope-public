"""CLI: headless operation of every command against the golden dataset."""

from __future__ import annotations

from pathlib import Path

import golden_dataset as gd
import pytest
from click.testing import CliRunner

from tescope import __version__
from tescope.cli.main import cli


@pytest.fixture()
def runner():
    return CliRunner()


@pytest.fixture(scope="module")
def cache_dir(golden_dir, tmp_path_factory):
    """A cache with the golden assembly installed via the CLI itself."""
    cache = tmp_path_factory.mktemp("clicache")
    r = CliRunner().invoke(
        cli,
        [
            "reference", "add",
            "--assembly", "testasm",
            "--repeatmasker", golden_dir["out"],
            "--chrom-sizes", golden_dir["chrom_sizes"],
            "--cache-dir", str(cache),
        ],
    )
    assert r.exit_code == 0, r.output
    return str(cache)


class TestBasics:
    def test_version(self, runner):
        r = runner.invoke(cli, ["--version"])
        assert r.exit_code == 0
        assert __version__ in r.output

    def test_help(self, runner):
        r = runner.invoke(cli, ["--help"])
        assert r.exit_code == 0
        for cmd in ("reference", "te", "run", "profile", "gui"):
            assert cmd in r.output


class TestReference:
    def test_add_output(self, cache_dir):
        # the fixture already ran `reference add`; verify idempotence guard
        r = CliRunner().invoke(
            cli,
            ["reference", "add", "--assembly", "testasm",
             "--repeatmasker", "/nonexistent", "--cache-dir", cache_dir],
        )
        assert r.exit_code != 0

    def test_list(self, runner, cache_dir):
        r = runner.invoke(cli, ["reference", "list", "--cache-dir", cache_dir])
        assert r.exit_code == 0
        assert "hg38" in r.output
        assert "testasm" in r.output
        assert "custom assembly" in r.output

    def test_fetch_unknown(self, runner, tmp_path):
        r = runner.invoke(cli, ["reference", "fetch", "nosuch", "--cache-dir", str(tmp_path)])
        assert r.exit_code != 0
        assert "not a known assembly" in r.output


class TestTe:
    def test_list(self, runner, cache_dir):
        r = runner.invoke(cli, ["te", "list", "testasm", "--cache-dir", cache_dir])
        assert r.exit_code == 0
        assert "TST1" in r.output and "TST2" in r.output

    def test_list_search(self, runner, cache_dir):
        r = runner.invoke(
            cli, ["te", "list", "testasm", "--search", "TST2", "--cache-dir", cache_dir]
        )
        assert "TST2" in r.output and "TST1 " not in r.output

    def test_info(self, runner, cache_dir):
        r = runner.invoke(cli, ["te", "info", "testasm", "TST1", "--cache-dir", cache_dir])
        assert r.exit_code == 0
        assert "consensus length:  1,000 bp" in r.output
        assert "full_length" in r.output

    def test_info_missing_repeat(self, runner, cache_dir):
        r = runner.invoke(cli, ["te", "info", "testasm", "TST9", "--cache-dir", cache_dir])
        assert r.exit_code != 0
        assert "not found" in r.output

    def test_not_installed_error_is_actionable(self, runner, tmp_path):
        r = runner.invoke(cli, ["te", "list", "hg38", "--cache-dir", str(tmp_path)])
        assert r.exit_code != 0
        assert "reference fetch hg38" in r.output


class TestProfileCommand:
    def test_full_profile_run(self, runner, cache_dir, golden_dir, tmp_path):
        out = tmp_path / "results"
        r = runner.invoke(
            cli,
            [
                "profile",
                "--genome", "testasm",
                "--te", "TST1",
                "--signal", golden_dir["bigwig"],
                "--peaks", golden_dir["peaks_bed"],
                "--bin-size", "100",
                "--min-length", "60",
                "--out", str(out),
                "--cache-dir", cache_dir,
                "--save-config", str(tmp_path / "saved.yaml"),
            ],
        )
        assert r.exit_code == 0, r.output
        assert (out / "qc.tsv").exists()
        assert (out / "analysis.yaml").exists()
        assert (out / "profile_stacked.png").exists()
        assert (tmp_path / "saved.yaml").exists()
        assert f"{len(gd.COPIES) - 1}/{len(gd.COPIES)} copies analyzed" in r.output

    def test_requires_a_track(self, runner, cache_dir):
        r = runner.invoke(cli, ["profile", "--genome", "testasm", "--te", "TST1",
                                "--cache-dir", cache_dir])
        assert r.exit_code != 0
        assert "at least one" in r.output

    def test_filter_error_message(self, runner, cache_dir, golden_dir, tmp_path):
        r = runner.invoke(
            cli,
            ["profile", "--genome", "testasm", "--te", "TST1",
             "--signal", golden_dir["bigwig"],
             "--min-length", "99999",
             "--out", str(tmp_path / "x"), "--cache-dir", cache_dir],
        )
        assert r.exit_code != 0
        assert "No TE annotations satisfy" in r.output


class TestRunCommand:
    def test_run_saved_config(self, runner, cache_dir, golden_dir, tmp_path):
        # Milestone H: a saved config runs headlessly through the CLI
        cfg = tmp_path / "analysis.yaml"
        cfg.write_text(
            f"""
assembly: testasm
repeat:
  name: TST1
tracks:
  - file: {golden_dir['bigwig']}
    name: signal
    group: ctrl
    replicate: 1
  - file: {golden_dir['peaks_bed']}
    name: peaks
aggregation:
  bin_size: 100
  statistic: mean
filters:
  min_consensus_coverage: 0.5
"""
        )
        out = tmp_path / "results"
        r = runner.invoke(
            cli, ["run", str(cfg), "--out", str(out), "--cache-dir", cache_dir]
        )
        assert r.exit_code == 0, r.output
        assert "2/7 copies analyzed" in r.output
        assert (out / "profile_signal.tsv").exists()

    def test_run_bad_config(self, runner, tmp_path):
        cfg = tmp_path / "bad.yaml"
        cfg.write_text("assembly: hg38\n")  # missing repeat/tracks
        r = runner.invoke(cli, ["run", str(cfg)])
        assert r.exit_code != 0
        assert "repeat" in r.output

    def test_relative_track_paths_resolve(self, runner, cache_dir, golden_dir, tmp_path):
        # track paths are relative to the config file's directory
        import shutil

        shutil.copy(golden_dir["bedgraph"], tmp_path / "sig.bedGraph")
        cfg = tmp_path / "analysis.yaml"
        cfg.write_text(
            "assembly: testasm\nrepeat:\n  name: TST1\n"
            "tracks:\n  - file: sig.bedGraph\naggregation:\n  bin_size: 100\n"
        )
        out = tmp_path / "res"
        r = runner.invoke(
            cli, ["run", str(cfg), "--out", str(out), "--no-figures",
                  "--cache-dir", cache_dir]
        )
        assert r.exit_code == 0, r.output
        assert Path(out, "qc.tsv").exists()
