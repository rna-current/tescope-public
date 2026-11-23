"""Installed-build About/Help resources: packaged copies exist, stay in
sync with the repository originals, and the dialog's resolver finds
them without the source checkout."""

from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "src/tescope/resources"


def test_packaged_resources_exist_and_match_repo_originals():
    assert (RES / "OFFLINE_HELP.md").stat().st_size > 500
    assert (RES / "LICENSE.txt").read_bytes() == (ROOT / "LICENSE").read_bytes()
    assert (RES / "CITATION.cff").read_bytes() == (ROOT / "CITATION.cff").read_bytes()


def test_importlib_resolution():
    pytest.importorskip("PySide6")
    from tescope.gui.widgets.about_dialog import _packaged

    for name in ("LICENSE.txt", "CITATION.cff", "OFFLINE_HELP.md"):
        p = _packaged(name)
        assert p is not None and p.exists(), name


def test_resolver_offers_core_buttons():
    pytest.importorskip("PySide6")
    from tescope.gui.widgets.about_dialog import _local_resources

    labels = set(_local_resources())
    assert {"License (MIT)", "How to cite", "Help"} <= labels
