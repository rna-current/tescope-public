"""Compact About dialog: modest top logo, full-width text, and
resource buttons that are shown ONLY when their destination resolves
(no nonfunctional links, no developer-machine paths displayed).

Resource resolution is layered so installed builds keep working
offline: packaged resources (``tescope.resources`` — license,
citation, offline help; present in wheels AND in frozen bundles via
the PyInstaller ``datas`` entry) → files shipped next to the frozen
executable (THIRD-PARTY-NOTICES.txt) → source-checkout extras
(docs/).
"""

from __future__ import annotations

import sys
from importlib import resources as importlib_resources
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
)

from tescope import __version__

from ..theme.icons import branding_logo_pixmap

REPO_URL = "https://github.com/rna-current/tescope-public"  # pyproject [project.urls]


def _packaged(name: str) -> Path | None:
    """A packaged tescope.resources file as a real filesystem path
    (wheel installs and PyInstaller onedir bundles ship real files)."""
    try:
        p = importlib_resources.files("tescope.resources") / name
        with importlib_resources.as_file(p) as real:
            if real.exists():
                return Path(real)
    except Exception:
        pass
    froot = _frozen_root()
    if froot is not None:  # frozen fallback: datas land under _MEIPASS
        cand = froot / "tescope/resources" / name
        if cand.exists():
            return cand
    return None


def _source_root() -> Path | None:
    """Repo root when running from a source checkout; None when frozen
    or installed without the repository files."""
    root = Path(__file__).resolve().parents[4]
    return root if (root / "pyproject.toml").exists() else None


def _frozen_root() -> Path | None:
    if getattr(sys, "frozen", False):
        return Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
    return None


def _local_resources() -> dict[str, Path]:
    """Resource label -> existing file, resolved in layers (packaged →
    frozen-bundle extras → source checkout). Only existing entries are
    returned, so no button is ever a dead end."""
    out: dict[str, Path] = {}
    for label, name in (("License (MIT)", "LICENSE.txt"),
                        ("How to cite", "CITATION.cff"),
                        ("Help", "OFFLINE_HELP.md")):
        p = _packaged(name)
        if p is not None:
            out[label] = p
    froot = _frozen_root()
    if froot is not None:
        for cand in (froot / "THIRD-PARTY-NOTICES.txt",
                     froot.parent / "THIRD-PARTY-NOTICES.txt"):
            if cand.exists():
                out["Third-party notices"] = cand
                break
    root = _source_root()
    if root is not None:
        if (root / "docs/index.md").exists():
            out["Documentation"] = root / "docs/index.md"
        if (root / "packaging/dist/TEscope/THIRD-PARTY-NOTICES.txt").exists() \
                and "Third-party notices" not in out:
            out["Third-party notices"] = (
                root / "packaging/dist/TEscope/THIRD-PARTY-NOTICES.txt")
    return out


class AboutDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About TEscope")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 16, 20, 12)
        layout.setSpacing(10)

        logo = branding_logo_pixmap()
        if not logo.isNull():
            label = QLabel()
            pm = logo.scaledToWidth(230, Qt.SmoothTransformation)
            pm.setDevicePixelRatio(1.0)
            label.setPixmap(pm)
            label.setAlignment(Qt.AlignHCenter)
            layout.addWidget(label)

        text = QLabel(
            f"<b>TEscope {__version__}</b> — fragment-aware visualization "
            "of genomic signals on transposable element consensus "
            "coordinates.<br><br>"
            "Developer: Hongseok Ha<br>"
            "License: MIT &nbsp;·&nbsp; "
            f'Repository: <a href="{REPO_URL}">github.com/rna-current/'
            "tescope-public</a><br><br>"
            "Your data stay on this machine; network access is used only "
            "to download public reference annotation and — when you "
            "request them — public Dfam consensus features/sequences."
        )
        text.setWordWrap(True)
        text.setOpenExternalLinks(True)
        text.setTextInteractionFlags(Qt.TextBrowserInteraction)
        layout.addWidget(text)

        resources = _local_resources()
        if resources:
            row = QHBoxLayout()
            row.setSpacing(8)
            for label, path in resources.items():
                btn = QPushButton(label)
                btn.setAutoDefault(False)
                btn.clicked.connect(
                    lambda _=False, p=path:
                    QDesktopServices.openUrl(QUrl.fromLocalFile(str(p))))
                row.addWidget(btn)
            row.addStretch(1)
            layout.addLayout(row)

        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.reject)
        buttons.accepted.connect(self.accept)
        layout.addWidget(buttons)
        self.setFixedWidth(430)
        self.adjustSize()
