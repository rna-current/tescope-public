#!/usr/bin/env python
"""Generate THIRD-PARTY-NOTICES.txt for binary bundles.

Collects the license metadata and full license texts of every distribution
installed in the build environment (which is what PyInstaller freezes into
the bundle) so LGPL/BSD/MIT notice obligations are met in one file.

Usage (after `uv sync --extra gui`, before/after PyInstaller):
    python packaging/collect_licenses.py dist/TEscope/THIRD-PARTY-NOTICES.txt
"""

from __future__ import annotations

import sys
from importlib.metadata import distributions
from pathlib import Path

EXCLUDE = {
    # build/dev-only tools, not frozen into the application bundle
    "pyinstaller",
    "pyinstaller-hooks-contrib",
    "altgraph",
    "pytest",
    "pytest-cov",
    "coverage",
    "ruff",
    "pluggy",
    "iniconfig",
}

HEADER = """\
TEscope third-party notices
===========================

TEscope itself is distributed under the MIT License (see LICENSE).

This binary distribution bundles the following third-party Python
packages and shared libraries. PySide6/Shiboken6 (Qt for Python) are
used under the terms of the LGPL-3.0 license via dynamic linking; the
complete corresponding license texts are reproduced below. Reference
annotation data is NOT bundled — it is downloaded by each user directly
from its public source (see the documentation, "Reference data").

The PyInstaller bootloader embedded in the executable is licensed under
GPL 2.0 with the Bootloader Exception, which permits distributing the
bundled application under its own license.
"""


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("THIRD-PARTY-NOTICES.txt")
    parts = [HEADER]
    for d in sorted(distributions(), key=lambda d: d.metadata["Name"].lower()):
        name = d.metadata["Name"]
        if name.lower() in EXCLUDE or name.lower() == "tescope":
            continue
        lic = d.metadata.get("License-Expression") or d.metadata.get("License") or ""
        classifiers = [
            v for k, v in d.metadata.items() if k == "Classifier" and v.startswith("License")
        ]
        if classifiers and (not lic or len(lic) > 100):
            lic = "; ".join(c.split("::")[-1].strip() for c in classifiers)
        parts.append(f"\n{'=' * 72}\n{name} {d.version}  —  {lic.splitlines()[0] if lic else 'see text below'}\n{'=' * 72}\n")
        texts = []
        if d.files:
            for f in d.files:
                fname = f.name.lower()
                if "license" in fname or "copying" in fname or fname == "notice":
                    try:
                        texts.append(Path(d.locate_file(f)).read_text(errors="replace"))
                    except OSError:
                        pass
        if texts:
            parts.append("\n\n".join(texts))
        elif lic and len(lic) > 100:
            parts.append(lic)
        else:
            parts.append(f"(license: {lic or 'see project metadata'})")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(parts), encoding="utf-8")
    n = sum(1 for p in parts if p.startswith("\n" + "=" * 72))
    print(f"wrote {out} ({n} packages)")


if __name__ == "__main__":
    main()
