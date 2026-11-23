"""TEscope GUI entry point (``tescope-gui`` / ``tescope gui``)."""

from __future__ import annotations

import sys


def main() -> int:
    try:
        from PySide6.QtCore import QSettings
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print(
            "TEscope GUI requires PySide6. Install with: pip install 'tescope[gui]'",
            file=sys.stderr,
        )
        return 1

    from .theme import DEFAULT_THEME, app_icon, apply_theme

    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("TEscope")
    app.setOrganizationName("TEscope")

    settings = QSettings("TEscope", "TEscope")
    theme_pref = str(settings.value("theme", DEFAULT_THEME))
    resolved = apply_theme(app, theme_pref)
    app.setWindowIcon(app_icon(dark=resolved == "dark"))

    from .main_window import MainWindow

    window = MainWindow(theme_pref=theme_pref)
    window.show()
    return app.exec()


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
