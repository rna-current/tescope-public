"""TEscope theme system: design tokens -> QPalette + generated QSS.

Visual skin only — no scientific GUI logic lives here, and publication
figure styling (:mod:`tescope.plotting.style`) is deliberately independent
of the application theme.

Architecture:
- ``LIGHT`` / ``DARK``: flat token dictionaries (colors, spacing, sizing).
- ``build_qss(tokens)``: one QSS template rendered from tokens, so both
  themes stay structurally identical.
- ``apply_theme(app, name)``: sets the Fusion base style (consistent
  cross-platform rendering), a token-matched ``QPalette`` (drives the
  native Fusion checkbox/combo indicators, which we intentionally do not
  repaint with images), and the generated stylesheet.
- ``name`` may be ``"light"``, ``"dark"`` or ``"system"`` (resolved via
  Qt's color-scheme hint, falling back to light).

Widget hooks used by the QSS (set in widget code, semantics unchanged):
- ``objectName``: ``primaryButton`` (Run), ``destructiveButton`` (Remove),
  ``headerBar``/``headerTitle``/``headerSubtitle``, ``panelCard`` (section
  group boxes), ``hintLabel`` (secondary hints), ``canvasFrame`` (plot
  container), ``statusMessage``.
- dynamic property ``state`` on the status label:
  ``running`` / ``ok`` / ``warn`` / ``error``.
"""

from __future__ import annotations

from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication, QStyleFactory

THEMES = ("light", "dark", "system")
DEFAULT_THEME = "light"

# --------------------------------------------------------------------------
# design tokens
# --------------------------------------------------------------------------

LIGHT: dict[str, str] = {
    # surfaces
    "bg": "#F3F2EF",            # warm near-white application background
    "surface": "#FCFCFB",       # cards / panels
    "surface_alt": "#FFFFFF",   # inputs, item views
    "surface_sunken": "#ECEAE6",# toolbars, header strip
    # strokes
    "border": "#DCDAD4",
    "border_strong": "#C3C0B8",
    "separator": "#E6E4DF",
    # text
    "text": "#22272B",
    "text_secondary": "#6A7076",
    "text_disabled": "#A9ADB1",
    "text_on_accent": "#FFFFFF",
    # accent (restrained scientific teal-blue)
    "accent": "#166E8F",
    "accent_hover": "#125C78",
    "accent_pressed": "#0E4A61",
    "accent_soft": "#E2EEF3",   # selection / active-row wash
    "accent_soft_border": "#B9D5E0",
    "focus": "#3C8DAD",
    # semantic
    "danger": "#A8423C",
    "danger_soft": "#F5E4E2",
    "success": "#2E7D5B",
    "warning": "#9A6A1F",
    # controls
    "control_bg": "#FFFFFF",
    "control_bg_hover": "#F4F6F7",
    "control_bg_pressed": "#E9EDEF",
    "control_bg_disabled": "#F1F0ED",
    # misc
    "scroll_handle": "#C8C5BE",
    "scroll_handle_hover": "#AFACA4",
    "table_header": "#F2F1EE",
    "radius": "4px",
    "control_h": "25px",
}

DARK: dict[str, str] = {
    "bg": "#1F2225",
    "surface": "#26292D",
    "surface_alt": "#2B2F33",
    "surface_sunken": "#1A1C1F",
    "border": "#393E43",
    "border_strong": "#4A5056",
    "separator": "#33373C",
    "text": "#E7E5E1",
    "text_secondary": "#9BA1A7",
    "text_disabled": "#62676D",
    "text_on_accent": "#FFFFFF",
    "accent": "#4C9FBF",
    "accent_hover": "#5FAECB",
    "accent_pressed": "#3E8AA8",
    "accent_soft": "#243B45",
    "accent_soft_border": "#33586A",
    "focus": "#5FAECB",
    "danger": "#C9706A",
    "danger_soft": "#3A2A29",
    "success": "#5CA98A",
    "warning": "#C79A55",
    "control_bg": "#2E3236",
    "control_bg_hover": "#353A3F",
    "control_bg_pressed": "#3C4147",
    "control_bg_disabled": "#26292C",
    "scroll_handle": "#4A4F55",
    "scroll_handle_hover": "#5C6269",
    "table_header": "#2B2E32",
    "radius": "4px",
    "control_h": "25px",
}


def resolve_theme(name: str) -> str:
    """Resolve 'system' to 'light'/'dark' using Qt's color-scheme hint."""
    if name not in THEMES:
        name = DEFAULT_THEME
    if name != "system":
        return name
    app = QApplication.instance()
    try:
        from PySide6.QtCore import Qt

        scheme = app.styleHints().colorScheme() if app else None
        if scheme == Qt.ColorScheme.Dark:
            return "dark"
    except (AttributeError, ImportError):  # pragma: no cover - old Qt
        pass
    return "light"


def tokens_for(name: str) -> dict[str, str]:
    return DARK if resolve_theme(name) == "dark" else LIGHT


# --------------------------------------------------------------------------
# palette (drives Fusion-native indicators, menus, item views)
# --------------------------------------------------------------------------


def build_palette(t: dict[str, str]) -> QPalette:
    p = QPalette()
    c = QColor
    for group in (QPalette.Active, QPalette.Inactive):
        p.setColor(group, QPalette.Window, c(t["bg"]))
        p.setColor(group, QPalette.WindowText, c(t["text"]))
        p.setColor(group, QPalette.Base, c(t["surface_alt"]))
        p.setColor(group, QPalette.AlternateBase, c(t["surface"]))
        p.setColor(group, QPalette.Text, c(t["text"]))
        p.setColor(group, QPalette.Button, c(t["control_bg"]))
        p.setColor(group, QPalette.ButtonText, c(t["text"]))
        p.setColor(group, QPalette.Highlight, c(t["accent"]))
        p.setColor(group, QPalette.HighlightedText, c(t["text_on_accent"]))
        p.setColor(group, QPalette.ToolTipBase, c(t["surface_alt"]))
        p.setColor(group, QPalette.ToolTipText, c(t["text"]))
        p.setColor(group, QPalette.PlaceholderText, c(t["text_secondary"]))
        p.setColor(group, QPalette.Link, c(t["accent"]))
    for role in (QPalette.WindowText, QPalette.Text, QPalette.ButtonText):
        p.setColor(QPalette.Disabled, role, c(t["text_disabled"]))
    p.setColor(QPalette.Disabled, QPalette.Base, c(t["control_bg_disabled"]))
    p.setColor(QPalette.Disabled, QPalette.Button, c(t["control_bg_disabled"]))
    p.setColor(QPalette.Disabled, QPalette.Highlight, c(t["border_strong"]))
    return p


# --------------------------------------------------------------------------
# stylesheet
# --------------------------------------------------------------------------


def build_qss(t: dict[str, str]) -> str:
    return f"""
/* ===== base ============================================================ */
QWidget {{
    color: {t["text"]};
    font-size: 10pt;
}}
QMainWindow, QDialog {{
    background: {t["bg"]};
}}
QToolTip {{
    background: {t["surface_alt"]};
    color: {t["text"]};
    border: 1px solid {t["border_strong"]};
    padding: 4px 7px;
    border-radius: {t["radius"]};
}}

/* ===== header bar ====================================================== */
#headerBar {{
    background: {t["surface_sunken"]};
    border-bottom: 1px solid {t["border"]};
}}
#headerTitle {{
    font-size: 13pt;
    font-weight: 600;
    letter-spacing: 0.5px;
    color: {t["text"]};
}}
#headerSubtitle {{
    font-size: 9pt;
    color: {t["text_secondary"]};
}}

/* ===== scroll containers =============================================== */
QScrollArea#leftScroll, QScrollArea#leftScroll > QWidget > QWidget {{
    background: transparent;
    border: none;
}}

/* ===== panels / cards ================================================== */
QGroupBox#panelCard {{
    background: {t["surface"]};
    border: 1px solid {t["border"]};
    border-radius: 6px;
    margin-top: 12px;
}}
QGroupBox#panelCard::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 2px;
    padding: 0px 4px;
    color: {t["text_secondary"]};
    font-size: 9pt;
    font-weight: 600;
    letter-spacing: 0.8px;
    text-transform: uppercase;
}}
#hintLabel {{
    color: {t["text_secondary"]};
    font-size: 9pt;
}}
#canvasFrame {{
    background: {t["surface"]};
    border: 1px solid {t["border"]};
    border-radius: 6px;
}}

/* ===== buttons ========================================================= */
QPushButton {{
    background: {t["control_bg"]};
    border: 1px solid {t["border_strong"]};
    border-radius: {t["radius"]};
    padding: 4px 14px;
    min-height: {t["control_h"]};
}}
QPushButton:hover {{ background: {t["control_bg_hover"]}; }}
QPushButton:pressed {{ background: {t["control_bg_pressed"]}; }}
QPushButton:disabled {{
    background: {t["control_bg_disabled"]};
    border-color: {t["border"]};
    color: {t["text_disabled"]};
}}
QPushButton:focus {{
    border: 1px solid {t["focus"]};
    outline: none;
}}
QPushButton#primaryButton {{
    background: {t["accent"]};
    border: 1px solid {t["accent"]};
    color: {t["text_on_accent"]};
    font-weight: 600;
    padding: 6px 18px;
}}
QPushButton#primaryButton:hover {{ background: {t["accent_hover"]}; border-color: {t["accent_hover"]}; }}
QPushButton#primaryButton:pressed {{ background: {t["accent_pressed"]}; border-color: {t["accent_pressed"]}; }}
QPushButton#primaryButton:disabled {{
    background: {t["control_bg_disabled"]};
    border-color: {t["border"]};
    color: {t["text_disabled"]};
}}
QPushButton#primaryButton:focus {{ border: 2px solid {t["focus"]}; }}
QPushButton#destructiveButton:hover {{
    background: {t["danger_soft"]};
    border-color: {t["danger"]};
    color: {t["danger"]};
}}

/* ===== inputs ========================================================== */
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{
    background: {t["control_bg"]};
    border: 1px solid {t["border_strong"]};
    border-radius: {t["radius"]};
    padding: 2px 8px;
    min-height: {t["control_h"]};
    selection-background-color: {t["accent"]};
    selection-color: {t["text_on_accent"]};
}}
QLineEdit:hover, QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover {{
    border-color: {t["scroll_handle_hover"]};
}}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
    border: 1px solid {t["focus"]};
}}
QLineEdit:disabled, QComboBox:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled {{
    background: {t["control_bg_disabled"]};
    color: {t["text_disabled"]};
    border-color: {t["border"]};
}}
QComboBox::drop-down {{
    border: none;
    width: 22px;
}}
QComboBox QAbstractItemView {{
    background: {t["surface_alt"]};
    border: 1px solid {t["border_strong"]};
    selection-background-color: {t["accent_soft"]};
    selection-color: {t["text"]};
    outline: none;
}}
QSpinBox::up-button, QDoubleSpinBox::up-button,
QSpinBox::down-button, QDoubleSpinBox::down-button {{
    width: 16px;
    border: none;
    background: transparent;
}}
QCheckBox {{ spacing: 7px; min-height: 22px; }}
QCheckBox:disabled {{ color: {t["text_disabled"]}; }}

/* ===== item views (TE tree, tables) ==================================== */
QTreeWidget, QTableWidget, QListWidget {{
    background: {t["surface_alt"]};
    alternate-background-color: {t["surface"]};
    border: 1px solid {t["border"]};
    border-radius: {t["radius"]};
    outline: none;
}}
QTreeWidget::item {{ min-height: 24px; }}
QTableWidget::item {{ padding: 2px 6px; }}
QTreeWidget::item:hover, QTableWidget::item:hover {{
    background: {t["control_bg_hover"]};
}}
QTreeWidget::item:selected, QTableWidget::item:selected {{
    background: {t["accent_soft"]};
    color: {t["text"]};
    border-left: 2px solid {t["accent"]};
}}
QHeaderView::section {{
    background: {t["table_header"]};
    color: {t["text_secondary"]};
    font-size: 9pt;
    font-weight: 600;
    border: none;
    border-bottom: 1px solid {t["border"]};
    border-right: 1px solid {t["separator"]};
    padding: 4px 8px;
}}
QTableCornerButton::section {{
    background: {t["table_header"]};
    border: none;
    border-bottom: 1px solid {t["border"]};
}}

/* ===== tabs ============================================================ */
QTabWidget::pane {{
    border: none;
    border-top: 1px solid {t["border"]};
    top: -1px;
}}
QTabBar::tab {{
    background: transparent;
    color: {t["text_secondary"]};
    border: none;
    border-bottom: 2px solid transparent;
    padding: 7px 16px 6px 16px;
    margin-right: 2px;
}}
QTabBar::tab:hover {{ color: {t["text"]}; }}
QTabBar::tab:selected {{
    color: {t["accent"]};
    font-weight: 600;
    border-bottom: 2px solid {t["accent"]};
}}
QTabBar::tab:focus {{ outline: none; }}

/* ===== menus / toolbar / status ======================================== */
QMenuBar {{
    background: {t["surface_sunken"]};
    border-bottom: 1px solid {t["border"]};
    padding: 1px 4px;
}}
QMenuBar::item {{
    background: transparent;
    padding: 4px 10px;
    border-radius: {t["radius"]};
}}
QMenuBar::item:selected {{ background: {t["accent_soft"]}; color: {t["text"]}; }}
QMenu {{
    background: {t["surface_alt"]};
    border: 1px solid {t["border_strong"]};
    padding: 4px;
}}
QMenu::item {{
    padding: 5px 26px 5px 22px;
    border-radius: {t["radius"]};
}}
QMenu::item:selected {{ background: {t["accent_soft"]}; color: {t["text"]}; }}
QMenu::item:disabled {{ color: {t["text_disabled"]}; }}
QMenu::separator {{
    height: 1px;
    background: {t["separator"]};
    margin: 4px 8px;
}}
QToolBar {{
    background: transparent;
    border: none;
    spacing: 2px;
    padding: 2px 0px;
}}
QToolButton {{
    background: transparent;
    border: 1px solid transparent;
    border-radius: {t["radius"]};
    padding: 3px;
}}
QToolButton:hover {{
    background: {t["control_bg_hover"]};
    border-color: {t["border"]};
}}
QToolButton:pressed, QToolButton:checked {{
    background: {t["accent_soft"]};
    border-color: {t["accent_soft_border"]};
}}
QStatusBar {{
    background: {t["surface_sunken"]};
    border-top: 1px solid {t["border"]};
    color: {t["text_secondary"]};
}}
QStatusBar::item {{ border: none; }}
#statusMessage {{ color: {t["text_secondary"]}; padding-left: 6px; }}
#statusMessage[state="running"] {{ color: {t["accent"]}; }}
#statusMessage[state="ok"] {{ color: {t["success"]}; }}
#statusMessage[state="warn"] {{ color: {t["warning"]}; }}
#statusMessage[state="error"] {{ color: {t["danger"]}; }}

/* ===== progress / splitter / scrollbars ================================ */
QProgressBar {{
    background: {t["control_bg_pressed"]};
    border: none;
    border-radius: 3px;
    max-height: 7px;
    text-align: center;
    color: transparent;
}}
QProgressBar::chunk {{
    background: {t["accent"]};
    border-radius: 3px;
}}
QSplitter::handle {{
    background: transparent;
    width: 7px;
}}
QSplitter::handle:hover {{ background: {t["accent_soft"]}; }}
QScrollBar:vertical {{
    background: transparent;
    width: 11px;
    margin: 2px;
}}
QScrollBar::handle:vertical {{
    background: {t["scroll_handle"]};
    border-radius: 3px;
    min-height: 30px;
}}
QScrollBar::handle:vertical:hover {{ background: {t["scroll_handle_hover"]}; }}
QScrollBar:horizontal {{
    background: transparent;
    height: 11px;
    margin: 2px;
}}
QScrollBar::handle:horizontal {{
    background: {t["scroll_handle"]};
    border-radius: 3px;
    min-width: 30px;
}}
QScrollBar::handle:horizontal:hover {{ background: {t["scroll_handle_hover"]}; }}
QScrollBar::add-line, QScrollBar::sub-line {{
    height: 0px; width: 0px;
}}
QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}
"""


# --------------------------------------------------------------------------
# application entry
# --------------------------------------------------------------------------


def apply_theme(app: QApplication, name: str = DEFAULT_THEME) -> str:
    """Apply the named theme; returns the resolved concrete theme name."""
    resolved = resolve_theme(name)
    t = tokens_for(resolved)
    if app.style().objectName().lower() != "fusion" and "Fusion" in QStyleFactory.keys():
        app.setStyle("Fusion")
    app.setPalette(build_palette(t))
    app.setStyleSheet(build_qss(t))
    return resolved


def repolish(widget) -> None:
    """Re-evaluate QSS after a dynamic property change (e.g. status state)."""
    style = widget.style()
    style.unpolish(widget)
    style.polish(widget)
