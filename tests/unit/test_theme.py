"""Theme system smoke tests (offscreen; no pixel snapshots by design)."""

from __future__ import annotations

import os

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402

from tescope.gui.theme import (  # noqa: E402
    DEFAULT_THEME,
    THEMES,
    app_icon,
    apply_theme,
    kind_icon,
    resolve_theme,
    tokens_for,
)
from tescope.gui.theme.theme import DARK, LIGHT, build_qss  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


class TestTokens:
    def test_light_and_dark_share_structure(self):
        assert set(LIGHT) == set(DARK)

    def test_default_is_light(self):
        assert DEFAULT_THEME == "light"
        assert resolve_theme("nonsense") == "light"

    def test_system_resolves_to_concrete(self, qapp):
        assert resolve_theme("system") in ("light", "dark")


class TestQss:
    @pytest.mark.parametrize("tokens", [LIGHT, DARK], ids=["light", "dark"])
    def test_builds_with_no_unresolved_placeholders(self, tokens):
        qss = build_qss(tokens)
        assert len(qss) > 1000
        assert "{t[" not in qss  # no unrendered template fields
        # key hooks present
        for hook in ("#primaryButton", "#panelCard", "#statusMessage", "#canvasFrame"):
            assert hook in qss

    @pytest.mark.parametrize("name", ["light", "dark", "system"])
    def test_apply_theme(self, qapp, name):
        resolved = apply_theme(qapp, name)
        assert resolved in ("light", "dark")
        assert qapp.styleSheet()
        assert qapp.palette().highlight().color().name().lower() == \
            tokens_for(resolved)["accent"].lower()
        # restore default for other tests
        apply_theme(qapp, DEFAULT_THEME)


class TestIcons:
    def test_app_icon_renders(self, qapp):
        icon = app_icon()
        assert not icon.isNull()
        assert not icon.pixmap(32, 32).isNull()

    @pytest.mark.parametrize("kind", ["continuous", "interval", "point"])
    def test_kind_icons(self, qapp, kind):
        assert not kind_icon(kind).isNull()

    def test_unknown_kind_is_null_not_error(self, qapp):
        assert kind_icon("mystery").isNull()


class TestThemedWindow:
    def test_window_constructs_under_each_theme(self, qapp, golden_store):
        from tescope.gui.main_window import MainWindow

        for name in THEMES:
            apply_theme(qapp, name)
            w = MainWindow(store=golden_store, theme_pref=name)
            # live counts render with the current palette
            w.selection.genome.set_current("testasm")
            w.selection.te.select_repeat("TST1")
            assert "7" in w.selection.filters.counts.text()
            w.close()
        apply_theme(qapp, DEFAULT_THEME)

    def test_theme_switch_rerenders_counts(self, qapp, golden_store):
        from tescope.gui.main_window import MainWindow

        apply_theme(qapp, "light")
        w = MainWindow(store=golden_store, theme_pref="light")
        w.selection.genome.set_current("testasm")
        w.selection.te.select_repeat("TST1")
        light_html = w.selection.filters.counts.text()
        w._set_theme("dark")
        qapp.processEvents()
        dark_html = w.selection.filters.counts.text()
        assert "7" in dark_html
        assert dark_html != light_html  # palette-derived colors re-rendered
        w._set_theme("light")
        w.close()
        apply_theme(qapp, DEFAULT_THEME)

    def test_publication_style_is_theme_independent(self, qapp):
        """UI theme must not leak into exported scientific figures."""
        from tescope.plotting.style import COLOR_CYCLE, STYLE

        before = (dict(STYLE), list(COLOR_CYCLE))
        apply_theme(qapp, "dark")
        from tescope.plotting.style import COLOR_CYCLE as C2
        from tescope.plotting.style import STYLE as S2

        assert (dict(S2), list(C2)) == before
        apply_theme(qapp, DEFAULT_THEME)

    def test_publication_font_fallback_covers_non_latin(self, qapp):
        """R-001: non-ASCII track labels must not render as boxes when the
        host has a CJK-capable font — the publication style carries a
        glyph-fallback family list. Since the 2026-09 typography policy
        the declared primary face is Arial (-> Liberation Sans ->
        DejaVu Sans), with the CJK candidates host-filtered behind it."""
        from tescope.plotting.style import (
            _FALLBACK_CANDIDATES,
            PREFERRED_FAMILIES,
            STYLE,
        )

        # must be the family list itself: the "sans-serif" alias +
        # font.sans-serif does not build a glyph-fallback chain
        fams = STYLE["font.family"]
        assert isinstance(fams, list) and fams[0] == "Arial"
        assert list(PREFERRED_FAMILIES) == fams[: len(PREFERRED_FAMILIES)]
        assert "DejaVu Sans" in fams  # guaranteed glyph backstop
        # candidate set covers CJK; the optional tail is host-filtered so
        # absent families do not spam findfont warnings on every draw
        assert any("Noto Sans" in f or "CJK" in f for f in _FALLBACK_CANDIDATES)
        assert set(fams[len(PREFERRED_FAMILIES):-1]) <= set(_FALLBACK_CANDIDATES)
        assert fams[-1] == "sans-serif"
