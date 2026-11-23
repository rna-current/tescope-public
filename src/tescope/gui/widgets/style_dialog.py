"""Figure-style dialog: edits the analysis' PlotStyleConfig in place.

Pure rendering controls (see plotting/plotstyle.py). The dialog mutates
the config's ``style`` section and asks the results view to redraw on
Apply — no analysis re-runs, no scientific values touched. Everything set
here serializes with the project config and reproduces headlessly.
"""

from __future__ import annotations

from matplotlib.colors import is_color_like, to_hex
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from ...plotting.plotstyle import (
    BANDS,
    FIGURE_WIDTH_PRESETS_CM,
    HEATMAP_CMAPS,
    HEATMAP_SORTS,
    LEGEND_POSITIONS,
    LINE_STYLES,
    PALETTES,
    PRESETS,
    REPRESENTATIONS,
    PlotStyleConfig,
    PlotStyleError,
    SeriesStyle,
)
from ...plotting.transforms import SMOOTHING_METHODS


def _optional_float(edit: QLineEdit, name: str):
    t = edit.text().strip()
    if not t:
        return None
    try:
        return float(t)
    except ValueError:
        raise PlotStyleError(f"{name}: not a number: {t!r}") from None


def _fmt_optional(v) -> str:
    return "" if v is None else f"{v:g}"


class _ColorButton(QPushButton):
    """Swatch button + picker
    empty value means 'palette default'."""

    def __init__(self, value: str = "", parent=None):
        super().__init__(parent)
        self.setFixedWidth(44)
        self._value = ""
        self.set_value(value)
        self.clicked.connect(self._pick)

    def set_value(self, value: str) -> None:
        self._value = value or ""
        if self._value:
            self.setStyleSheet(f"background-color: {to_hex(self._value)};")
            self.setText("")
        else:
            self.setStyleSheet("")
            self.setText("auto")

    def value(self) -> str:
        return self._value

    def _pick(self) -> None:
        initial = QColor(to_hex(self._value)) if self._value else QColor("#0072B2")
        c = QColorDialog.getColor(initial, self, "Series color")
        if c.isValid():
            self.set_value(c.name())


class PlotStyleDialog(QDialog):
    def __init__(self, style: PlotStyleConfig, series_labels: list[str],
                 on_apply=None, parent=None, input_caps: dict | None = None):
        super().__init__(parent)
        self.setWindowTitle("Figure style")
        self.setMinimumWidth(480)
        self._style = style
        self._on_apply = on_apply
        self._labels = list(series_labels)
        # capability summary of the CURRENT analysis inputs (results view):
        # has_continuous gates smoothing/representation; has_matrix gates
        # per-copy distribution bands (grouped-only analyses have none)
        self._caps = input_caps or {"has_continuous": True, "has_matrix": True}

        lay = QVBoxLayout(self)
        tabs = QTabWidget()
        tabs.addTab(self._general_tab(), "General")
        tabs.addTab(self._signal_tab(), "Signal")
        tabs.addTab(self._text_axes_tab(), "Text && axes")
        tabs.addTab(self._series_tab(), "Series")
        tabs.addTab(self._ann_heat_tab(), "Annotation && heatmap")
        lay.addWidget(tabs)

        btns = QDialogButtonBox(
            QDialogButtonBox.Apply | QDialogButtonBox.Reset | QDialogButtonBox.Close
        )
        btns.button(QDialogButtonBox.Apply).clicked.connect(self._apply)
        btns.button(QDialogButtonBox.Reset).setText("Reset to preset")
        btns.button(QDialogButtonBox.Reset).clicked.connect(self._reset)
        btns.rejected.connect(self.reject)
        lay.addWidget(btns)
        self._load()

    # -- tabs --------------------------------------------------------------

    def _general_tab(self) -> QWidget:
        w = QWidget()
        form = QFormLayout(w)
        self.preset_combo = QComboBox()
        self.preset_combo.addItems(list(PRESETS))
        form.addRow("Preset:", self.preset_combo)
        self.palette_combo = QComboBox()
        self.palette_combo.addItems(list(PALETTES))
        form.addRow("Palette:", self.palette_combo)

        row = QHBoxLayout()
        self.width_combo = QComboBox()
        self.width_combo.addItem("preset default", 0.0)
        for name, cm in FIGURE_WIDTH_PRESETS_CM.items():
            self.width_combo.addItem(f"{name} (~{cm} cm)", cm)
        self.width_combo.addItem("custom", -1.0)
        row.addWidget(self.width_combo, 1)
        self.width_spin = QDoubleSpinBox()
        self.width_spin.setRange(0.0, 60.0)
        self.width_spin.setSuffix(" cm")
        self.width_spin.setSpecialValueText("auto")
        row.addWidget(self.width_spin)
        form.addRow("Figure width:", row)
        self.width_combo.currentIndexChanged.connect(self._width_preset_changed)

        self.height_spin = QDoubleSpinBox()
        self.height_spin.setRange(0.0, 60.0)
        self.height_spin.setSuffix(" cm")
        self.height_spin.setSpecialValueText("auto (per panel)")
        form.addRow("Figure height:", self.height_spin)

        self.title_edit = QLineEdit()
        self.title_edit.setPlaceholderText("default: repeat, assembly, copy counts")
        form.addRow("Title:", self.title_edit)
        self.xlabel_edit = QLineEdit()
        self.xlabel_edit.setPlaceholderText("default: consensus position axis")
        form.addRow("X label:", self.xlabel_edit)
        self.ylabel_edit = QLineEdit()
        self.ylabel_edit.setPlaceholderText("default: statistic (overlay view)")
        form.addRow("Y label:", self.ylabel_edit)
        return w

    def _signal_tab(self) -> QWidget:
        w = QWidget()
        form = QFormLayout(w)
        self.repr_combo = QComboBox()
        self.repr_combo.addItems(list(REPRESENTATIONS))
        self.repr_combo.setToolTip(
            "'filled' adds a translucent area under each profile line; the "
            "plotted values are identical.")
        form.addRow("Representation:", self.repr_combo)
        self.fill_alpha = QDoubleSpinBox()
        self.fill_alpha.setRange(0.0, 1.0)
        self.fill_alpha.setSingleStep(0.05)
        form.addRow("Fill alpha:", self.fill_alpha)

        self.smooth_combo = QComboBox()
        self.smooth_combo.addItems(list(SMOOTHING_METHODS))
        self.smooth_combo.setToolTip(
            "DISPLAY-ONLY smoothing: exported profile tables keep the "
            "unsmoothed values; smoothed values export separately as "
            "profile_smoothed_*.tsv with method and parameters.")
        form.addRow("Smoothing:", self.smooth_combo)
        self.smooth_window = QDoubleSpinBox()
        self.smooth_window.setRange(3, 201)
        self.smooth_window.setDecimals(0)
        self.smooth_window.setSingleStep(2)
        form.addRow("Window (bins, odd):", self.smooth_window)
        self.smooth_sigma = QDoubleSpinBox()
        self.smooth_sigma.setRange(0.1, 50.0)
        self.smooth_sigma.setSingleStep(0.5)
        form.addRow("Gaussian sigma (bins):", self.smooth_sigma)
        self.smooth_poly = QDoubleSpinBox()
        self.smooth_poly.setRange(1, 9)
        self.smooth_poly.setDecimals(0)
        form.addRow("Sav.-Golay order:", self.smooth_poly)
        self.show_raw_check = QCheckBox("Show raw profile under smoothed (thin line)")
        form.addRow(self.show_raw_check)

        self.band_combo = QComboBox()
        self.band_combo.addItems(list(BANDS))
        self.band_combo.setToolTip(
            "sem: uncertainty of the group/track mean (SEM). iqr: variation "
            "ACROSS TE COPIES (25–75% of per-copy values; ungrouped tracks "
            "only — grouped series have no single per-copy matrix). These "
            "are different quantities; see docs.")
        form.addRow("Band:", self.band_combo)
        self.outer_band_check = QCheckBox("Add 10–90% outer ribbon (iqr band)")
        form.addRow(self.outer_band_check)

        self.shading_check = QCheckBox("Shade annotated feature regions in signal panels")
        form.addRow(self.shading_check)
        self.fragment_check = QCheckBox("Also export fragment-comparison figure "
                                        "(all / full-length / 5'-truncated)")
        form.addRow(self.fragment_check)

        if not self._caps.get("has_continuous", True):
            for c in (self.repr_combo, self.fill_alpha, self.smooth_combo,
                      self.smooth_window, self.smooth_sigma, self.smooth_poly,
                      self.show_raw_check, self.band_combo, self.outer_band_check):
                c.setEnabled(False)
                c.setToolTip("Unavailable: the current analysis has no "
                             "continuous signal track (interval/point "
                             "occurrence output cannot be smoothed).")
        elif not self._caps.get("has_matrix", True):
            i = self.band_combo.findText("iqr")
            item = self.band_combo.model().item(i)
            item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
            item.setToolTip("Unavailable: all continuous series are grouped; "
                            "per-copy distribution bands need an ungrouped "
                            "track's per-copy matrix.")
        return w

    def _text_axes_tab(self) -> QWidget:
        w = QWidget()
        form = QFormLayout(w)
        self.font_spins = {}
        for key, label in (("title", "Title size"), ("label", "Axis-label size"),
                           ("tick", "Tick-label size"), ("legend", "Legend size"),
                           ("annotation", "Annotation-label size")):
            sp = QDoubleSpinBox()
            sp.setRange(0.0, 32.0)
            sp.setSuffix(" pt")
            sp.setSpecialValueText("preset")
            self.font_spins[key] = sp
            form.addRow(f"{label}:", sp)

        self.legend_check = QCheckBox("Show legend (overlay view)")
        form.addRow(self.legend_check)
        self.legend_combo = QComboBox()
        self.legend_combo.addItems(list(LEGEND_POSITIONS))
        form.addRow("Legend position:", self.legend_combo)

        self.grid_check = QCheckBox("Show grid")
        form.addRow(self.grid_check)
        self.grid_alpha = QDoubleSpinBox()
        self.grid_alpha.setRange(0.05, 1.0)
        self.grid_alpha.setSingleStep(0.05)
        form.addRow("Grid alpha:", self.grid_alpha)

        self.ymin_edit = QLineEdit()
        self.ymin_edit.setPlaceholderText("auto")
        self.ymax_edit = QLineEdit()
        self.ymax_edit.setPlaceholderText("auto")
        yrow = QHBoxLayout()
        yrow.addWidget(self.ymin_edit)
        yrow.addWidget(QLabel("to"))
        yrow.addWidget(self.ymax_edit)
        form.addRow("Y display range:", yrow)
        self.xmin_edit = QLineEdit()
        self.xmin_edit.setPlaceholderText("0")
        self.xmax_edit = QLineEdit()
        self.xmax_edit.setPlaceholderText("consensus end")
        xrow = QHBoxLayout()
        xrow.addWidget(self.xmin_edit)
        xrow.addWidget(QLabel("to"))
        xrow.addWidget(self.xmax_edit)
        form.addRow("X view window (bp):", xrow)
        note = QLabel("Display ranges are view windows only — no data are "
                      "filtered or transformed.")
        note.setObjectName("hintLabel")
        note.setWordWrap(True)
        form.addRow(note)
        return w

    def _series_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        self.series_table = QTableWidget(len(self._labels), 5)
        self.series_table.setHorizontalHeaderLabels(
            ["Series", "Color", "Hex", "Line style", "Width"])
        self.series_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.series_table.verticalHeader().setVisible(False)
        self._series_widgets = []
        for r, label in enumerate(self._labels):
            item = QTableWidgetItem(label)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            self.series_table.setItem(r, 0, item)
            btn = _ColorButton()
            hexedit = QLineEdit()
            hexedit.setPlaceholderText("auto")
            ls = QComboBox()
            ls.addItem("preset", "")
            [ls.addItem(x, x) for x in LINE_STYLES]
            width = QDoubleSpinBox()
            width.setRange(0.0, 8.0)
            width.setSingleStep(0.1)
            width.setSpecialValueText("preset")
            btn.clicked.connect(lambda _=False, b=btn, h=hexedit: h.setText(b.value()))
            self.series_table.setCellWidget(r, 1, btn)
            self.series_table.setCellWidget(r, 2, hexedit)
            self.series_table.setCellWidget(r, 3, ls)
            self.series_table.setCellWidget(r, 4, width)
            self._series_widgets.append((label, btn, hexedit, ls, width))
        lay.addWidget(self.series_table)
        hint = QLabel("Colors: pick or type any hex value, or leave empty "
                      "for the palette. Saved choices are never reassigned.")
        hint.setObjectName("hintLabel")
        hint.setWordWrap(True)
        lay.addWidget(hint)
        return w

    def _ann_heat_tab(self) -> QWidget:
        w = QWidget()
        form = QFormLayout(w)
        self.ann_alpha = QDoubleSpinBox()
        self.ann_alpha.setRange(0.1, 1.0)
        self.ann_alpha.setSingleStep(0.05)
        form.addRow("Annotation alpha:", self.ann_alpha)
        self.ann_border = QDoubleSpinBox()
        self.ann_border.setRange(0.0, 3.0)
        self.ann_border.setSingleStep(0.1)
        form.addRow("Annotation border width:", self.ann_border)
        self.db_lane_btn = _ColorButton()
        form.addRow("Database lane color:", self.db_lane_btn)
        self.user_lane_btn = _ColorButton()
        form.addRow("User lane color:", self.user_lane_btn)
        note = QLabel("Origin stays visible whatever the colors: the "
                      "user-defined lane keeps its hatch and both lanes keep "
                      "their labels.")
        note.setObjectName("hintLabel")
        note.setWordWrap(True)
        form.addRow(note)

        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems(list(HEATMAP_CMAPS))
        form.addRow("Heatmap colormap:", self.cmap_combo)
        self.hvmin_edit = QLineEdit()
        self.hvmin_edit.setPlaceholderText("auto")
        self.hvmax_edit = QLineEdit()
        self.hvmax_edit.setPlaceholderText("auto")
        hrow = QHBoxLayout()
        hrow.addWidget(self.hvmin_edit)
        hrow.addWidget(QLabel("to"))
        hrow.addWidget(self.hvmax_edit)
        form.addRow("Heatmap color limits:", hrow)
        self.cbar_check = QCheckBox("Show colorbar")
        form.addRow(self.cbar_check)
        self.hsort_combo = QComboBox()
        self.hsort_combo.addItems(list(HEATMAP_SORTS))
        self.hsort_combo.setToolTip("Row ORDER only — values never change.")
        form.addRow("Heatmap row sort:", self.hsort_combo)
        self.hprofile_check = QCheckBox("Aggregate profile panel above heatmap")
        form.addRow(self.hprofile_check)
        note2 = QLabel("Color limits change display only — the underlying "
                       "values are never rescaled (this is not normalization).")
        note2.setObjectName("hintLabel")
        note2.setWordWrap(True)
        form.addRow(note2)
        return w

    # -- load / collect ----------------------------------------------------

    def _width_preset_changed(self) -> None:
        v = self.width_combo.currentData()
        if v is not None and v >= 0:
            self.width_spin.setValue(float(v))

    def _load(self) -> None:
        st = self._style
        self.preset_combo.setCurrentText(st.preset)
        self.palette_combo.setCurrentText(st.palette)
        self.width_spin.setValue(st.fig_width_cm)
        self.height_spin.setValue(st.fig_height_cm)
        self.title_edit.setText(st.title)
        self.xlabel_edit.setText(st.x_label)
        self.ylabel_edit.setText(st.y_label)
        for key, sp in self.font_spins.items():
            sp.setValue(getattr(st, f"{key}_size"))
        self.legend_check.setChecked(st.show_legend)
        self.legend_combo.setCurrentText(st.legend_position)
        self.grid_check.setChecked(st.show_grid)
        self.grid_alpha.setValue(st.grid_alpha)
        self.repr_combo.setCurrentText(st.representation)
        self.fill_alpha.setValue(st.fill_alpha)
        self.smooth_combo.setCurrentText(st.smoothing_method)
        self.smooth_window.setValue(st.smoothing_window)
        self.smooth_sigma.setValue(st.smoothing_sigma)
        self.smooth_poly.setValue(st.smoothing_polyorder)
        self.show_raw_check.setChecked(st.smoothing_show_raw)
        self.band_combo.setCurrentText(st.band)
        self.outer_band_check.setChecked(st.outer_band)
        self.shading_check.setChecked(st.annotation_shading)
        self.fragment_check.setChecked(st.fragment_panels)
        self.hsort_combo.setCurrentText(st.heatmap_sort)
        self.hprofile_check.setChecked(st.heatmap_show_profile)
        self.ymin_edit.setText(_fmt_optional(st.y_min))
        self.ymax_edit.setText(_fmt_optional(st.y_max))
        self.xmin_edit.setText(_fmt_optional(st.x_min))
        self.xmax_edit.setText(_fmt_optional(st.x_max))
        for label, btn, hexedit, ls, width in self._series_widgets:
            s = st.series_style(label)
            btn.set_value(s.color)
            hexedit.setText(s.color)
            ls.setCurrentIndex(max(0, ls.findData(s.line_style)))
            width.setValue(s.line_width)
        self.ann_alpha.setValue(st.annotation_alpha)
        self.ann_border.setValue(st.annotation_border_width)
        self.db_lane_btn.set_value(st.database_lane_color)
        self.user_lane_btn.set_value(st.user_lane_color)
        self.cmap_combo.setCurrentText(st.heatmap_cmap)
        self.hvmin_edit.setText(_fmt_optional(st.heatmap_vmin))
        self.hvmax_edit.setText(_fmt_optional(st.heatmap_vmax))
        self.cbar_check.setChecked(st.heatmap_colorbar)

    def _collect(self) -> None:
        st = self._style
        st.preset = self.preset_combo.currentText()
        st.palette = self.palette_combo.currentText()
        st.fig_width_cm = self.width_spin.value()
        st.fig_height_cm = self.height_spin.value()
        st.title = self.title_edit.text()
        st.x_label = self.xlabel_edit.text()
        st.y_label = self.ylabel_edit.text()
        for key, sp in self.font_spins.items():
            setattr(st, f"{key}_size", sp.value())
        st.show_legend = self.legend_check.isChecked()
        st.legend_position = self.legend_combo.currentText()
        st.show_grid = self.grid_check.isChecked()
        st.grid_alpha = self.grid_alpha.value()
        st.representation = self.repr_combo.currentText()
        st.fill_alpha = self.fill_alpha.value()
        st.smoothing_method = self.smooth_combo.currentText()
        st.smoothing_window = int(self.smooth_window.value())
        st.smoothing_sigma = self.smooth_sigma.value()
        st.smoothing_polyorder = int(self.smooth_poly.value())
        st.smoothing_show_raw = self.show_raw_check.isChecked()
        st.band = self.band_combo.currentText()
        st.outer_band = self.outer_band_check.isChecked()
        st.annotation_shading = self.shading_check.isChecked()
        st.fragment_panels = self.fragment_check.isChecked()
        st.heatmap_sort = self.hsort_combo.currentText()
        st.heatmap_show_profile = self.hprofile_check.isChecked()
        st.y_min = _optional_float(self.ymin_edit, "Y min")
        st.y_max = _optional_float(self.ymax_edit, "Y max")
        st.x_min = _optional_float(self.xmin_edit, "X min")
        st.x_max = _optional_float(self.xmax_edit, "X max")
        series: dict = dict(st.series or {})
        for label, btn, hexedit, ls, width in self._series_widgets:
            color = hexedit.text().strip() or btn.value()
            if color and not is_color_like(color):
                raise PlotStyleError(f"series '{label}': invalid color {color!r}")
            entry = SeriesStyle(
                color=to_hex(color) if color else "",
                line_style=ls.currentData() or "",
                line_width=width.value(),
            )
            if entry.color or entry.line_style or entry.line_width > 0:
                series[label] = entry.to_dict()
            else:
                series.pop(label, None)
        st.series = series
        st.annotation_alpha = self.ann_alpha.value()
        st.annotation_border_width = self.ann_border.value()
        st.database_lane_color = self.db_lane_btn.value()
        st.user_lane_color = self.user_lane_btn.value()
        st.heatmap_cmap = self.cmap_combo.currentText()
        st.heatmap_vmin = _optional_float(self.hvmin_edit, "heatmap min")
        st.heatmap_vmax = _optional_float(self.hvmax_edit, "heatmap max")
        st.heatmap_colorbar = self.cbar_check.isChecked()
        st.validate()

    # -- actions -----------------------------------------------------------

    def _apply(self) -> None:
        try:
            self._collect()
        except PlotStyleError as e:
            QMessageBox.warning(self, "Invalid style", str(e))
            return
        if self._on_apply:
            self._on_apply()

    def _reset(self) -> None:
        preset = self.preset_combo.currentText()
        fresh = PlotStyleConfig.preset_defaults(preset)
        for f in fresh.__dataclass_fields__:
            setattr(self._style, f, getattr(fresh, f))
        self._load()
        if self._on_apply:
            self._on_apply()
