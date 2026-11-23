"""Capture the manuscript figure panels and the User Guide screenshots.

Real application, real display, real data; window size / splitter / tab
arrangement only — the data are untouched. Unedited QWidget.grab().

Usage (from examples/h3k9me3_l1hs, after fetch_data.sh):
    python figure/capture_figure_panels.py <out_dir>
Environment: TESCOPE_CACHE_DIR (optional) selects the reference cache.
"""

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

EXAMPLE = Path(__file__).resolve().parents[1]
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
os.chdir(EXAMPLE)  # configuration paths are relative to the example dir

from PySide6.QtWidgets import QApplication, QSplitter  # noqa: E402

app = QApplication(sys.argv)

from tescope import __version__  # noqa: E402
from tescope.config.io import load_config  # noqa: E402
from tescope.core.analysis import run_analysis  # noqa: E402
from tescope.core.reference.store import ReferenceStore  # noqa: E402
from tescope.gui.main_window import MainWindow  # noqa: E402
from tescope.gui.widgets.about_dialog import AboutDialog  # noqa: E402
from tescope.gui.widgets.inspector import LocusInspector  # noqa: E402
from tescope.gui.widgets.style_dialog import PlotStyleDialog  # noqa: E402

cache = os.environ.get("TESCOPE_CACHE_DIR")
store = ReferenceStore(cache_dir=cache) if cache else ReferenceStore()
cfg = load_config(EXAMPLE / "h3k9me3_l1hs_subset.yaml")
cfg.style.preset = "presentation"

win = MainWindow(store=store)
win.resize(1360, 880)
win.show()
app.processEvents()
win.apply_config(cfg)
app.processEvents()
split = win.findChild(QSplitter)
if split is not None:  # narrow the sidebar: display-only arrangement
    split.setSizes([330, 1030])
app.processEvents()

print("running analysis…", flush=True)
result = run_analysis(cfg, store=store)
win.results.show_result(result)
app.processEvents()

# (a) main window, Profile tab, fragments metagene view
view = win.results
view.setCurrentIndex(0)
view.mode_combo.setCurrentText("fragments")
app.processEvents()
win.grab().save(str(OUT / "fig1_a_main_fragments.png"))

# (b) Loci tab: table + applied subset comparison
loci_i = next(i for i in range(view.count()) if view.tabText(i) == "Loci")
view.setCurrentIndex(loci_i)
app.processEvents()
tab = view.loci_tab
tab._apply_subset(compare=True)
app.processEvents()
for ch in tab.findChildren(QSplitter):
    ch.setSizes([300, 330])
app.processEvents()
# The matplotlib toolbar shows a transient "(x, y) = …" readout whenever
# the desktop pointer happens to rest over the canvas; clear it so the
# capture does not depend on where the pointer is.
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT  # noqa: E402

for tb in win.findChildren(NavigationToolbar2QT):
    tb.set_message("")
view.grab().save(str(OUT / "fig1_b_loci_compare.png"))

# (c) Locus Inspector on the top subset locus
ids = tab.selected_locus_ids()
insp = LocusInspector(win)
insp.resize(540, 600)
insp.show_locus(result, ids, ids[0])
insp.show()
app.processEvents()
insp.grab().save(str(OUT / "fig1_c_inspector.png"))

# User Guide dialogs
about = AboutDialog(win)
about.show()
app.processEvents()
about.grab().save(str(OUT / "about_dialog.png"))
about.close()
style = PlotStyleDialog(cfg.style, [t.name for t in cfg.tracks], parent=win)
style.show()
app.processEvents()
style.grab().save(str(OUT / "style_dialog.png"))
style.close()

prov = {
    "generated": datetime.now(timezone.utc).isoformat(),
    "tescope_version": __version__,
    "display": "real X11 display; unedited QWidget.grab(); window/splitter/"
               "tab arrangement only (transient pointer-coordinate readout "
               "cleared) — no data alteration",
    "data": "ENCODE ENCSR000APE: ENCFF601JGK GRCh38 K562 H3K9me3 fold change "
            "over control (bigWig), ENCFF963GZJ replicated peaks; hg38 "
            "RepeatMasker L1HS, 1661 records",
    "config": "examples/h3k9me3_l1hs/h3k9me3_l1hs_subset.yaml",
    "style_preset": "presentation",
    "n_selected": int(result.n_selected),
    "panels": {"a": "fig1_a_main_fragments.png",
               "b": "fig1_b_loci_compare.png",
               "c": "fig1_c_inspector.png"},
}
(OUT / "fig1_captures.provenance.json").write_text(json.dumps(prov, indent=2))
print("captured:", sorted(p.name for p in OUT.glob("*.png")), flush=True)
app.quit()
