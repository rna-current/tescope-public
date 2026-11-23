"""Background workers: analysis and reference downloads run off the UI thread."""

from __future__ import annotations

import traceback

from PySide6.QtCore import QObject, QThread, Signal

from ..config.model import AnalysisConfig
from ..core.analysis import AnalysisResult, run_analysis
from ..core.reference.store import ReferenceStore


class AnalysisWorker(QObject):
    progress = Signal(str, int, int)  # message, step, total
    finished = Signal(object)  # AnalysisResult
    failed = Signal(str)

    def __init__(self, config: AnalysisConfig, store: ReferenceStore):
        super().__init__()
        self.config = config
        self.store = store

    def run(self) -> None:
        try:
            result: AnalysisResult = run_analysis(
                self.config, self.store, progress=self.progress.emit
            )
        except Exception as e:  # noqa: BLE001 - boundary: everything becomes a dialog
            detail = "" if _is_user_error(e) else f"\n\n{traceback.format_exc(limit=5)}"
            self.failed.emit(f"{e}{detail}")
            return
        self.finished.emit(result)


class FetchWorker(QObject):
    progress = Signal(str, int, int)  # stage, done, total (bytes)
    finished = Signal(str)  # assembly name
    failed = Signal(str)

    def __init__(self, store: ReferenceStore, assembly: str,
                 custom_path: str | None = None, chrom_sizes: str | None = None):
        super().__init__()
        self.store = store
        self.assembly = assembly
        self.custom_path = custom_path
        self.chrom_sizes = chrom_sizes

    def run(self) -> None:
        try:
            if self.custom_path:
                self.store.add_custom(
                    self.assembly,
                    self.custom_path,
                    chrom_sizes_path=self.chrom_sizes,
                    force=True,
                    progress=self.progress.emit,
                )
            else:
                self.store.fetch(self.assembly, progress=self.progress.emit)
        except Exception as e:  # noqa: BLE001
            self.failed.emit(str(e))
            return
        self.finished.emit(self.assembly)


class AnnotationLookupWorker(QObject):
    """Fetch/refresh database consensus features off the UI thread."""

    finished = Signal(object)  # FeatureSet
    failed = Signal(str)

    def __init__(self, ann_cfg, repeat_name: str, cons_length: int, cache_dir,
                 refresh: bool = False):
        super().__init__()
        self.ann_cfg = ann_cfg
        self.repeat_name = repeat_name
        self.cons_length = cons_length
        self.cache_dir = cache_dir
        self.refresh = refresh

    def run(self) -> None:
        try:
            from ..core.features.cache import drop_cached
            from ..core.features.resolve import resolve_features

            if self.refresh:
                drop_cached(self.repeat_name, self.cons_length, "dfam", self.cache_dir)
            fs = resolve_features(
                self.ann_cfg, self.repeat_name, self.cons_length,
                cache_dir=self.cache_dir, allow_network=True,
            )
        except Exception as e:  # noqa: BLE001 - boundary
            self.failed.emit(str(e))
            return
        self.finished.emit(fs)


class SequenceLookupWorker(QObject):
    """Resolve the consensus sequence off the UI thread (cache-first)."""

    finished = Signal(object)  # SequenceResolution
    failed = Signal(str)

    def __init__(self, seq_cfg, repeat_name: str, cons_length: int, cache_dir,
                 refresh: bool = False):
        super().__init__()
        self.seq_cfg = seq_cfg
        self.repeat_name = repeat_name
        self.cons_length = cons_length
        self.cache_dir = cache_dir
        self.refresh = refresh

    def run(self) -> None:
        try:
            from ..core.sequence.cache import drop_cached
            from ..core.sequence.resolve import resolve_sequence

            if self.refresh and self.seq_cfg.source == "dfam":
                drop_cached(self.repeat_name, self.cons_length, "dfam", self.cache_dir)
            res = resolve_sequence(
                self.seq_cfg, self.repeat_name, self.cons_length,
                cache_dir=self.cache_dir, allow_network=True,
            )
        except Exception as e:  # noqa: BLE001 - boundary
            self.failed.emit(str(e))
            return
        self.finished.emit(res)


def start_in_thread(worker: QObject, parent=None) -> QThread:
    """Move ``worker`` to a fresh QThread and start it; caller keeps refs."""
    thread = QThread(parent)
    worker.moveToThread(thread)
    thread.started.connect(worker.run)
    worker.finished.connect(thread.quit)
    worker.failed.connect(thread.quit)
    return thread


def _is_user_error(e: Exception) -> bool:
    from ..config.model import ConfigError
    from ..core.analysis import AnalysisError
    from ..core.chromnames import ChromMismatchError
    from ..core.reference.store import ReferenceError

    return isinstance(e, (AnalysisError, ReferenceError, ConfigError, ChromMismatchError))
