"""Local reference store: download, import, cache and query TE annotations.

Layout (under the TEscope cache directory):

```
<cache>/references/<assembly>/
├── source/<original file>      retained original annotation (provenance)
├── annotation.npz              parsed canonical frame (frame_cache format)
├── hierarchy.npz               class/family/name copy-count table
├── cons_lengths.tsv            repeat name -> consensus model length
├── chrom.sizes                 reference chromosome sizes (when available)
└── provenance.json             source URL/path, retrieval date, md5, versions
```

The cache directory resolves from, in order: an explicit argument, the
``TESCOPE_CACHE_DIR`` environment variable, the platform user-cache dir
(platformdirs). User signal data is never written here — only public or
user-supplied reference annotation.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import urllib.error
import urllib.request
from collections.abc import Callable
from pathlib import Path

import pandas as pd

from ..annotation.frame_cache import FrameCacheError, load_frame, save_frame
from ..annotation.repeatmasker import parse_any
from ..annotation.schema import (
    AnnotationProvenance,
    consensus_lengths,
    repeat_hierarchy,
)
from .registry import KNOWN_ASSEMBLIES, AssemblySpec

ProgressCallback = Callable[[str, int, int], None]  # (stage, done_bytes, total_bytes)


class ReferenceError(RuntimeError):
    pass


class ReferenceNotInstalledError(ReferenceError):
    pass


def default_cache_dir() -> Path:
    env = os.environ.get("TESCOPE_CACHE_DIR")
    if env:
        return Path(env)
    import platformdirs

    return Path(platformdirs.user_cache_dir("tescope"))


def _md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _download(url: str, dest: Path, progress: ProgressCallback | None = None,
              retries: int = 3, timeout: int = 60) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err: Exception | None = None
    for _attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "TEscope"})
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as out:
                total = int(resp.headers.get("Content-Length") or 0)
                done = 0
                while True:
                    block = resp.read(1 << 20)
                    if not block:
                        break
                    out.write(block)
                    done += len(block)
                    if progress:
                        progress("download", done, total)
            os.replace(tmp, dest)
            return
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last_err = e
            tmp.unlink(missing_ok=True)
    raise ReferenceError(f"download failed after {retries} attempts: {url} ({last_err})")


class ReferenceStore:
    def __init__(self, cache_dir: str | os.PathLike | None = None):
        self.cache_dir = Path(cache_dir) if cache_dir else default_cache_dir()
        self.refs_dir = self.cache_dir / "references"
        self._annotation_cache: dict[str, pd.DataFrame] = {}

    # -- discovery ---------------------------------------------------------

    def assembly_dir(self, assembly: str) -> Path:
        return self.refs_dir / assembly

    def is_installed(self, assembly: str) -> bool:
        d = self.assembly_dir(assembly)
        return (d / "annotation.npz").exists() and (d / "provenance.json").exists()

    def installed_assemblies(self) -> list[str]:
        if not self.refs_dir.exists():
            return []
        return sorted(d.name for d in self.refs_dir.iterdir() if self.is_installed(d.name))

    def known_assemblies(self) -> dict[str, AssemblySpec]:
        return dict(KNOWN_ASSEMBLIES)

    def provenance(self, assembly: str) -> AnnotationProvenance:
        p = self.assembly_dir(assembly) / "provenance.json"
        if not p.exists():
            raise ReferenceNotInstalledError(self._not_installed_msg(assembly))
        return AnnotationProvenance.from_dict(json.loads(p.read_text()))

    def _not_installed_msg(self, assembly: str) -> str:
        if assembly in KNOWN_ASSEMBLIES:
            return (
                f"Reference annotation for '{assembly}' is not yet downloaded. "
                f"Run: tescope reference fetch {assembly}"
            )
        return (
            f"Assembly '{assembly}' is not installed. Import a RepeatMasker file with: "
            f"tescope reference add --assembly {assembly} --repeatmasker <file.out>"
        )

    # -- installation ------------------------------------------------------

    def fetch(self, assembly: str, progress: ProgressCallback | None = None,
              force: bool = False) -> AnnotationProvenance:
        """Download and install a known public assembly annotation."""
        spec = KNOWN_ASSEMBLIES.get(assembly)
        if spec is None:
            raise ReferenceError(
                f"'{assembly}' is not a known assembly "
                f"({', '.join(sorted(KNOWN_ASSEMBLIES))}); use reference add for custom ones"
            )
        if self.is_installed(assembly) and not force:
            return self.provenance(assembly)
        d = self.assembly_dir(assembly)
        src_dir = d / "source"
        src_dir.mkdir(parents=True, exist_ok=True)
        src_file = src_dir / os.path.basename(spec.rmsk_url)
        if progress:
            progress("download", 0, 0)
        _download(spec.rmsk_url, src_file, progress)
        chrom_sizes = d / "chrom.sizes"
        try:
            _download(spec.chromsizes_url, chrom_sizes, progress)
        except ReferenceError:
            chrom_sizes = None  # non-fatal: chrom sizes are advisory
        prov = AnnotationProvenance(
            assembly=assembly,
            source=spec.source,
            source_url=spec.rmsk_url,
            source_path=str(src_file),
            retrieved_at=AnnotationProvenance.now(),
            md5=_md5(src_file),
        )
        self._install_from_source(assembly, src_file, prov, progress)
        return prov

    def add_custom(
        self,
        assembly: str,
        repeatmasker_path: str | os.PathLike,
        *,
        chrom_sizes_path: str | os.PathLike | None = None,
        notes: str = "",
        progress: ProgressCallback | None = None,
        force: bool = False,
    ) -> AnnotationProvenance:
        """Install a user-provided RepeatMasker annotation for any assembly."""
        if self.is_installed(assembly) and not force:
            raise ReferenceError(
                f"assembly '{assembly}' already installed; pass force=True to replace"
            )
        src = Path(repeatmasker_path)
        if not src.exists():
            raise FileNotFoundError(src)
        d = self.assembly_dir(assembly)
        src_dir = d / "source"
        src_dir.mkdir(parents=True, exist_ok=True)
        kept = src_dir / src.name
        if kept.resolve() != src.resolve():
            shutil.copy2(src, kept)
        if chrom_sizes_path:
            shutil.copy2(chrom_sizes_path, d / "chrom.sizes")
        prov = AnnotationProvenance(
            assembly=assembly,
            source="user-file",
            source_path=str(src),
            retrieved_at=AnnotationProvenance.now(),
            md5=_md5(kept),
            notes=notes,
        )
        self._install_from_source(assembly, kept, prov, progress)
        return prov

    def _install_from_source(
        self,
        assembly: str,
        source_file: Path,
        prov: AnnotationProvenance,
        progress: ProgressCallback | None = None,
    ) -> None:
        if progress:
            progress("parse", 0, 0)
        df, fmt = parse_any(source_file)
        prov.source_format = fmt
        prov.extra["n_hits"] = int(len(df))
        prov.extra["n_malformed_lines"] = int(df.attrs.get("n_malformed_lines", 0))
        prov.extra["tescope_version"] = _version_string()
        d = self.assembly_dir(assembly)
        d.mkdir(parents=True, exist_ok=True)
        if progress:
            progress("index", 0, 0)
        save_frame(df, d / "annotation.npz")
        save_frame(repeat_hierarchy(df), d / "hierarchy.npz")
        consensus_lengths(df).rename("cons_length").to_csv(
            d / "cons_lengths.tsv", sep="\t", header=True
        )
        (d / "provenance.json").write_text(json.dumps(prov.to_dict(), indent=2))
        self._annotation_cache.pop(assembly, None)

    def remove(self, assembly: str) -> None:
        d = self.assembly_dir(assembly)
        if d.exists():
            shutil.rmtree(d)
        self._annotation_cache.pop(assembly, None)

    # -- queries -----------------------------------------------------------

    def load_annotation(self, assembly: str) -> pd.DataFrame:
        if assembly in self._annotation_cache:
            return self._annotation_cache[assembly]
        d = self.assembly_dir(assembly)
        path = d / "annotation.npz"
        if not path.exists():
            raise ReferenceNotInstalledError(self._not_installed_msg(assembly))
        try:
            df = load_frame(path)
        except FrameCacheError:
            # cache format drift: rebuild from the retained source file
            prov = self.provenance(assembly)
            src = Path(prov.source_path)
            if not src.exists():
                src_dir = d / "source"
                candidates = list(src_dir.iterdir()) if src_dir.exists() else []
                if not candidates:
                    raise
                src = candidates[0]
            self._install_from_source(assembly, src, prov)
            df = load_frame(path)
        self._annotation_cache[assembly] = df
        return df

    def hierarchy(self, assembly: str) -> pd.DataFrame:
        d = self.assembly_dir(assembly)
        p = d / "hierarchy.npz"
        if p.exists():
            try:
                return load_frame(p)
            except FrameCacheError:
                pass
        return repeat_hierarchy(self.load_annotation(assembly))

    def consensus_lengths(self, assembly: str) -> pd.Series:
        p = self.assembly_dir(assembly) / "cons_lengths.tsv"
        if p.exists():
            t = pd.read_csv(p, sep="\t", index_col=0)
            return t["cons_length"]
        return consensus_lengths(self.load_annotation(assembly))

    def chrom_sizes(self, assembly: str) -> dict[str, int] | None:
        p = self.assembly_dir(assembly) / "chrom.sizes"
        if not p.exists():
            return None
        out: dict[str, int] = {}
        for line in p.read_text().splitlines():
            f = line.split()
            if len(f) >= 2:
                out[f[0]] = int(f[1])
        return out

    def select_repeat(self, assembly: str, repeat_name: str) -> pd.DataFrame:
        """All hits for one repeat name (subfamily), canonical frame subset."""
        df = self.load_annotation(assembly)
        sub = df[df["repeat_name"] == repeat_name]
        if len(sub) == 0:
            hier = self.hierarchy(assembly)
            names = hier["repeat_name"].astype(str)
            close = names[names.str.contains(repeat_name, case=False, regex=False)].tolist()[:5]
            hint = f" Did you mean: {', '.join(close)}?" if close else ""
            raise ReferenceError(
                f"repeat '{repeat_name}' not found in {assembly} annotation.{hint}"
            )
        return sub


def _version_string() -> str:
    from tescope import __version__

    return __version__
