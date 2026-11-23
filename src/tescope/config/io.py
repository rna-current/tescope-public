"""Load/save analysis configs as YAML (or JSON, by extension)."""

from __future__ import annotations

import json
import os
from pathlib import Path

import yaml

from .model import AnalysisConfig, ConfigError


def load_config(path: str | os.PathLike) -> AnalysisConfig:
    p = Path(path)
    if not p.exists():
        raise ConfigError(f"config file not found: {p}")
    text = p.read_text(encoding="utf-8")
    try:
        if p.suffix.lower() == ".json":
            data = json.loads(text)
        else:
            data = yaml.safe_load(text)
    except (yaml.YAMLError, json.JSONDecodeError) as e:
        raise ConfigError(f"{p}: cannot parse config: {e}") from e
    if not isinstance(data, dict):
        raise ConfigError(f"{p}: config must be a mapping")
    cfg = AnalysisConfig.from_dict(data)
    _resolve_track_paths(cfg, p.parent)
    return cfg


def save_config(cfg: AnalysisConfig, path: str | os.PathLike) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    data = cfg.to_dict()
    if p.suffix.lower() == ".json":
        p.write_text(json.dumps(data, indent=2), encoding="utf-8")
    else:
        p.write_text(
            yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8"
        )


def _resolve_track_paths(cfg: AnalysisConfig, base: Path) -> None:
    """Track paths in a config file are relative to the config's directory."""
    for t in cfg.tracks:
        tp = Path(t.file)
        if not tp.is_absolute():
            t.file = str((base / tp).resolve())
    uf = cfg.annotations.user_features_file
    if uf and not Path(uf).is_absolute():
        cfg.annotations.user_features_file = str((base / uf).resolve())
