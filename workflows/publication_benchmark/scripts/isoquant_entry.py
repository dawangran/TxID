#!/usr/bin/env python3
"""Run the installed IsoQuant while placing its mutable config in a declared path.

IsoQuant 3.13.0 hard-codes ``$HOME/.config/IsoQuant``. Publication workflows
often run with a read-only home directory, so this entry point changes only the
four cache/config file paths and delegates all analysis to IsoQuant itself.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import isoquant


def _set_configs_directory(args) -> None:
    raw = os.environ.get("TXID_ISOQUANT_CONFIG_DIR")
    if not raw:
        raise RuntimeError("TXID_ISOQUANT_CONFIG_DIR is required")
    config_dir = Path(raw)
    config_dir.mkdir(parents=True, exist_ok=True)
    args.db_config_path = str(config_dir / "db_config.json")
    args.index_config_path = str(config_dir / "index_config.json")
    args.bed_config_path = str(config_dir / "bed_config.json")
    args.alignment_config_path = str(config_dir / "alignment_config.json")
    for raw_path in (
        args.db_config_path,
        args.index_config_path,
        args.bed_config_path,
        args.alignment_config_path,
    ):
        path = Path(raw_path)
        if not path.exists():
            path.write_text(json.dumps({}, indent=2) + "\n", encoding="utf-8")


isoquant.set_configs_directory = _set_configs_directory
isoquant.main_entry()
