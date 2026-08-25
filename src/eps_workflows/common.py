"""Shared validation and script-rendering utilities."""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Iterable


SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("The configuration root must be a JSON object")
    return data


def require_keys(data: dict[str, Any], keys: Iterable[str], context: str) -> None:
    missing = [key for key in keys if data.get(key) in (None, "")]
    if missing:
        raise ValueError(f"Missing required {context} keys: {', '.join(missing)}")


def validate_identifier(value: str, field: str) -> str:
    if not SAFE_IDENTIFIER.fullmatch(value):
        raise ValueError(f"{field} contains unsafe characters: {value!r}")
    return value


def write_script(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")
    if os.name != "nt":
        path.chmod(path.stat().st_mode | 0o111)


def execute_script(path: Path) -> None:
    subprocess.run(["bash", str(path)], check=True)
