#!/usr/bin/env python3
"""Public top-level entry point for the genome annotation pipeline."""

from __future__ import annotations

import sys
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "src"))

from eps_workflows.annotation_workflow import main


if __name__ == "__main__":
    raise SystemExit(main())
