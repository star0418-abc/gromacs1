#!/usr/bin/env python3
"""
Direct CLI launcher for the MD simulation pipeline.

Usage:
    python cli.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage all
"""

import sys
from pathlib import Path

# Ensure pipeline package is importable when launching from repo root.
sys.path.insert(0, str(Path(__file__).parent))

from pipeline.cli import main


if __name__ == "__main__":
    sys.exit(main())
