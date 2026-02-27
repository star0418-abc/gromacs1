#!/usr/bin/env python3
"""
MD Simulation Pipeline Launcher

PACKMOL → HTPOLYNET → GROMACS → RDF/CN

Usage:
    python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --stage all
    
See --help for all options.
"""

import sys
from pathlib import Path

# Ensure pipeline package is importable
sys.path.insert(0, str(Path(__file__).parent))

from pipeline.cli import main

if __name__ == "__main__":
    sys.exit(main())
