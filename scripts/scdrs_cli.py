"""
Run the `scdrs` command-line tool with a NumPy 2 compatibility shim.
scdrs 1.0.2 uses `np.float_` (removed in NumPy 2.0) in gearys_c(), which
breaks the heterogeneity test of `perform-downstream`.
Usage: python scdrs_cli.py <scdrs subcommand> [args...]
"""
import runpy
import shutil
import sys

import numpy as np

if not hasattr(np, "float_"):
    np.float_ = np.float64

exe = shutil.which("scdrs")
if exe is None:
    sys.exit("scdrs executable not found on PATH")
sys.argv = [exe] + sys.argv[1:]
runpy.run_path(exe, run_name="__main__")
