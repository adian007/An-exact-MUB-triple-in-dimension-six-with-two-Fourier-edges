#!/usr/bin/env python3
"""Environment check for Agent 6 optimization work."""
import sys
print(f"Python: {sys.version}")
print(f"Python executable: {sys.executable}")

# Check numpy
try:
    import numpy as np
    print(f"numpy: {np.__version__}")
except ImportError as e:
    print(f"numpy: NOT AVAILABLE ({e})")

# Check scipy
try:
    import scipy
    print(f"scipy: {scipy.__version__}")
except ImportError as e:
    print(f"scipy: NOT AVAILABLE ({e})")

# Check CVXPY
try:
    import cvxpy
    print(f"cvxpy: {cvxpy.__version__}")
except ImportError as e:
    print(f"cvxpy: NOT AVAILABLE ({e})")

# Check other useful packages
for pkg in ["matplotlib", "pandas", "jax", "torch", "sympy", "mpmath"]:
    try:
        mod = __import__(pkg)
        ver = getattr(mod, "__version__", "unknown")
        print(f"{pkg}: {ver}")
    except ImportError:
        print(f"{pkg}: NOT AVAILABLE")

# Check Julia
import subprocess
try:
    result = subprocess.run(
        ["/home/adian/julia-1.12.6/bin/julia", "--version"],
        capture_output=True, text=True, timeout=10
    )
    print(f"Julia: {result.stdout.strip()}")
except Exception as e:
    print(f"Julia: NOT AVAILABLE ({e})")