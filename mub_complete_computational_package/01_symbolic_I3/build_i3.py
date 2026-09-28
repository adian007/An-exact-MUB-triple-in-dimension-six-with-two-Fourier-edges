"""Compatibility entry point for the project's canonical exact I3 builder."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from dita_i3 import main


if __name__ == "__main__":
    main()
