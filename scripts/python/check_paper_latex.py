#!/usr/bin/env python3
"""Verify all \\input{} targets exist for the active multifile paper sources."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper"
SOURCES = [
    PAPER / "main_theorems_multifile.tex",
    PAPER / "methods_audit_multifile.tex",
]


def find_inputs(tex: Path) -> list[str]:
    text = tex.read_text(encoding="utf-8")
    return re.findall(r"(?<!%)\\input\{([^}]+)\}", text)


def resolve_input(base: Path, target: str) -> Path:
    p = Path(target)
    if p.suffix != ".tex":
        p = p.with_suffix(".tex")
    return (base.parent / p).resolve()


def main() -> int:
    missing: list[str] = []
    for main in SOURCES:
        if not main.exists():
            missing.append(f"(source missing) {main.relative_to(ROOT)}")
            continue
        for inp in find_inputs(main):
            # skip commented lines roughly handled by lookbehind; also skip if line is comment
            path = resolve_input(main, inp)
            if not path.exists():
                missing.append(str(path.relative_to(ROOT)))

    if missing:
        print("MISSING LaTeX inputs:")
        for m in missing:
            print(f"  - {m}")
        return 1

    print(
        "OK: all \\input{} targets found for "
        "main_theorems_multifile.tex and methods_audit_multifile.tex"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
