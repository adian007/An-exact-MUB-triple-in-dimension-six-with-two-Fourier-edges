#!/usr/bin/env python3
"""Verify all \\input{} targets exist before LaTeX compile (Overleaf/local)."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper"
MAIN = PAPER / "mub6_karlsson_ieee.tex"


def find_inputs(tex: Path) -> list[str]:
    text = tex.read_text(encoding="utf-8")
    return re.findall(r"\\input\{([^}]+)\}", text)


def resolve_input(base: Path, target: str) -> Path:
    p = Path(target)
    if p.suffix != ".tex":
        p = p.with_suffix(".tex")
    return (base.parent / p).resolve()


def check_math_mode_escaped_underscore(tex: Path) -> list[str]:
    """Flag \\_{...} inside $...$ (invalid in math mode)."""
    issues: list[str] = []
    text = tex.read_text(encoding="utf-8")
    for m in re.finditer(r"\$[^$]+\$", text):
        chunk = m.group(0)
        if r"\_" in chunk and r"\texttt" not in chunk:
            line = text[: m.start()].count("\n") + 1
            issues.append(f"{tex.relative_to(ROOT)}:{line}: math-mode \\\\_: {chunk[:60]}...")
    return issues


def main() -> int:
    missing: list[str] = []
    underscore_issues: list[str] = []

    for main in (PAPER / "main.tex", MAIN):
        if not main.exists():
            continue
        for inp in find_inputs(main):
            path = resolve_input(main, inp)
            if not path.exists():
                missing.append(str(path.relative_to(ROOT)))

    for tex in [MAIN, *sorted((PAPER / "proofs").glob("*.tex"))]:
        underscore_issues.extend(check_math_mode_escaped_underscore(tex))

    if missing:
        print("MISSING LaTeX inputs (upload these to Overleaf):")
        for m in missing:
            print(f"  - {m}")
        print("\nUpload the entire paper/ folder, including paper/proofs/*.tex")
        print("Or use paper/main_standalone.tex (single file).")
        return 1

    if underscore_issues:
        print("WARN: possible math-mode underscore issues:")
        for u in underscore_issues:
            print(f"  - {u}")

    print("OK: all \\input{} targets found for paper/main.tex and mub6_karlsson_ieee.tex")
    return 0 if not underscore_issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
