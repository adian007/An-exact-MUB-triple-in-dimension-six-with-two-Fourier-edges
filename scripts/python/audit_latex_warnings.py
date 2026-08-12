#!/usr/bin/env python3
"""Parse mub6_karlsson_ieee.log and emit NDJSON diagnostics to debug-388a52.log."""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOG = ROOT / "paper" / "mub6_karlsson_ieee.log"
DEBUG = ROOT / "debug-388a52.log"
SESSION = "388a52"


def main() -> int:
    if not LOG.exists():
        payload = {
            "sessionId": SESSION,
            "runId": "audit",
            "hypothesisId": "baseline",
            "location": "audit_latex_warnings.py",
            "message": "log file missing",
            "data": {"path": str(LOG)},
            "timestamp": int(time.time() * 1000),
        }
        DEBUG.parent.mkdir(parents=True, exist_ok=True)
        with DEBUG.open("a", encoding="utf-8") as f:
            f.write(json.dumps(payload) + "\n")
        print(f"MISSING: {LOG}")
        return 1

    text = LOG.read_text(encoding="utf-8", errors="replace")
    patterns = {
        "hyperref_pdf_string": r"Package hyperref Warning: Token not allowed",
        "float_h_to_ht": r"float specifier changed to `ht'",
        "underfull_hbox": r"Underfull \\hbox",
    }
    counts = {k: len(re.findall(v, text)) for k, v in patterns.items()}

    # #region agent log
    with DEBUG.open("a", encoding="utf-8") as f:
        f.write(
            json.dumps(
                {
                    "sessionId": SESSION,
                    "runId": "audit",
                    "hypothesisId": "summary",
                    "location": "audit_latex_warnings.py:main",
                    "message": "warning counts",
                    "data": counts,
                    "timestamp": int(time.time() * 1000),
                }
            )
            + "\n"
        )
        for m in re.finditer(
            r"^(Package hyperref Warning:.*|LaTeX Warning: `h' float specifier.*|Underfull \\hbox.*)$",
            text,
            re.MULTILINE,
        ):
            f.write(
                json.dumps(
                    {
                        "sessionId": SESSION,
                        "runId": "audit",
                        "hypothesisId": "detail",
                        "location": "mub6_karlsson_ieee.log",
                        "message": m.group(1)[:200],
                        "data": {},
                        "timestamp": int(time.time() * 1000),
                    }
                )
                + "\n"
            )
    # #endregion

    print(json.dumps(counts, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
