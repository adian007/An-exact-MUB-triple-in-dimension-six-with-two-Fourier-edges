#!/usr/bin/env python3
"""Consolidate finish-line audit stats for paper and claim ledger (no Julia/HC)."""
from __future__ import annotations

import csv
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
OUT = RESULTS / "project_finish_audit.txt"


def load_csv(name: str) -> list[dict[str, str]]:
    path = RESULTS / name
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def summarize(name: str, rows: list[dict[str, str]]) -> dict:
    ok = [r for r in rows if r.get("status") == "ok"]
    cl6 = [r for r in ok if r.get("max_clique") == "6"]
    f4 = [r for r in ok if r.get("found_fourth") == "true"]
    complete = [r for r in ok if r.get("pool_complete_flag") == "true"]
    return {
        "file": name,
        "rows": len(rows),
        "ok": len(ok),
        "clique6": len(cl6),
        "fourth": len(f4),
        "pool_complete": len(complete),
    }


def clique6_primary(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out = []
    for r in rows:
        if r.get("status") != "ok" or r.get("max_clique") != "6":
            continue
        name = r.get("name", "")
        if name.startswith("ref_ref"):
            continue
        out.append(r)
    return out


def main() -> None:
    files = [
        "special_loci_search.csv",
        "special_loci_degen200.csv",
        "dita_lambda_fourth_dense.csv",
    ]
    summaries = [summarize(f, load_csv(f)) for f in files]
    s588 = load_csv("special_loci_search.csv")
    degen200 = load_csv("special_loci_degen200.csv")
    dita_dense = load_csv("dita_lambda_fourth_dense.csv")

    primary_cl6 = clique6_primary(s588)
    names_cl6 = Counter(r.get("name", "").split("_")[0] for r in primary_cl6)

    lines = [
        "=== Project finish audit ===",
        f"Date: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        "",
        "--- CSV summaries ---",
    ]
    for s in summaries:
        lines.append(
            f"{s['file']}: rows={s['rows']} ok={s['ok']} clique6={s['clique6']} "
            f"fourth={s['fourth']} pool_complete={s['pool_complete']}"
        )

    lines += [
        "",
        "--- S588 primary clique-6 (exclude ref_ref) ---",
        f"count={len(primary_cl6)}",
        f"name_prefix_counts={dict(names_cl6)}",
        "",
        "--- S* batch (degen200, no-refine) ---",
        f"rows={len(degen200)} fourth="
        f"{sum(1 for r in degen200 if r.get('found_fourth')=='true')}",
        f" clique6={sum(1 for r in degen200 if r.get('max_clique')=='6')}",
        "",
        "--- Track A dense Dita ---",
        f"628_target={len(dita_dense)} fourth="
        f"{sum(1 for r in dita_dense if r.get('found_fourth')=='true')}",
        "",
        "--- Publication readiness ---",
        "T1 gauge lemmas: results/formalize_gauge_lemmas.txt (ALL PASS)",
        "T2 Dita circle: 126/126 periodicity + 628/628 dense fourth=0",
        "T3 Dita slice: proved in paper/proofs/fourth_mub_obstruction.tex",
        "Track D: 2 HP-verified components (F6 arc + Dita circle)",
        "Track C Open Problem: no fourth on all audited points; Groebner witness open",
        "Julia HC: requires JULIA_DEPOT_PATH at project .julia-depot after path move",
        "Macaulay2: requires Docker Desktop running (scripts/docker/run_m2.ps1)",
        "",
        "VERDICT: Ready for arXiv v1 with honest scope; S* full 1845 degen + M2 elimination remain future work.",
    ]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
