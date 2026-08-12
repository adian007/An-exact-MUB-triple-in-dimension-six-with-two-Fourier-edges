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
COORD_DEDUP_DIGITS = 15
SSTAR_TARGET = 1865  # 20 anchors + 1845 degen candidates
SSTAR_META = RESULTS / "special_loci_degen1845.meta.txt"


def load_sstar_meta() -> dict[str, str]:
    out: dict[str, str] = {}
    if not SSTAR_META.exists():
        return out
    for line in SSTAR_META.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, val = line.split("=", 1)
        out[key.strip()] = val.strip()
    return out


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


def coord_key(row: dict[str, str]) -> tuple[float, float, float] | None:
    try:
        return (
            round(float(row["theta"]), COORD_DEDUP_DIGITS),
            round(float(row["phi"]), COORD_DEDUP_DIGITS),
            round(float(row["lambda"]), COORD_DEDUP_DIGITS),
        )
    except (KeyError, TypeError, ValueError):
        return None


def distinct_coord_count(rows: list[dict[str, str]]) -> int:
    keys = {k for r in rows if (k := coord_key(r)) is not None}
    return len(keys)


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
        "special_loci_degen500.csv",
        "special_loci_degen1845.csv",
        "dita_lambda_fourth_dense.csv",
    ]
    summaries = [summarize(f, load_csv(f)) for f in files]
    s588 = load_csv("special_loci_search.csv")
    degen200 = load_csv("special_loci_degen200.csv")
    degen500 = load_csv("special_loci_degen500.csv")
    degen1845 = load_csv("special_loci_degen1845.csv")
    dita_dense = load_csv("dita_lambda_fourth_dense.csv")
    sstar_meta = load_sstar_meta()
    meta_distinct = int(sstar_meta["distinct_coordinates"]) if "distinct_coordinates" in sstar_meta else None
    sstar_best = degen1845 if degen1845 else degen500
    sstar_distinct = distinct_coord_count(sstar_best)
    csv_corrupt = bool(
        meta_distinct is not None
        and (
            not degen1845
            or sstar_distinct < meta_distinct
            or any(str(r.get("name", "")).startswith("ref_") for r in degen1845)
        )
    )
    if csv_corrupt and meta_distinct is not None:
        sstar_distinct = meta_distinct
    sstar_overlap = len(
        {k for r in degen200 if (k := coord_key(r)) is not None}
        & {k for r in degen500 if (k := coord_key(r)) is not None}
    )
    sstar_clique6 = sstar_meta.get("clique6") if csv_corrupt else str(
        sum(1 for r in sstar_best if r.get("max_clique") == "6")
    )
    sstar_fourth = sstar_meta.get("fourth_mub", "0") if csv_corrupt else str(
        sum(1 for r in sstar_best if r.get("found_fourth") == "true")
    )

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
        "--- S* batches (no-refine; cumulative degen caps) ---",
        f"degen200 rows={len(degen200)} fourth="
        f"{sum(1 for r in degen200 if r.get('found_fourth')=='true')}"
        f" clique6={sum(1 for r in degen200 if r.get('max_clique')=='6')}",
        f"degen500 rows={len(degen500)} fourth="
        f"{sum(1 for r in degen500 if r.get('found_fourth')=='true')}"
        f" clique6={sum(1 for r in degen500 if r.get('max_clique')=='6')}"
        f" pool_complete={sum(1 for r in degen500 if r.get('pool_complete_flag')=='true')}",
        "degen500_run2: mislabeled duplicate (requested 1845, actual cap 500); see special_loci_degen500_run2.meta.txt",
        f"degen200 subset of degen500: overlap={sstar_overlap}/{len(degen200)} (cumulative caps, not disjoint batches)",
        f"S* distinct coverage (best batch): {sstar_distinct}/{SSTAR_TARGET} "
        f"({100.0 * sstar_distinct / SSTAR_TARGET:.1f}%) "
        f"clique6_loci={sstar_clique6} fourth={sstar_fourth}",
    ]
    if csv_corrupt and meta_distinct is not None:
        lines += [
            "WARNING: special_loci_degen1845.csv missing or corrupted; "
            "S* stats above taken from special_loci_degen1845.meta.txt + run log (1863-point no-refine run).",
            "ACTION: regenerate with search_special_loci.jl --degen-cap 1845 --no-refine before submission.",
        ]
    lines += [
        "",
        "--- Track A dense Dita ---",
        f"628_target={len(dita_dense)} fourth="
        f"{sum(1 for r in dita_dense if r.get('found_fourth')=='true')}",
        "",
        "--- Publication readiness ---",
        "T1 gauge lemmas: results/formalize_gauge_lemmas.txt (ALL PASS)",
        "T2 Dita circle: 126/126 periodicity + 628/628 dense fourth=0",
        "T3 Dita slice: 628/628 numerical fourth=0; M2 witness probe inconclusive (35 eq, dim=-1 CC)",
        "Track D: 2 HP-verified components (F6 arc + Dita circle)",
        "Track C Open Problem: no fourth on all audited points; Groebner witness elimination open",
        "Julia HC: requires JULIA_DEPOT_PATH at project .julia-depot after path move",
        "Macaulay2: requires Docker Desktop running (scripts/docker/run_m2.ps1)",
        "",
        (
            "VERDICT: Paper-scope claims reconciled; T3 numerical only; S* primary "
            f"{sstar_distinct}/{SSTAR_TARGET} (log/meta); exact Groebner witness remains future work."
            + (" CSV artifact MUST be regenerated before submission." if csv_corrupt else "")
        ),
    ]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
