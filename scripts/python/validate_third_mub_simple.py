#!/usr/bin/env python3
"""Validate third-MUB candidates from special_loci_search.csv (REFactored 2026-09-29)"""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
CSV_PATH = RESULTS / "special_loci_search.csv"
DITA = math.acos(1 / math.sqrt(3))

def main():
    with open(CSV_PATH, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    
    ok_rows = [r for r in rows if r.get("status") == "ok"]
    mc6 = [r for r in ok_rows if int(r.get("max_clique") or 0) >= 6]
    no_ref = [r for r in mc6 if "ref_ref" not in r.get("name", "")]
    
    print(f"Total rows: {len(rows)}")
    print(f"OK rows: {len(ok_rows)}")
    print(f"max_clique>=6: {len(mc6)}")
    print(f"No ref_ref in name (Python exclude): {len(no_ref)}")
    print(f"Has ref_ref in name: {len(mc6) - len(no_ref)}")
    
    # Write summary
    (RESULTS / "validation_summary_py.txt").write_text(
        f"n_csv_rows={len(rows)}\n"
        f"n_ok_rows={len(ok_rows)}\n"
        f"n_candidates_exclude_ref_ref={len(no_ref)}\n"
        f"n_candidates_include_ref_ref={len(mc6)}\n",
        encoding="utf-8"
    )
    print("Wrote validation_summary_py.txt")

if __name__ == "__main__":
    main()