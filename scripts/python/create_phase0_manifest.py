#!/usr/bin/env python3
"""Create a hash-bound, read-only snapshot for the Gate 0 handoff."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results" / "phase0_manifest.json"

ARTIFACTS = [
    "Project.toml",
    "Manifest.toml",
    "PLAN_next_steps_4th_mub.md",
    "INTERFACES.md",
    "docs/scientific_status.md",
    "docs/results/final_honest_status.md",
    "docs/results/phase1b_certified_clique_audit_2026-09-29.md",
    "docs/results/fold_sweep_phase2_2026-09-29.md",
    "b3_a4_i4_results/summary.json",
    "b3_a4_i4_results/pool_pi_over_3.npz",
    "scripts/python/w1_pi3_groebner.py",
    "scripts/python/refine_pi3_vector_mp.py",
    "docs/research/reports/phase2_pi3_reconstruction_probe_2026-09-30.md",
    "results/phase2_pi3_B3_mp100.json",
    "results/pi3_phase_groebner_basis.txt",
    "symbolic_export/pi3_exact_B3_ansatz.md",
    "symbolic_export/w1_pi3_exact_ansatz.m2",
    "scripts/python/audit_fold_critical_match.py",
    "docs/research/reports/phase2_fold_critical_match_audit_2026-09-30.md",
    "results/phase2_fold_critical_match_audit.json",
    "scripts/python/fit_exact_pi3_candidate.py",
    "scripts/python/prove_pi3_fourier_symbolic.py",
    "results/phase2_exact_pi3_fourier_fit.json",
    "scripts/python/export_pi3_exact_abc_w1.py",
    "symbolic_export/w1_pi3_exact_abc.m2",
    "scripts/python/verify_pi3_exact_remainders.py",
    "results/phase2_pi3_exact_remainder_verification.json",
    "docs/research/reports/phase2_pi3_exact_remainder_verification_2026-10-01.md",
    "scripts/python/analyze_pool_graph_pi3.py",
    "results/phase2_pi3_pool_graph_audit.json",
    "docs/research/reports/phase2_pi3_pool_graph_audit_2026-10-01.md",
    "scripts/python/w1_pi3_groebner.py",
    "scripts/python/export_pi3_exact_abc_w1.py",
    "symbolic_export/w1_pi3_exact_abc.m2",
    "results/phase2_w1_setup.json",
    "results/track_c_elimination/m2_w1_pi3_exact_abc.log",
    "results/track_c_elimination/m2_w1_pi3_exact_abc_wsl.log",
    "docs/research/reports/phase2_w1_groebner_run_2026-10-01.md",
    "symbolic_export/w1_pi3_coefficient_field_preflight.sing",
    "symbolic_export/w1_pi3_quadratic_square_classes.sing",
    "symbolic_export/w1_pi3_m2_field_smoke.m2",
    "results/phase2_w1_coefficient_field_preflight.log",
    "results/phase2_w1_quadratic_square_classes.log",
    "docs/research/reports/phase2_w1_coefficient_field_2026-10-01.md",
    "scripts/python/audit_certified_fold_singular_roots.py",
    "results/phase2_fold_certified_singular_roots.json",
    "docs/research/reports/phase2_fold_certified_singular_roots_2026-10-01.md",
    "results/campaigns/i3_singular_locus/root_certification.json",
    "scripts/python/audit_julia_pi3_rebuild.py",
    "results/phase2_pi3_julia_rebuild_audit.json",
    "results/campaigns/i3_singular_locus/third_mub_cliques.json",
    "docs/research/reports/phase2_pi3_julia_rebuild_2026-10-01.md",
    "docs/research/reports/alphaxiv_consultation_2026-10-01.md",
    "scripts/python/check_pi3_phase_field_ideal.py",
    "results/phase2_pi3_phase_field_ideal.json",
    "docs/research/reports/phase2_pi3_phase_field_consistency_2026-10-01.md",
    "docs/research/reports/alphaxiv_consultation_2026-10-01.md",
    "scripts/python/fit_pi3_abc_candidate.py",
    "results/phase2_abc_fourier_fit.json",
    "scripts/python/validate_pi3_exact_b3_field.py",
    "results/phase2_exact_pi3_b3_validation.json",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(*args: str) -> str:
    result = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", *args],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() or result.stderr.strip()


def main() -> None:
    files = {}
    for relative in ARTIFACTS:
        path = ROOT / relative
        files[relative] = {
            "exists": path.is_file(),
            "sha256": sha256(path) if path.is_file() else None,
            "size_bytes": path.stat().st_size if path.is_file() else None,
        }

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "repository": str(ROOT),
        "git_head": git("rev-parse", "HEAD"),
        "git_status_short": git("status", "--short", "--branch").splitlines(),
        "python": sys.version,
        "platform": platform.platform(),
        "artifacts": files,
        "execution_contract": {
            "pool_target": "lambda=pi/3; summary key 1.0471975511965976",
            "w1_script": "scripts/python/w1_pi3_groebner.py",
            "next_external_step": "exact or certified-high-precision B3 reconstruction before Groebner elimination",
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
