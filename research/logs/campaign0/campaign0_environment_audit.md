# Campaign 0 — environment and repository audit

**Timestamp:** 2026-09-17T22:14:03+02:00  
**Scope:** environment, inventory, provenance, and available executables only. No mathematical solver, sweep, elimination, or numerical-discovery calculation was run.

## Commands executed

```powershell
$PSVersionTable.PSVersion
Get-Command julia,python,python3,git,docker,macauley2,M2,sage,mathematica,wolframscript,phc
python --version
julia --version
docker --version
Get-FileHash Project.toml,Manifest.toml,src\MubSearch.jl,test\runtests.jl,paper\main_theorems.tex,paper\methods_audit.tex -Algorithm SHA256
wsl.exe -l -q
git status --short --branch
```

## Environment observed

| Item | Result | Status / limitation |
|---|---|---|
| Host shell | PowerShell 5.1.26100.9444 | available |
| Python | 3.13.14 at `C:\Users\adian\AppData\Local\Programs\Python\Python313\python.exe` | available |
| Git executable | 2.45.1.1 | available, but repository ownership protection prevents commands |
| Docker executable | 29.7.2 | executable available; its user config is inaccessible in this sandbox |
| WolframScript | 1.4.0.0 | executable available |
| Julia on host `PATH` | not found | BLOCKED; repository documentation attributes this to Windows Smart App Control |
| WSL enumeration | access denied (`Wsl/EnumerateDistros/Service/E_ACCESSDENIED`) | BLOCKED in this execution environment |
| Macaulay2, Sage, PHCpack, Bertini | not found on `PATH` | NOT_RUN / unavailable from host PATH |
| OS and CPU via CIM | access denied | NOT_RUN; hardware capacity not inferred |

The pinned project is `MubSearch` version `0.2.0`, with compat constraints including HomotopyContinuation 2.22, Graphs 1.14, Nemo 0.56, CSV 0.10, and DataFrames 1.8. The repository contains 180 source/document files with extensions `.jl`, `.py`, `.m2`, `.wl`, `.tex`, or `.md`.

## Provenance hashes at audit time

| Artifact | SHA-256 |
|---|---|
| `Project.toml` | `D5DB9E4CFAD3C1B75FD84BB7555C038889E4B967870347049CDD55822EF04B17` |
| `Manifest.toml` | `3DD635622DAEA226E25CA75EF6E12BC433049B13511189055EE797347B3FE229` |
| `src/MubSearch.jl` | `891BD6E4E33F7FF72A525013FB20911E991025B232E559C59D5F9E13309B69E6` |
| `test/runtests.jl` | `8EF9F339FBC84183B2282990BB77F6FEF1090BE8AA31CCABEDB0D7C517B744D1` |
| `paper/main_theorems.tex` | `E77E1766005274B27AFBDE00E43EE49FD48638B23A0ADDC02F0D9A7816A9F6E3` |
| `paper/methods_audit.tex` | `110219E7C447A28831F4C3CC124A36A9356306A149F8CFAC5BFEB7C55A090F84` |

`git status` did not run because Git reports dubious ownership: the workspace owner SID differs from the sandbox user SID. No global `safe.directory` setting was changed. A pre-existing benchmark artifact records commit `7bf60d4aaf2881030b5d376c1757706f30107ff1`, Julia 1.12.6, HomotopyContinuation 2.22.1, Linux, and seed 20260914; this is recorded repository provenance, not a Campaign 0 computation.

## Research structure initialized

Agent-specific directories were created under `research/` for parameter spaces, methods, a possible four-parameter model, exact systems, interval certification, optimization, symmetry, independent reproduction, literature, and referee review. The campaign ledger will use only the statuses specified in the research directive.

## Strongest conclusion

Campaign 0 is complete as an environment audit. Mathematical computations are **NOT_RUN**. Current execution cannot run Julia, WSL, or a full Git provenance check from the host sandbox; any later run must record the actual execution environment and raw command output.

## Next action

Launch the requested literature, parameter-space, numerical-algebraic-method, exact-formulation, and symmetry/gauge agents. They must not launch broad sweeps or large Gröbner calculations. Their reports require adversarial review before Campaign 1.
