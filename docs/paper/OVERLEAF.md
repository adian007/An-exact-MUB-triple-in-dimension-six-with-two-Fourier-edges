# Overleaf upload

## Paper 1 (theorems)

Upload **`docs/paper/main_theorems.tex`** as Overleaf `main.tex` (self-contained).

Edit locally: `main_theorems_multifile.tex` + `preamble_common.tex` + `proofs/gauge_structure.tex` + `proofs/fourth_mub_theorems.tex`, then:

```powershell
python scripts/python/bundle_main_theorems.py
```

## Paper 2 (methods/audit)

Upload **`paper/methods_audit.tex`** as Overleaf `main.tex` (self-contained).

Edit locally: `methods_audit_multifile.tex` + `preamble_common.tex` + `proofs/locus_geometry.tex` + `proofs/dita_third_mub_methods.tex`, then:

```powershell
python scripts/python/bundle_methods_audit.py
```

## Local multifile compile

```powershell
cd paper
pdflatex main_theorems_multifile.tex
pdflatex methods_audit_multifile.tex
```
