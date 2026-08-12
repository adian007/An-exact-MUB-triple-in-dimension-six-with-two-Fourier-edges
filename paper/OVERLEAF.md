# Overleaf upload

## Recommended (single file)

1. Regenerate after editing the source paper:
   ```powershell
   python scripts/python/bundle_paper_tex.py
   ```
2. Upload **`paper/main.tex`** to Overleaf (Overleaf compiles `main.tex` by default).
3. Recompile **twice** (first pass writes `.aux` labels/citations; second pass resolves refs).

No `proofs/` folder required — proof appendices are inlined automatically.

`main_standalone.tex` is a duplicate of `main.tex` for convenience.

## Multi-file (local development)

Compile from the `paper/` directory:

```powershell
cd paper
pdflatex mub6_karlsson_ieee.tex
pdflatex mub6_karlsson_ieee.tex
```

Or use `main_multifile.tex` (wrapper around `mub6_karlsson_ieee.tex`).

Requires `paper/proofs/*.tex` on disk. Verify:

```powershell
python scripts/python/check_paper_latex.py
```

## Common errors

| Error | Cause | Fix |
|-------|-------|-----|
| `File 'proofs/gauge_and_locus.tex' not found` | Old `main.tex` wrapper without `proofs/` uploaded | Run `bundle_paper_tex.py` and re-upload **`main.tex`** |
| Many `undefined reference` / `undefined citation` on pass 1 | Normal first-pass warnings | Recompile a **second time** |
| `output.pdf` exists | Rare Overleaf naming conflict | Rename/delete `output.pdf` in the project |

## Source files (edit these, not main.tex)

- `mub6_karlsson_ieee.tex` — main body
- `paper/proofs/gauge_and_locus.tex`
- `paper/proofs/dita_third_mub.tex`
- `paper/proofs/fourth_mub_obstruction.tex`
