import re, glob, json, html as H
for f in glob.glob("_build_page_*.html"):
    s = open(f, encoding="utf-8", errors="replace").read()
    print("="*10, f, len(s))
    for pat in [r'nemotron[^"\\]{0,60}', r'htc-ocr[^"\\]{0,40}', r'classification-ocr[^"\\]{0,40}']:
        hits = list(dict.fromkeys(re.findall(pat, s, re.I)))
        if hits: print("   ", pat, "->", hits[:6])
