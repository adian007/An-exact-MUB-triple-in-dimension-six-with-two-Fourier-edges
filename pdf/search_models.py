import io, sys, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
t = open('model_report.txt', encoding='utf-8', errors='replace').read()
toks = [x.strip().strip('"' ).strip(chr(39)).strip() for x in re.split(r'\s*,\s*|\s+OR\s+|\s+or\s+', t)]
toks = [x for x in toks if x]
print('TOTAL_TOKENS', len(toks))
seen=set()
for x in toks:
    u = x.lower()
    if any(k in u for k in ['ocr','docx','parse','paddle','tesseract','vision','kosmos','deplot','kosmos-2']):
        if x not in seen:
            seen.add(x); print('VISION/OCR:', x)
print('---3B---')
for x in toks:
    u=x.lower()
    if u.endswith('3b-instruct') or '3b-instruct' in u or '-3b' in u or '3.2-3b' in u or '3.1-3b' in u:
        if x not in seen:
            seen.add(x); print('3B:', x)
