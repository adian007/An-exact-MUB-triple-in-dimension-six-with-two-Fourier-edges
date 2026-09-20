import io, sys, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
t = open('model_report.txt', encoding='utf-8', errors='replace').read()
toks = [x.strip().strip("\"'").strip() for x in re.split(r'\s*,\s*|\s+OR\s+|\s+or\s+', t)]
toks = [x for x in toks if x]
print("COUNT", len(toks))
for i, x in enumerate(toks):
    print(i, repr(x))
print("---GREP---")
for x in toks:
    u = x.lower()
    for kw in ['nvidia','nemo','docx','ocrm','mii','nvd','llama-3.2-3b','3b-instruct','-3b']:
        if kw in u:
            print("MATCH", kw, "->", x); break
