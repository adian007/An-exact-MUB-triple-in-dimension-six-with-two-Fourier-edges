import re
s = open("_build_page_nvidia__nemotron-parse-2.0.html", encoding="utf-8", errors="replace").read()
i = s.find("parse")
while i != -1:
    seg = s[max(0,i-200):i+300]
    if any(x in seg.lower() for x in ["3b", "2b", "billion", "params", "model id", "model_id", "OCR/parse", "1.5", "2025"]):
        print("======= CTX near", i, "=======")
        print(re.sub(r'<[^>]+>',' ', seg)[:600])
    i = s.find("parse", i+1)
