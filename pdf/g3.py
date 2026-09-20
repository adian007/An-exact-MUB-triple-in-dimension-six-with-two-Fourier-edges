import re, html, json
s = open("_build_page_nvidia__nemotron-parse-2.0.html", encoding="utf-8", errors="replace").read()
# find contexts around "3B" / "3b" / "billion"
for pat in [r'3\s*[bB]', r'billion', r'1\.5\s*[bB]', r'2\s*[bB]']:
    for m in re.finditer(pat, s):
        i = m.start()
        seg = re.sub(r'<[^>]+>',' ', s[max(0,i-260):i+260])
        seg = html.unescape(seg).replace("\\n"," ").replace("\\u003c","<").replace("\\u003e",">")
        seg = re.sub(r'\s+',' ', seg)
        print("### CTX ###", repr(seg))
