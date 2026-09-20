import re
s = open("_build_page_nvidia__nemotron-parse-2.0.html", encoding="utf-8", errors="replace").read()
pats = [
  r'https://[a-zA-Z0-9./\-]{0,90}integrate[a-zA-Z0-9./\-]{0,60}',
  r'https://[a-zA-Z0-9./\-]{0,60}api\.nvidia[a-zA-Z0-9./\-]{0,60}',
  r'[A-Za-z0-9_./\-]{0,30}chat/completions[A-Za-z0-9_./\-?=&]{0,40}',
  r'[A-Za-z0-9_./\-]{0,30}ocr/[a-z0-9./\-]{3,50}',
]
seen=set()
out=[]
for p in pats:
  for m in re.findall(p, s):
    if m not in seen:
      seen.add(m); out.append(m)
for m in out[:60]:
  print(m)
print("---- model id near nemotron-parse ----")
for mm in list(dict.fromkeys(re.findall(r'["\x27]([A-Za-z0-9_./-]+parse[A-Za-z0-9_./-]*)["\x27]', s)))[:10]:
  print(repr(mm))
print("---- any 'v1' base url contexts ----")
for m in list(dict.fromkeys(re.findall(r'https://[^"\x27 ]{10,80}(?:v1|integrate)[^"\x27 ]{0,40}', s)))[:30]:
  print(repr(m))
