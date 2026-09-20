import glob, re, html as H
s = open("_build_page_nvidia__nemotron-parse-2.0.html", encoding="utf-8", errors="replace").read()
# find api endpoint strings
eps = list(dict.fromkeys(re.findall(r'https://[a-zA-Z0-9.]*api\.nvidia\.com[^"\\ ]{0,80}', s)))
for e in eps[:20]: print("EP:", e)
print("----- model strings -----")
for m in list(dict.fromkeys(re.findall(r'["'"'](nvidia|meta|microsoft|google|Qwen|adept)[^"'"']{0,70}parse[^"'"']{0,40}["'"']', s, re.I)))[:20]:
    print("MODEL:", m)
# find references to 'v1/'
for v in list(dict.fromkeys(re.findall(r'v1/[a-z/]{3,60}[^"\\ ]{0,40}', s)))[:30]:
    print("V1:", v)
