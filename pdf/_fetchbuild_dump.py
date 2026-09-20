import requests, re, json, html
H = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124 Safari/537.36"}
for mid in ["nvidia/nemotron-parse-2.0","nvidia/htc-ocr-2","nvidia/llava-3.2-8b-vit","nvidia/mumu-0.5-8b","nvidia/ai-classification-ocr-2","nvidia/phi-3.5-vision-instruct"]:
    try:
        r = requests.get(f"https://build.nvidia.com/{mid}", headers=H, timeout=30)
        print("="*16, mid, "http", r.status_code, "len", len(r.text))
        m = re.search(r'"(base_url|path|inferenceUrl|endpoint|api")', r.text)
        m2 = re.search(r'nemotron-parse|htc-ocr|classification-ocr|api\.nvidia\.com', r.text)
        print("  stub:", (m and m.group(0)), "| hit:", (m2 and m2.group(0)))
        # print first occurrence of api.nvidia.com endpoint string
        for mm in re.finditer(r'https://[a-z0-9.\-]*nvidia\.com[^"\\ ]*', r.text):
            s = mm.group(0)
            if 'api' in s or 'org' in s:
                print("   URL:", s)
        open(f"_build_page_{mid.replace('/','__')}.html","w",encoding="utf-8",errors="replace").write(r.text)
    except Exception as e:
        print(mid, "ERR", e)

