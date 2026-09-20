import os, sys, concurrent.futures
import requests
try:
    from dotenv import load_dotenv
    load_dotenv()
except Exception:
    pass
api_key = os.environ.get("NVIDIA_OCR") or os.environ.get("NVIDIA_API_KEY") or ""
H = {"Authorization": "Bearer " + api_key, "Accept": "application/json"}

def get(id):
    try:
        r = requests.get(f"https://ai.api.nvidia.com/v1/models/{id}", headers=H, timeout=30)
        if r.status_code != 200:
            return f"{id:50s} http={r.status_code}", None
        j = r.json()
        return f"{id:50s} http=200 class={j.get('class')} created={j.get('created')} media={j.get('media_types')}", None
    except Exception as e:
        return f"{id:50s} EXC {e}", None

ids = [x.strip() for x in sys.argv[1:]]
out = []
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
    for line, _ in ex.map(get, ids):
        if line:
            out.append(line)
            print(line, flush=True)
with open("_probe_results.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(out))
