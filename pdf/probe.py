import concurrent.futures, json, os, sys, time
import requests

try:
    from dotenv import load_dotenv
    load_dotenv()
except Exception as e:
    print('dotenv not available:', e)

api_key = os.environ.get("NVIDIA_OCR") or os.environ.get("NVIDIA_API_KEY") or os.environ.get("NVIDIA_NIM_API_KEY") or ""
print("KEY_SET =", bool(api_key), "len =", len(api_key))

def get(id):
    r = requests.get(f"https://ai.api.nvidia.com/v1/models/{id}", headers={"Authorization": "Bearer "+api_key}, timeout=30)
    if r.status_code != 200:
        return id, r.status_code, None
    j = r.json()
    md = j.get("media_types", [])
    c = j.get("created", "")
    mfn = j.get("class")
    return id, r.status_code, {"class": mfn, "created": c, "media_types": md}

def main(ids):
    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        for id, c, info in ex.map(get, ids):
            rows.append((id, c, info))
            print(f"{id:55s} http={c} class={info and info['class']} created={info and info['created']} media={info and info['media_types']}", flush=True)
    return rows

if __name__ == "__main__":
    ids = sys.argv[1:]
    if not ids:
        print("usage: python probe.py <model_id> [more ids...]")
    main(ids)
