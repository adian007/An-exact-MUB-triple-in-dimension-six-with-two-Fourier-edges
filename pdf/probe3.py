import base64, io, os, requests
try:
    from dotenv import load_dotenv
    load_dotenv()
except Exception: pass
key = os.environ.get("NVIDIA_OCR") or ""
b64 = base64.b64encode(b"") # placeholder replaced below
# tiny 1x1 png
png = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg==")
data_url = "data:image/png;base64," + base64.b64encode(png).decode()
for host in ["https://ai.api.nvidia.com","https://integrate.api.nvidia.com"]:
    for m in ["nvidia/llava-3.2-8b-vit","meta/llama-3.2-11b-vision-instruct"]:
        try:
            body = {"model": m, "max_tokens": 64, "messages":[{"role":"user","content":[
                {"type":"text","text":"What color is this image?"},
                {"type":"image_url","image_url":{"url":data_url}}]}]}
            r = requests.post(f"{host}/v1/chat/completions", headers={"Authorization":"Bearer "+key,"Content-Type":"application/json"}, json=body, timeout=60)
            print(f"{host.split('/')[-1]:22s} {m:42s} -> {r.status_code} {r.text[:220]}")
        except Exception as e:
            print(m, host, "EXC", e)
