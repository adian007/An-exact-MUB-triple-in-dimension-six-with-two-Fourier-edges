import os, json, urllib.request
def load_env():
    env = {}
    p = os.path.join(os.getcwd(), '.env')
    if os.path.exists(p):
        for line in open(p, encoding='utf-8', errors='replace'):
            line = line.strip()
            if line and not line.startswith('#') and '=' in line:
                k, _, v = line.partition('=')
                env[k.strip()] = v.strip().strip('"').strip("'")
    for k, v in env.items():
        os.environ.setdefault(k, v)
load_env()
key = os.environ.get('NVIDIA_OCR', '')
url = 'https://integrate.api.nvidia.com/v1/models'
req = urllib.request.Request(url, headers={'Authorization': 'Bearer ' + key})
body = urllib.request.urlopen(req, timeout=60).read().decode('utf-8', 'replace')
data = json.loads(body)
for m in sorted(data.get('data', []), key=lambda x: x.get('id','')):
    print(m.get('id'))
