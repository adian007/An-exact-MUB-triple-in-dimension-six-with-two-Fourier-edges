import urllib.request, json, os
p = '.env'
d = {}
for l in open(p, encoding='utf-8'):
    k, _, v = l.partition('=')
    d[k.strip()] = v.strip().strip('"').strip("'").strip()
key = d.get('NVIDIA_OCR', '')
req = urllib.request.Request('https://integrate.api.nvidia.com/v1/models', headers={'Authorization': 'Bearer ' + key})
body = urllib.request.urlopen(req, timeout=60).read().decode('utf-8', 'replace')
ids = [m.get('id', '') for m in json.loads(body).get('data', [])]
print('TOTAL', len(ids))
for i in sorted(i for i in ids if any(x in i.lower() for x in ('ocr', 'tesseract', 'paddle'))):
    print('OCR  ', i)
for i in sorted(i for i in ids if 'sky-lizard' in i.lower()):
    print('SKY  ', i)
