import io, json, urllib.request
d = {}
for l in open('.env', encoding='utf-8'):
    k, _, v = l.partition('=')
    d[k.strip()] = v.strip().strip('"').strip("'").strip()
key = d.get('NVIDIA_OCR', '')
req = urllib.request.Request('https://integrate.api.nvidia.com/v1/models',
                             headers={'Authorization': 'Bearer ' + key})
body = urllib.request.urlopen(req, timeout=60).read().decode('utf-8', 'replace')
data = json.loads(body)
ids = [m.get('id', '') for m in data.get('data', []) if m.get('id')]
need = dict(ocr=0, tesseract=0, paddle=0, sky=0)
lots = []
for i in sorted(ids):
    low = i.lower()
    if 'ocr' in low or 'tesseract' in low or 'paddle' in low:
        lots.append(('OCR_EXACT', i))
    elif 'sky-lizard' in low:
        lots.append(('SKY_LIZARD', i))
    else:
        continue
out = io.StringIO()
out.write('TOTAL %d\n' % len(ids))
if lots:
    for tag, i in lots:
        out.write('%s  %s\n' % (tag, i))
    # show a bit more context for OCR-ish ids only (may be alias/date suffixes)
    for tag, i in lots:
        out.write('   meta: %s\n' % filtered_meta(data, i))
else:
    out.write('NO EXACT OCR MATCHES. Sample ids: %s\n' % ', '.join(ids[:60]))
def filtered_meta(d, mid):
    for m in d.get('data', []):
        if m.get('id') == mid:
            return {k: m[k] for k in ('object','created','owned_by') if k in m}
    return {}
open('model_report.txt', 'w', encoding='utf-8').write(out.getvalue())
