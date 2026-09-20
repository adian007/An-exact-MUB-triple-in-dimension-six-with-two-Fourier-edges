import io, json
p = 'build_manifest.json'
c = io.open(p,'r',encoding='utf-8',errors='replace').read()
try:
    d = json.loads(c)
    print('TOP KEYS:', list(d.keys())[:50])
    import re as _r
    def walk(o, pre=''):
        if isinstance(o, dict):
            for k,v in o.items():
                if isinstance(v,(dict,list)):
                    print(pre+str(k)); walk(v, pre+'  ')
                else:
                    print(pre+str(k), '=', str(v)[:80])
        elif isinstance(o, list):
            for i,v in enumerate(o):
                print(pre+'[%d]'%i); walk(v, pre+'  ')
    walk(d)
except Exception as e:
    print('JSON ERR', e)
    print(c[:8000])
