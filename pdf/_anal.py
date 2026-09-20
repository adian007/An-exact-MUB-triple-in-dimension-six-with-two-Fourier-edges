import io, re, glob, os
def rd(p):
    with io.open(p,'r',encoding='utf-8',errors='replace') as f: return f.read()
def report(p):
    c = rd(p)
    out = []
    def hits(pat, each=2, ctx=60, cap=6):
        mvs = list(re.finditer(pat, c))
        out.append('  PAT %-46s -> %d' % (pat, len(mvs)))
        seen = set(); n = 0
        for m in mvs:
            s = max(0, m.start()-ctx); e = min(len(c), m.end()+ctx)
            v = c[s:e].replace('\r','\\r').replace('\n','\\n')
            if v in seen: continue
            seen.add(v); n += 1
            if len(v) > 150: v = v[:150] + '...'
            out.append('     ' + v)
            if n >= each: break
    print('='*40)
    print('FILE:', os.path.basename(p))
    for pat in [
        r'render_html_content', r'render_is_visible_content',
        r'"rows"\s*:', r'output_markdown', r'predict_bbox',
        r'predict_boundary_boxes', r'predict_text_in_pic', r'predict_page_annotation',
        r'predict_xf_heading', r'predict_is_visible_content', r'predict_read_order',
        r'end_org1', r'render_show_media', r'render_limit_concurrent_reads',
        r'"input"\s*:\s*{', r'"image"\s*:\s*"', r'application/pdf',
        r'base64', r'use_ocr', r'img_ann', r'"prompt"\s*:\s*"', r'xsd_string',
        r'user_boolean', r'output_mean_char_and_line', r'with_image', r'image_url',
        r'"receipt"', r'quotation',
    ]:
        hits(pat, each=2, ctx=50)
    return out
if __name__ == '__main__':
    files = sorted(glob.glob('D:\\MUBs in 6-dimension\\pdf\\*.html'), key=lambda s: os.path.getmtime(s))
    for f in files:
        print('FILE-LIST', os.path.basename(f))
    for f in files:
        report(f)
