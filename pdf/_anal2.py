import io, glob, os, re
def rd(p):
    with io.open(p,'r',encoding='utf-8',errors='replace') as f: return f.read()
def counts(p):
    c = rd(p); L = []
    for pat in [
        'render_html_content','render_is_visible_content','"rows"','output_markdown',
        'predict_bbox','predict_boundary_boxes','predict_text_in_pic','predict_page_annotation',
        'predict_xf_heading','predict_is_visible_content','predict_read_order',
        'render_show_media','render_limit_concurrent_reads','"input"','application/pdf','base64',
        'use_ocr','img_ann','"prompt"','xsd_string','user_boolean','output_mean_char_and_line',
        'with_image','image_url','"receipt"','quotation','nemotron','Phi-3.5','llama','mumu',
    ]:
        n = len(re.findall(re.escape(pat), c))
        if n: L.append('%s=%d' % (pat, n))
    return ', '.join(L)
for p in sorted(glob.glob('*.html'), key=lambda s: os.path.getmtime(s)):
    print(os.path.basename(p))
    print('   ', counts(p))
    print()
