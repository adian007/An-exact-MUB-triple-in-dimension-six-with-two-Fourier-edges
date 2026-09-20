import sys, glob, os
import fitz
files = sorted(glob.glob('*.pdf'))
print('total_pdf_files =', len(files))
rows=[]
for p in files:
    try:
        d = fitz.open(p)
        npg = d.page_count
        txt_pages = 0
        total = 0
        for i in range(npg):
            t = d[i].get_text('text')
            total += len(t.strip())
            if len(t.strip()) > 40:
                txt_pages += 1
        rows.append((total, npg, txt_pages, os.path.getsize(p)/1e3, p))
    except Exception as e:
        rows.append((-1, -1, -1, os.path.getsize(p)/1e3, p + '  ERR ' + str(e)))
print(f'{"textchars":>10} {"pages":>5} {"pgs>40ch":>8} {"sizeKB":>8}  file')
for total,npg,tp,kb,p in sorted(rows):
    print(f'{total:10d} {npg:5d} {tp:8d} {kb:8.1f}  {p}')
