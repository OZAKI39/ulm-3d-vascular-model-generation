"""Compare visual reading with a separate primary-PDF text extraction."""
import csv
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
a = [s.split() for s in (ROOT/'raw/CF2003_TABLE17_TRANSCRIPTION_A_VISUAL.txt').read_text().splitlines()
     if s and not s.startswith('#')]
text = (ROOT/'raw/CF2003_p407_PDF_TEXT.txt').read_text().replace('·', '.').replace('−', '-')
lines = text.splitlines()
start = lines.index('0')
b = [lines[start+3*i:start+3*i+3] for i in range(20)]
assert [int(x[0]) for x in a] == [int(x[0]) for x in b] == list(range(20))
for row in b:
    assert all(re.fullmatch(r'-?\d+\.\d+(?:e-\d+)?', v) for v in row[1:])
matches = [x == y for x, y in zip(a, b)]
target = ROOT/'reference/CF2003_TABLE17_VERIFICATION.csv'
shutil.copyfile(target, ROOT/'raw/CF2003_TABLE17_PENDING_BEFORE_PRIMARY_PDF.csv')
with target.open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['index','u_A','u_B','u_final','omega_A','omega_B','omega_final','match','source_page'])
    for x,y,ok in zip(a,b,matches):
        w.writerow([x[0],x[1],y[1],x[1] if ok else '',x[2],y[2],x[2] if ok else '',
                    'YES' if ok else 'NO','CF2003 primary PDF p407, PDF page 27'])
record = {'table17_verified':all(matches),'coefficients_per_column':20,'degree':19,
          'all_literal_strings_match':all(matches),
          'A':'Manual visual reading of primary page image, recorded before primary PDF text extraction',
          'B':'PyMuPDF 1.24.10 primary PDF text extraction; ASCII typographic normalization only',
          'independence':'Two extraction methods from one paper, same reviewing agent; not two physical theories',
          'negative_u_indices':[int(x[0]) for x in a if x[1].startswith('-')],
          'negative_omega_indices':[int(x[0]) for x in a if x[2].startswith('-')],
          'primary_pdf_sha256':'2d572961b10cc7b46ae2a666f8dfc80123ee303b118d964c510ae86a739dc3e5',
          'corrigendum_checked':True,'corrigendum_impact':'NONE','normalization_visually_verified':True}
(ROOT/'validation/CF2003_SOURCE_VERIFICATION.json').write_text(json.dumps(record,indent=2)+'\n')
assert all(matches), 'STOP_TABLE17_AMBIGUITY'
print(json.dumps(record))
