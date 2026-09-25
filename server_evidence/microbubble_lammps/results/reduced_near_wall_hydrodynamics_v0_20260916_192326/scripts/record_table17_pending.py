"""Record unverified text evidence; deliberately does not evaluate a closure."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
raw = ROOT / 'raw/CF2003_TABLE17_TRANSCRIPTION_B_UNVERIFIED.txt'
rows = [line.split() for line in raw.read_text().splitlines()
        if line and not line.startswith('#')]
assert [int(row[0]) for row in rows] == list(range(20))
target = ROOT / 'reference/CF2003_TABLE17_VERIFICATION.csv'
assert not target.exists(), 'Evidence already exists; do not overwrite'
with target.open('w', newline='') as stream:
    out = csv.writer(stream)
    out.writerow(['index', 'u_A', 'u_B', 'u_final', 'omega_A', 'omega_B',
                  'omega_final', 'match', 'source_page'])
    for index, u, omega in rows:
        out.writerow([index, '', u, '', '', omega, '', 'NOT_VERIFIED',
                      'CF2003 p407; Doczz web text only; primary image unavailable'])
print('20 pending rows written; 0 visually verified; no final coefficients')
