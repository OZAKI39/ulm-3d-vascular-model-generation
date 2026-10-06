"""Check exported files against archived input; no scientific solve or mutation."""
import base64
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
import plotly
from PIL import Image, __version__ as pillow_version

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'output'
manifest = json.loads((OUT / 'manifests/SOURCE_AND_VALIDATION.json').read_text())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

for item in manifest['source_files']:
    assert sha(ROOT / item['path']) == item['sha256'], item['path']

formats = []
for item in manifest['figures']:
    stem = OUT / 'figures' / item['stem']
    with Image.open(stem.with_suffix('.png')) as im:
        assert list(im.size) == item['dimensions_px']
        assert all(abs(value - 600) < .1 for value in im.info['dpi'])
        im.verify()
    assert stem.with_suffix('.pdf').read_bytes().startswith(b'%PDF-')
    assert ET.parse(stem.with_suffix('.svg')).getroot().tag.endswith('svg')
    formats.append(dict(stem=item['stem'],png_dpi=600,dimensions_px=item['dimensions_px'],formats=['png','pdf','svg']))

row_counts = {}
for name, expected in [('particles_and_internal_force.csv',2197), ('lower_half_signed_layers.csv',7),
                       ('actual_dependency_path.csv',3), ('selected_core_force_channels.csv',4), ('evidence_summary.csv',1)]:
    with (OUT / 'tables' / name).open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == expected
    row_counts[name] = len(rows)

html = (OUT / manifest['interactive_html']).read_text()
assert 'plotly.js v' in html
assert not re.search(r'<script[^>]+src\s*=', html, re.I)
tail = html[html.rindex('Plotly.newPlot(') + len('Plotly.newPlot('):]
decoder = json.JSONDecoder()
args = []
for _ in range(4):
    tail = tail.lstrip()
    value, end = decoder.raw_decode(tail)
    args.append(value)
    tail = tail[end:].lstrip()
    if tail.startswith(','):
        tail = tail[1:]
_, traces, layout, config = args

def array(value):
    if isinstance(value, dict) and 'bdata' in value:
        data = np.frombuffer(base64.b64decode(value['bdata']), dtype=value['dtype'])
        return data.reshape(tuple(map(int,value['shape'].split(',')))) if 'shape' in value else data
    return np.asarray(value)

counts = [len(array(t['x'])) for t in traces]
assert counts == [1859,338,27,2197,3]
assert len(layout['updatemenus'][0]['buttons']) == 3
assert traces[0]['marker']['reversescale'] is True
for button in layout['updatemenus'][0]['buttons']:
    assert len(button['args'][0]['visible']) == 5
state_path = next(item['path'] for item in manifest['source_files'] if item['role'] == 'state')
with np.load(ROOT / state_path, allow_pickle=False) as state:
    X, x, force = state['X'], state['x'], state['force']
    display = (X + 50 * (x-X)) * 1e6
    for trace in (traces[0], traces[1], traces[2], traces[4]):
        custom = array(trace['customdata'])
        ids = custom[:,0].astype(int)
        np.testing.assert_array_equal(custom[:,1:4], force[ids])
        np.testing.assert_array_equal(np.column_stack([array(trace[c]) for c in ('x','y','z')]), display[ids])
    ids_all = np.concatenate([array(traces[j]['customdata'])[:,0] for j in (0,1)]).astype(int)
    np.testing.assert_array_equal(np.sort(ids_all), np.arange(len(X)))
    np.testing.assert_array_equal(array(traces[0]['marker']['color']), array(traces[0]['customdata'])[:,3]*1e9)

for target in re.findall(r'\]\(([^)]+)\)', (HERE / 'README.md').read_text()):
    if target == 'output/manifests/ARTIFACT_CHECK.json':
        continue  # This report is written after all checks have passed.
    if not target.startswith(('https:', 'http:', '#')):
        assert (HERE / target.split('#')[0]).exists(), target

files = [p for p in OUT.rglob('*') if p.is_file() and p.name != 'ARTIFACT_CHECK.json']
report = dict(status='PASS',git_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
    source_hashes_unchanged=True,source_files_checked=len(manifest['source_files']),export_checks=formats,table_rows=row_counts,
    interactive_trace_counts=counts,interactive_force_vectors_and_display_coordinates_exact=True,interactive_buttons=3,
    interactive_offline_bundle=True,browser_interaction_test='NOT_RUN: browser unavailable; schema and embedded data checked',
    versions=dict(numpy=np.__version__,matplotlib=matplotlib.__version__,plotly=plotly.__version__,pillow=pillow_version),
    new_simulations=0,source_and_result_files_modified=0,
    scripts=[dict(path=p.name,sha256=sha(p)) for p in [HERE/'make_schematics.py',HERE/'verify_artifacts.py']],
    outputs=[dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)])
(OUT / 'manifests/ARTIFACT_CHECK.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ('status','source_files_checked','table_rows','interactive_trace_counts','browser_interaction_test')},indent=2))
