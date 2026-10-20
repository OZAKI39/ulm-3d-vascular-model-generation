"""WSL standard-library readback of the sealed, selectively synchronized delivery."""
from pathlib import Path
import hashlib
import json
import os
import re

root = Path(__file__).resolve().parents[2]
task = 'clot_multifragment_crushing_mechanism_v1'
run = '20261016T222934Z'
out = root / 'verification' / task / run
fig = root / 'results' / task / run


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(4 * 1024**2), b''):
            digest.update(data)
    return digest.hexdigest()


def read(path):
    return json.loads(path.read_text())


manifest = read(out / 'DELIVERY_MANIFEST.json')
assert manifest['status'] == 'SEALED'
for item in manifest['files']:
    path = root / item['path']
    assert path.is_file(), path
    if item['kind'] == 'symlink':
        assert path.is_symlink() and os.readlink(path) == item['target'], path
        assert sha(path) == item['target_sha256'], path
    else:
        assert path.stat().st_size == item['bytes'], path
        assert sha(path) == item['sha256'], path

required = [
    'SPECIMEN_AND_MATERIAL_FREEZE.json', 'INDENTER_GEOMETRY.json',
    'INDENTER_LOADING_PROTOCOL.json', 'HETEROGENEITY_FIELD.npz',
    'HETEROGENEITY_FREEZE.json', 'UNIFORM_CRUSHING_RESULTS.json',
    'HETEROGENEOUS_CRUSHING_RESULTS.json', 'FRACTURE_EVENTS.json',
    'CONNECTED_COMPONENT_HISTORY.json', 'MULTIFRAGMENT_SIZE_DISTRIBUTION.csv',
    'POST_UNLOAD_FRAGMENT_STABILITY.json', 'DYNAMIC_CONTACT_AND_ENERGY_LEDGER.json',
    'RUN_RESOURCES.json', 'FINAL_STATUS.json',
    'REVIEW_CLOT_MULTIFRAGMENT_CRUSHING_MECHANISM_V1.md',
    'REPRODUCE.md', 'OPEN_RESULTS.html', 'SOURCE_SHA256.json',
    'SOURCE_DIFF.patch', 'ENVIRONMENT.json', 'DELIVERY_VALIDATION.json',
    'WSL_DOCUMENT_UPDATE.json',
]
assert all((out / name).is_file() for name in required)
assert read(out / 'DELIVERY_VALIDATION.json')['status'] == 'PASS'
assert len(read(out / 'FINAL_STATUS.json')) == 10
report = (out / 'REVIEW_CLOT_MULTIFRAGMENT_CRUSHING_MECHANISM_V1.md').read_text()
assert len(re.findall(r'^\d+\. \*\*', report, re.M)) == 16
figures = read(fig / 'FIGURE_MANIFEST.json')['figures']
assert len(figures) == 10
for item in figures:
    for extension in ['.png', '.svg']:
        assert (fig / (item['name'] + extension)).is_file()
links = re.findall(r'(?:href|src)="([^"]+)"', (fig / 'OPEN_RESULTS.html').read_text())
for link in links:
    assert not link.startswith(('http:', 'https:')), link
    assert (fig / link).resolve().is_file(), link
document = read(out / 'WSL_DOCUMENT_UPDATE.json')
assert document['history_preserved_exactly'] and document['readback_equal']
assert sha(Path(document['path'])) == document['after_sha256']
receipt = dict(
    status='PASS', root=str(root), checked_manifest_files=len(manifest['files']),
    manifest_sha256=sha(out / 'DELIVERY_MANIFEST.json'),
    checked_html_links=len(links), figure_pairs=len(figures), scientific_answers=16,
    formal_cases=2, final_status_keys=10,
    wsl_document_sha256=document['after_sha256'],
    server_only_builds=manifest['server_only_builds'],
    scope='Original trajectories and evidence synchronized; reproducible build binaries excluded.',
)
(out / 'LOCAL_READBACK.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n')
print('LOCAL_READBACK_PASS', len(manifest['files']), 'files;', len(figures), 'figure pairs')
