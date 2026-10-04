"""Restore only snapshot-hash-verified missing historical linked artifacts."""
from pathlib import Path
import hashlib
import json
import shutil
import xml.etree.ElementTree as ET
from audit_inventory import AUDIT, PROJECT, ROOT, dump_csv, dump_json

ORIGINAL = Path('/home/lzy/projects/clot_clearance')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    snapshot = json.loads((PROJECT / 'github_sync/SNAPSHOT_MANIFEST.json').read_text())
    identity = {r['path']: r for r in snapshot['files']}
    references = json.loads((AUDIT / 'REFERENCE_SCAN_BEFORE.json').read_text())
    targets = set()
    for item in references['unresolved']:
        if item['kind'] != 'markdown_link':
            continue
        source = ROOT / item['source']
        assert sha(source) == sha(ORIGINAL / source.relative_to(PROJECT)), 'Report differs from original; do not infer identity.'
        target = source.parent / item['token']
        targets.add(str(target.relative_to(PROJECT)))
    # The companion display array is also explicitly named in its original README.
    targets.add('visualization/streaming_fragmentation_demo/streaming_traction_v4/display_samples.npz')
    for name in list(targets):
        if name.endswith('.pvd'):
            for node in ET.parse(ORIGINAL / name).iter('DataSet'):
                targets.add(str(Path(name).parent / node.attrib['file']))
    matched = []
    for name in ['REPORT_ZH.md', 'results/straight_pipe_demo/CONFIG.json', 'results/straight_pipe_demo/states.npz',
                 'runs/fragmentation_pilot_003/CONFIG.json', 'runs/fragmentation_pilot_003/states.npz']:
        actual = sha(PROJECT / name)
        assert actual == sha(ORIGINAL / name)
        matched.append({'path': name, 'sha256': actual})
    rows = []
    for name in sorted(targets):
        source, target = ORIGINAL / name, PROJECT / name
        assert name in identity and identity[name]['sha256'] == sha(source), name
        assert not target.exists(), name
        rows.append(dict(path=str(target.relative_to(ROOT)), source=str(source), size_bytes=source.stat().st_size,
                         hash=identity[name]['sha256'], reason='Pre-existing broken report reference or its VTK/display dependency; restore exact snapshot original, no recomputation.'))
    dump_csv('REFERENCE_RESTORATION_MANIFEST.csv', rows, ['path', 'source', 'size_bytes', 'hash', 'reason'])
    dump_json('REFERENCE_RESTORATION_IDENTITY.json', {'matched_current_inputs': matched, 'original_snapshot': str(PROJECT/'github_sync/SNAPSHOT_MANIFEST.json'),
                                                   'files': len(rows), 'bytes': sum(r['size_bytes'] for r in rows)})
    for row in rows:
        target = ROOT / row['path']
        target.parent.mkdir(parents=True, exist_ok=True)
        with Path(row['source']).open('rb') as src, target.open('xb') as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        assert sha(target) == row['hash']
        row['status'] = 'RESTORED_VERIFIED'
    dump_csv('REFERENCE_RESTORATION_LOG.csv', rows, ['path', 'source', 'size_bytes', 'hash', 'status', 'reason'])
    print('Restored', len(rows), 'files,', sum(r['size_bytes'] for r in rows), 'bytes from verified snapshot originals.')


if __name__ == '__main__':
    main()
