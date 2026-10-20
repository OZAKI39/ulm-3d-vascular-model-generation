from pathlib import Path
import json,hashlib,re,os,urllib.parse
R=Path(__file__).resolve().parents[2];TASK='clot_spatially_heterogeneous_variational_lps_v1';RUN='20261017T214917Z';O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()
manifest=json.loads((O/'DELIVERY_MANIFEST.json').read_text())
for p,item in manifest['files'].items():
    q=R/p;assert q.exists(),p;assert q.stat().st_size==item['bytes'] and sha(q)==item['sha256'],p
    if 'symlink_target' in item:assert q.is_symlink() and os.readlink(q)==item['symlink_target']
for p,target in manifest['directory_symlinks'].items():
    q=R/p;assert q.is_symlink() and os.readlink(q)==target and q.resolve().is_dir()
links=0
for p in [F/'OPEN_RESULTS.html',O/'OPEN_RESULTS.html',R/'analysis'/TASK/'OPEN_RESULTS.html']:
    for link in re.findall(r'(?:href|src)="([^"]+)"',p.read_text()):
        if '://' in link or link.startswith('#'):continue
        assert (p.parent/urllib.parse.unquote(link)).resolve().exists(),link;links+=1
doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text())
assert sha(Path(doc['path']))==doc['after_sha256']
entry=json.loads((O/'ENTRYPOINT_UPDATE.json').read_text());assert sha(R/entry['path'])==entry['after_sha256']
result=dict(status='PASS',files_checked=len(manifest['files']),bytes_checked=manifest['total_bytes'],manifest_sha256=sha(O/'DELIVERY_MANIFEST.json'),HTML_links_checked=links,WSL_document_sha256=doc['after_sha256'],entrypoint_sha256=entry['after_sha256'])
(O/'LOCAL_READBACK.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
