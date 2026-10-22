"""Seal files after real trajectories, diagnostics and both animations exist."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];MB=ROOT/'microbubble'
summary=json.loads((MB/'data/final_summary.json').read_text())
assert summary['integration_finished'] and summary['numerical_gate'] and summary['GPU_diagnostics_pass']
for pose in ['0','15']:
    render=json.loads((MB/('candidate_'+pose)/'data/render_manifest.json').read_text())
    assert render['count']==1500 and not render['smoke_only'] and render['text_language']=='English'
files={}
for p in sorted(MB.rglob('*')):
    if not p.is_file() or p.name in ['DELIVERY_COMPLETE.json','LOCAL_VERIFICATION.json'] or '__pycache__' in p.parts or 'logs' in p.parts or 'tracks' in p.parts or 'pilot' in p.parts:continue
    with p.open('rb') as f:files[str(p.relative_to(MB))]=hashlib.file_digest(f,'sha256').hexdigest()
record=dict(status='COMPLETE',flow_sha256=summary['flow_sha256'],files=files,
            tracks='Each COMPLETE.json contains and verifies its own raw trajectory, audit and receipt hashes; verify_collected.py checks all 1500')
(MB/'DELIVERY_COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
