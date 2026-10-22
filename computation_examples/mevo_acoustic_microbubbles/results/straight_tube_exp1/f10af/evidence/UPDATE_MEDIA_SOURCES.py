from pathlib import Path
import json,hashlib,difflib
p=Path('/home/lzy/projects/microbubble_aggregation_assets/media_sources.json')
b=Path(__file__).resolve().parent;root=b.parent;old=p.read_text();(b/'media_sources_before.json').write_text(old);d=json.loads(old)
for name,source in [('f10af_actual_pattern_evolution.png',next((root/'figures').glob('Fig_F10AF_11*.png'))),('f10af_pattern_attempt.mp4',root/'videos/Video_F10AF_1_stripe_and_pattern_formation.mp4')]:
 dest=p.parent/name;digest=hashlib.sha256(dest.read_bytes()).hexdigest();assert digest==hashlib.sha256(source.read_bytes()).hexdigest();assert not any(r['file']==name for r in d['files']);d['files'].append(dict(file=name,source=str(source),sha256=digest,bytes=dest.stat().st_size,scope='Actual F10AF saved dynamics through 50 ms; no confirmed stripe; exploratory homogeneous-water proxy'))
d['last_updated']='2026-10-19';new=json.dumps(d,ensure_ascii=False,indent=2)+'\n';p.write_text(new);(b/'media_sources_after.json').write_text(new);(b/'media_sources.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='media_sources_before',tofile='media_sources_after')))
