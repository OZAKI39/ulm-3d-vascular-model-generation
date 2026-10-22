"""Validate full-cohort outputs, rendered files and local HTML links before delivery."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlparse,unquote
import json,sys
from PIL import Image
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import REL,DT,canonical,digest,load_births,completion_matches
R=ROOT/REL;D=R/'data';s=json.loads((D/'final_summary.json').read_text());co=json.loads((D/'FINAL_FORMAL_COHORT.json').read_text());rows=json.loads((D/'formal_metrics.json').read_text())
assert DT==s['dt_s']==.001 and co['count']==len(rows)==s['final_formal_N']>=500
assert co['events']==load_births(ROOT)['events'][:len(rows)]
assert s['core_cohort_sha']=='851c3a4893031334526942fed8b8d430c95245da324eafc7ac7b8600033cfdb9'
regression=json.loads((R/'logs/final_regression/summary.json').read_text())
assert s['old_files_changed']==0 and s['final_tests']['tests']>0
assert regression['mode']=='full' and regression['returncode']==0
assert regression['excluded_full_trajectory_tests']==[]
assert regression['tests']==s['final_tests']['tests']
assert all(s['final_tests'][k]==0 for k in ['failures','errors','skipped'])
bench=json.loads((D/'worker_scaling_results.json').read_text());identity=bench['identity']
for e,r in zip(co['events'],rows):
 folder=R/r['track_relative_path'];assert completion_matches(folder,identity,e)
 meta=json.loads((folder/'trajectory.json').read_text());assert meta['integration_config']['dt_s']==.001
assert sum(s['full_formal'][k] for k in ['O1','O2','O3','stationary','censored','failed'])==len(rows)
visual=json.loads((D/'visualization_manifest.json').read_text());assert visual['all_drawn_particle_ids']==list(range(1,len(rows)+1))
figures=[]
for i in range(1,13 if len(rows)>500 else 12):
 pngs=list((R/'figures').glob(f'Figure_{i:02d}_*.png'));assert len(pngs)==1
 p=pngs[0];pdf=p.with_suffix('.pdf');assert pdf.read_bytes().startswith(b'%PDF-')
 with Image.open(p) as im:
  assert min(im.info.get('dpi',(0,0)))>=299
  im.verify()
 assert digest(p)==visual['output_files'][p.name] and digest(pdf)==visual['output_files'][pdf.name]
 figures.append(dict(png=p.name,pdf=pdf.name,verified_300dpi=True))
animations=json.loads((D/'animation_manifest.json').read_text());assert animations['formal_cohort_sha256']==s['final_formal_cohort_sha']
for a in animations['videos']:
 assert a['full_decode_pass'] and digest(R/a['file'])==a['sha256']
 if a['file'].endswith('formal_overview.mp4'):assert a['all_paths_drawn']==len(rows)
for o in ['O1','O2','O3']:
 expected=min((r['particle_id'] for r in rows if r['status']=='COMPLETED' and r['outlet']==o),default=None)
 assert animations['representatives'][o]==expected
class Links(HTMLParser):
 def __init__(self):super().__init__();self.links=[]
 def handle_starttag(self,tag,attrs):
  self.links.extend(v for k,v in attrs if k in ['href','src'] and v)
parser=Links();parser.feed((R/'OPEN_RESULTS.html').read_text());checked=[]
for link in parser.links:
 u=urlparse(link)
 if u.scheme or u.netloc or not u.path:continue
 p=R/unquote(u.path);assert p.exists(),str(p);checked.append(link)
assert (R/'P9A5_FORMAL_TRAJECTORY_REVIEW_ZH.md').is_file()
receipt=dict(formal_N=len(rows),all_nominal_dt_s=.001,all_formal_hashes_verified=True,cohort_membership_verified=True,core500_unchanged=True,figures=figures,animations_verified=len(animations['videos']),html_local_links_verified=checked,final_tests=s['final_tests'],old_files_changed=0)
(D/'delivery_verification.json').write_bytes(canonical(receipt))
manifest={str(p.relative_to(R)):digest(p) for p in sorted(R.rglob('*')) if p.is_file() and p!=D/'delivery_manifest.json' and '__pycache__' not in p.parts}
(D/'delivery_manifest.json').write_bytes(canonical(manifest))
print(json.dumps(dict(formal_N=len(rows),files=len(manifest),figures=len(figures),animations=len(animations['videos']),all_checks_pass=True),indent=2))
