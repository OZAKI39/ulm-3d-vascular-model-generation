"""Stdlib-only read-only payload audit; works on both server and WSL copies."""
from pathlib import Path
import json,hashlib,re
from html.parser import HTMLParser
R=Path(__file__).resolve().parents[2]
TASK='lps_native_softening_full_protocol_timestep_admission_v1';RUN='20261013T224838Z'
O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(4194304),b''):h.update(block)
 return h.hexdigest()
def main():
 m=json.loads((O/'PAYLOAD_MANIFEST.json').read_text());assert sha(O/'PAYLOAD_MANIFEST.json')==(O/'PAYLOAD_MANIFEST.sha256').read_text().strip()
 errors=[]
 for rel,info in m['files'].items():
  p=R/rel
  if not p.is_file() or p.stat().st_size!=info['bytes'] or sha(p)!=info['sha256']:errors.append(rel)
 assert not errors,errors
 dep=m['library_dependency'];lib=Path(dep['server_path'])
 if not lib.is_file():lib=R/dep['local_prior_binary_relative']
 assert lib.is_file() and sha(lib)==dep['sha256'],'Reused native library missing or changed'
 st=json.loads((O/'FINAL_STATUS.json').read_text());assert len(st['items'])==30 and st['main_summary']['executed_steps']==2562 and st['full_protocol_safe'] is False
 a=json.loads((O/'IMPLEMENTATION_CHECKS.json').read_text());assert a['formal_cases']==1 and a['all_sparse_d_histories_reconstruct_exactly'] and a['identity']
 freeze=O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.json';assert sha(freeze)==(O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.sha256').read_text().strip()
 for n,h in json.loads(freeze.read_text())['driver_sha256'].items():assert sha(R/'analysis/full_protocol_timestep_v1'/n)==h
 report=(O/'REVIEW_LPS_NATIVE_SOFTENING_FULL_PROTOCOL_TIMESTEP_ADMISSION_V1.md').read_text();assert len(re.findall(r'^### \d+\.',report,flags=re.M))==36
 assert len(list(F.glob('Fig*.png')))==len(list(F.glob('Fig*.svg')))==6
 class Links(HTMLParser):
  def handle_starttag(self,tag,attrs):
   for k,v in attrs:
    if k in ('src','href') and not v.startswith(('http','data:','#')):assert (self.base/v).resolve().is_file(),v
 parser=Links()
 for p in [F/'OPEN_RESULTS.html',O/'OPEN_RESULTS.html']:parser.base=p.parent;parser.feed(p.read_text())
 doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text());assert sha(O/'clot_rupture.before.md')==doc['before_sha256'] and sha(O/'clot_rupture.after.md')==doc['after_sha256']
 localdoc=Path('/home/lzy/projects/clot_rupture.md');doc_current=sha(localdoc)==doc['after_sha256'] if localdoc.exists() else 'NOT_LOCAL_TO_THIS_HOST'
 print(json.dumps(dict(status='PASS',verified_files=len(m['files']),verified_bytes=m['total_bytes'],manifest_sha256=sha(O/'PAYLOAD_MANIFEST.json'),reused_library_verified=True,six_png_six_svg=True,html_links_valid=True,report_questions=36,status_items=30,formal_cases=1,actual_steps=2562,document_current_matches=doc_current),ensure_ascii=False))
if __name__=='__main__':main()
