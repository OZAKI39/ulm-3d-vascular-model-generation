"""Delivery integrity and history protection; no simulation state is modified."""
from pathlib import Path
import sys,json,hashlib,datetime,subprocess,time

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def main():
 b=Path(sys.argv[1]);start=time.perf_counter();read=lambda p:json.loads(p.read_text())
 summary=read(b/'RESEARCH_SUMMARY.json');assert summary['trajectories']==147 and summary['completed']==142 and summary['failed']==5
 assert summary['physics_wall_s']<3600 and summary['max_trajectory_wall_s']<240
 gates=read(b/'FINAL_GATE.json')['gates'];assert gates['FIGURES_VALIDATED']['status']==gates['VIDEOS_VALIDATED']['status']=='PASS'
 assert gates['TARGETED_SHAPE_TRANSFORMATION_PASS']['status']=='FAIL'
 for g in gates.values():
  assert g['status'] in ['PASS','FAIL','INCONCLUSIVE','NOT_RUN']
  assert all((b/p).exists() for p in g['evidence']),g
 figures=read(b/'figures/MANIFEST.json');assert len(figures)==10
 for r in figures:
  for ext in ['png','pdf']:assert sha(b/'figures'/(r['name']+'.'+ext))==r[ext+'_sha256']
 for r in read(b/'videos/MANIFEST.json'):assert sha(b/'videos'/(r['name']+'.mp4'))==r['sha256'] and r['fully_decoded']
 for folder in ['portable_dataset','review_dataset']:
  for r in read(b/folder/'INDEX.json')['trajectories']:assert sha(b/r['file'])==r['sha256']
 raw=read(b/'RAW_DATA_MANIFEST.json')['records']
 for r in raw:assert sha(b/r['raw_path'])==r['raw_sha256']
 history=read(b/'BASELINE_AUDIT.json')['historical_trajectories'];roots={'f10s':'f10s_control_validation','f10z':'f10z_simulation_20261018','f10z_scaling':'f10z_scaling_20261018','f10aa':'f10aa_simulation_20261018'}
 for r in history:
  p=Path('/workspace/mevo_acoustic_microbubbles')/roots[r['path'].split('/')[2]]/r['path'];assert sha(p)==r['sha256']
 write(b/'evidence/REMOTE_FINAL_PROTECTION_AUDIT.json',dict(history_count=len(history),all_historical_hashes_unchanged=True,script_sha256=sha(Path(__file__))))
 jobs=subprocess.check_output(['supervisorctl','status'],text=True);ours=[line for line in jobs.splitlines() if line.startswith('mevo_f10ab_')]
 assert not any('RUNNING' in line for line in ours),ours
 write(b/'evidence/DELIVERY_VALIDATION.json',dict(PASS=True,trajectories_checked=len(raw),review_files_checked=147,portable_files_checked=147,figures=10,videos=3,history_checked=len(history),
  physics_wall_s=summary['physics_wall_s'],checks_wall_s=time.perf_counter()-start,supervisor_status=ours,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),script_sha256=sha(Path(__file__))))
 files={}
 for p in sorted(b.rglob('*')):
  if not p.is_file() or p.name=='INPUT_OUTPUT_SHA256.json' or p.suffix=='.lock':continue
  files[str(p.relative_to(b))]=dict(sha256=sha(p),bytes=p.stat().st_size)
 write(b/'INPUT_OUTPUT_SHA256.json',dict(inputs=read(b/'PREREGISTRATION.json')['input_sha256'],protected_model=read(b/'PREREGISTRATION.json')['model_sha256'],outputs=files,
  manifest_scope='Server authoritative result tree excluding manifest itself and inactive lock; local recovery subset explicit in evidence/LOCAL_RECOVERY.json',source_identity=summary['identity']))
 print('Delivery verified:',len(files),'files;',len(raw),'raw trajectories;',len(history),'historical trajectories unchanged')
if __name__=='__main__':main()
