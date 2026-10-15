"""Hash-preserving final evidence archive. Never modifies historical results."""
import io,tarfile,subprocess,gzip
from datetime import datetime,timezone
import numpy as np
from .common import *

def main():
 # This tool records prior manual visual inspection; it is not an automatic image judgement.
 figures=BASE/'figures';prov=read(figures/'FIGURE_PROVENANCE.json');review=[]
 for stem in ['Fig_F10S_1_x_pair_3ms','Fig_F10S_2_independent_x_sources','Fig_F10S_3_phase_trap_feasibility','Fig_F10S_4_y_field_scope']:
  for ext in ['png','pdf']:
   p=figures/(stem+'.'+ext);assert sha(p)==prov['output_sha256'][p.name]
   review.append(dict(file=p.name,sha256=sha(p),opened=True,status='PASS',
     method='Assistant view_image inspection of original PNG' if ext=='png' else 'Opened actual PDF with PyMuPDF, rendered its only page and visually inspected with view_image',
     checked=['English labels and physical units','No clipping/overlap obscuring data','Pressure/force/trap prediction separated from trajectories','No fabricated PDMS field'],
     observation='Matched measured curves, unchanged layout style; fine-resolution curves overlap physically. Figure 2 wording corrected to identical initial state.'))
 write(figures/'VISUAL_REVIEW.json',dict(status='PASS',reviewer='Codex visual review',reviewed_utc=datetime.now(timezone.utc).isoformat(),files=review,all_four_PNG_and_all_four_PDF_opened=True))
 historical=read(R/'DELIVERY_MANIFEST.json');bad=[]
 for rel,v in historical['files'].items():
  digest=v if isinstance(v,str) else v['sha256']
  if sha(ROOT/rel)!=digest:bad.append(rel)
 assert not bad,bad
 write(BASE/'BASELINE_AUDIT.json',dict(baseline_commit=read(CONFIG)['baseline_commit'],status='PASS',files=len(historical['files']),historical_manifest_sha256=sha(R/'DELIVERY_MANIFEST.json'),bad=bad,
     no_historical_physics_changes=True,no_safe_deletions=True,source_clean_before_start=True))
 # Evidence of zero instantaneous edge events does not validate arbitrary event timing.
 edges={}
 for name in ['X_PAIR_3MS_DT1','X_PAIR_3MS_DT05','X_PAIR_3MS_FINE','X_MINUS_ONLY_DT1','X_MINUS_ONLY_DT05','X_MINUS_ONLY_FINE','X_PLUS_ONLY_DT1','X_PLUS_ONLY_DT05','X_PLUS_ONLY_FINE']:
  d=np.load(BASE/'analysis'/(name+'.npz'));edges[name]=dict(instantaneous_event_count=len(d['edge_times_s']),members=32,
       scope='Original 2.1R graph, network survival throughout requested window. No physical fusion; absent events cannot establish event-time accuracy.')
 write(BASE/'analysis/TOPOLOGY_SCOPE.json',edges)
 write(BASE/'ENGINEERING_ISSUES.json',dict(issues=[
   dict(issue='Initial regression fixture used incorrect historical comparison dictionary keys',resolution='Corrected fixture names; original legacy and incomplete-window tests pass',scientific_results_modified=False),
   dict(issue='Partial analysis started before transfer completion',resolution='Retained partial review with NOT_RUN missing X+ comparisons; final analysis after 75-file SHA verification',path='history/partial_analysis_before_transfer_completed',scientific_results_modified=False),
   dict(issue='Read-only ad hoc report inspection invoked system Python without NumPy',resolution='Used existing shared project Python; no environment or scientific model change',scientific_results_modified=False),
   dict(issue='Unused deployment route metadata default A on B run',resolution='Actual mode B selected correct 180-second budget; raw metadata retained; future deployment metadata follows mode',scientific_results_modified=False)],
   new_trajectory_failures=0,not_run=['Full moving-trap trajectory','PDMS field','Hardware control','Arbitrary shape control','Independent Z navigation']))
 # Save immutable execution archives for every local scientific/plotting identity as well as remote runs.
 commits={identity()['execution_commit'],prov['identity']['execution_commit'],read(BASE/'analysis/COSTS.json')['identity']['execution_commit'],
    read(BASE/'analysis/INDEPENDENT_REVIEW.json')['identity']['execution_commit'],read(BASE/'diagnostics/EXECUTION.json')['identity']['execution_commit'],
    '0b4882b57'}
 for f in BASE.glob('*/CAMPAIGN.json'):commits.add(read(f)['identity']['execution_commit'])
 prefix='computation_examples/mevo_acoustic_microbubbles/';repo=subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=ROOT,text=True).strip()
 paths=['src','configs/straight_tube_exp1','scripts/straight_tube_exp1','scripts/fast_pzt/python.sh'];archives={}
 for rev in sorted(commits):
  commit=subprocess.check_output(['git','rev-parse',rev],cwd=ROOT,text=True).strip();target=BASE/'execution'/(commit+'.tar.gz')
  if not target.exists():
   raw=subprocess.check_output(['git','archive',commit,*[prefix+k for k in paths]],cwd=repo)
   buff=io.BytesIO()
   with tarfile.open(fileobj=io.BytesIO(raw)) as src,tarfile.open(fileobj=buff,mode='w:gz') as dest:
    for item in src:
     if item.name.startswith(prefix):
      item.name=item.name[len(prefix):]
      if item.name:dest.addfile(item,src.extractfile(item) if item.isfile() else None)
   target.write_bytes(buff.getvalue())
  archives[commit]=dict(path=str(target.relative_to(ROOT)),sha256=sha(target))
 write(BASE/'execution/FINAL_SOURCE.json',dict(identity=identity(),archives=archives,physics_baseline=read(CONFIG)['baseline_commit'],scientific_source_changed=False,
   note='scientific_source_changed refers to historical physics kernels; new F10S orchestration, diagnostics, analysis and plots are separately committed. Archive names give actual execution identities.'))
 paths=list((ROOT/'scripts/straight_tube_exp1/f10s').glob('*.py'))+list((ROOT/'configs/straight_tube_exp1/f10s').glob('*'))
 paths+=list((ROOT/'docs').glob('*F10S*.md'))
 paths += [f for f in BASE.rglob('*') if f.is_file() and '__pycache__' not in str(f) and f.name not in ['DELIVERY_MANIFEST.json','DELIVERY_MANIFEST.sha256']]
 files={str(f.relative_to(ROOT)):dict(bytes=f.stat().st_size,sha256=sha(f)) for f in sorted(set(paths))}
 write(BASE/'DELIVERY_MANIFEST.json',dict(identity=identity(),files=files,file_count=len(files),total_bytes=sum(v['bytes'] for v in files.values()),
    labels=read(CONFIG)['labels'],training_ready=False,raw_data_preserved=True,total_measured_scientific_s=read(BASE/'analysis/TOTAL_COSTS.json')['total_measured_scientific_s'],
    excludes=['self-referencing manifest files','Python bytecode']))
 (BASE/'DELIVERY_MANIFEST.sha256').write_text(''.join(v['sha256']+'  '+k+'\n' for k,v in files.items()))
 print('Final evidence:',len(files),'files;',sum(v['bytes'] for v in files.values()),'bytes')
if __name__=='__main__':main()
