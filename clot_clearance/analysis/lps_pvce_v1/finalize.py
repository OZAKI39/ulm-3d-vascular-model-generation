"""Read-only scientific crosschecks and task-only manifests; never starts a solve."""
from analysis.lps_pvce_v1.core import *
import sys,subprocess,time

def main():
 results=ROOT/'results/lps_pvce_energy_release_feasibility_v1/20261011T132103Z';checks=json.loads((BASE/'IMPLEMENTATION_CHECKS.json').read_text());independent=[]
 for j,tip in enumerate([-4,-2,0]):
  z=np.load(BASE/f'A{j}/reference.npz');X=z['X'];b=z['bonds'];q=np.rint(2*X/DX).astype(np.int64);ai=q[b[:,0]];bi=q[b[:,1]];den=bi[:,2]-ai[:,2];num=ai[:,0]*bi[:,2]-ai[:,2]*bi[:,0];cross=ai[:,2]*bi[:,2]<0;cut=cross&np.where(den>0,num<tip*den,num>tip*den);assert np.array_equal(cut,z['alive']==0)
  end=np.load(BASE/f'A{j}/equilibrium.npz');mask=z['mask'];dx=end['x_m']-X;assert np.max(abs(dx[z['high'],2]-1e-6))<1e-19 and np.max(abs(dx[z['low'],2]+1e-6))<1e-19;assert np.array_equal(end['alive'],z['alive']);assert mask.sum()==288
  independent.append(dict(case=f'A{j}',integer_rational_intersection_matches=True,actual_deleted=int(cut.sum()),prescribed_opening_exact=True,final_topology_unchanged=True))
 checks['integer_geometry_crosscheck']=independent;checks['postprocess_only_corrections']='JSON NumPybool cast; empty plottingcollection + layout; document-script string syntax. No solver fix or equilibrium rerun.';save(BASE/'IMPLEMENTATION_CHECKS.json',checks)
 protected=json.loads((BASE/'provenance/PROTECTED_START.json').read_text());changed=[n for n,h in protected.items() if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h];assert not changed;save(BASE/'provenance/PROTECTED_FINAL.json',dict(checked=len(protected),changed=changed,utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
 # Record exact lock identities, state-file shapes and source files for replay.
 required=['REVIEW_LPS_PVCE_ENERGY_RELEASE_FEASIBILITY_V1.md','FINAL_STATUS.json','PVCE_DESIGN_FREEZE.json','PVCE_TOPOLOGY_AUDIT.json','PVCE_EQUILIBRIUM_AUDIT.json','PVCE_ENERGY_GRADIENT_AUDIT.json','PVCE_ENERGY_RELEASE_RESULTS.json','PVCE_CASE_COMPARISON.csv','IMPLEMENTATION_CHECKS.json','OPEN_RESULTS.html']
 assert all((BASE/n).is_file() for n in required);assert len(list(results.glob('Fig*.png')))==4 and len(list(results.glob('Fig*.svg')))==4
 state=json.loads((BASE/'FINAL_STATUS.json').read_text());assert state['DOCUMENT_UPDATE']=='ACTUAL_WSL_UPDATED_AND_READ_BACK'
 # Native-step counter counts force checks, not physical loading steps.
 status=subprocess.run(['supervisorctl','status','y_clot_pvce_20261011T132103Z'],capture_output=True,text=True).stdout
 save(BASE/'provenance/COMPLETION_CHECK.json',dict(required_artifacts=required,four_figures=True,supervisor_status=status,formal_equilibrium_calls=3,equilibrium_retests=0,force_preflight=1,production_unchanged=True,document_receipt=state['DOCUMENT_UPDATE'],native_checks='One initialization run0 + two nonintegrating force evaluations per formal state; separate affine preflight; iteration histories come from actualL-BFGS, no physical-time claim.'))
 source=[]
 for p in sorted((ROOT/'analysis/lps_pvce_v1').glob('*.py')):source.append(dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
 save(BASE/'provenance/TASK_SOURCE_SHA.json',source)
 paths=[p for d in [BASE,results] for p in d.rglob('*') if p.is_file() and p.name not in ['FILE_MANIFEST.json','MANIFEST_VERIFICATION.json','WSL_FINAL_VERIFICATION.json','WSL_COMPACT_COLLECTION.json'] and 'pdf_reader' not in p.parts]
 paths += list((ROOT/'analysis/lps_pvce_v1').glob('*.py'))
 records=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(paths)];save(BASE/'FILE_MANIFEST.json',dict(files=records,total_bytes=sum(q['bytes'] for q in records),exclusions='Self and downstream receipts; local temporarypdf_reader package notscientificartifact. All full states and serverbody.data included.'))
 print('FINALIZED',len(records),'files',sum(q['bytes'] for q in records),'bytes',status.strip())
if __name__=='__main__':main()
