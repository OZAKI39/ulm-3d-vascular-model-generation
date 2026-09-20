"""One authorized run, guarded by immutable plans and prior health gates."""
import json,subprocess,sys
from remote import ROOT,upload
R=ROOT/'reports/sv1_3q';C=ROOT/'configs/sv1_3q';S=ROOT/'scripts/sv13q'
candidate,mode=sys.argv[1:3];assert candidate in ('R2','R3','R5','RA');assert mode in ('smoke','window','early')
if mode!='smoke':assert json.loads((R/(candidate+'_SMOKE_acceptance.json')).read_text())['status']=='PASS'
if mode=='early':
 selection=json.loads((R/'late_ranking.json').read_text());assert candidate in selection['top_two']
 assert json.loads((R/(candidate+'_WINDOW_acceptance.json')).read_text())['status']=='PASS'
name=candidate+'_'+mode.upper();policy=json.loads((C/'policy.json').read_text())
options=policy['baseline_options']+(' -sv_pc_adaptive_rebuild true' if candidate=='RA' else ' -sv_pc_rebuild_interval '+candidate[1:])
start,end={'smoke':(0,3),'window':(60,70),'early':(10,20)}[mode]
plan=dict(name=name,mode=mode,candidate=candidate,MPI_ranks=1,GPUs=1,OMP_NUM_THREADS=1,start_step=start,end_step=end,PETSC_OPTIONS=options,build_report='svmp_reuse_build')
p=C/'runplans'/(name+'.json');p.parent.mkdir(exist_ok=True);assert not p.exists();p.write_text(json.dumps(plan,indent=2)+'\n');upload(p,'configs/runplans/'+p.name)
native=subprocess.run([sys.executable,'-B',S/'invoke.py','flow_remote.py',name])
subprocess.run([sys.executable,'-B',S/'fetch_case.py',name],check=True)
accept=subprocess.run([sys.executable,'-B',S/'accept_case.py',name])
raise SystemExit(0 if native.returncode==accept.returncode==0 else 1)
