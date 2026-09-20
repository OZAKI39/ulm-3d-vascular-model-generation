"""Prepare and execute one explicit authorized run; never hides a failed candidate."""
import json,re,shlex,subprocess,sys
from pathlib import Path
from remote import ROOT,REMOTE,upload,ssh_prefix
R=ROOT/'reports/sv1_3o';C=ROOT/'configs/sv1_3o';S=ROOT/'scripts/sv13o'
name,mode,candidate=sys.argv[1:4];assert mode in ('official','window','profile')
policy=json.loads((C/'policy.json').read_text())
if candidate=='D':options=json.loads((C/'candidate_D.json').read_text())['PETSC_OPTIONS']
else:options=policy['candidate_options'][candidate]
if mode!='official':assert json.loads((R/'OFFICIAL_OUTPUT_STOP_acceptance.json').read_text())['status']=='PASS'
plan=dict(name=name,mode=mode,candidate=candidate,MPI_ranks=1,GPUs=1,OMP_NUM_THREADS=1,start_step=0 if mode=='official' else 60,end_step=20 if mode=='official' else 70,PETSC_OPTIONS=options)
p=C/'runplans'/(name+'.json');p.parent.mkdir(exist_ok=True);assert not p.exists();p.write_text(json.dumps(plan,indent=2)+'\n');upload(p,'configs/runplans/'+p.name)
r=subprocess.run([sys.executable,'-B',S/'invoke.py','flow_remote.py',name])
subprocess.run([sys.executable,'-B',S/'sync_remote.py'],check=True)
ssh=ssh_prefix();case=ROOT/'outputs/sv1_3o'/name;case.mkdir(exist_ok=True)
subprocess.run(['rsync','-a','-e',shlex.join(ssh[:-1]),ssh[-1]+':'+REMOTE+'/outputs/'+name+'/',str(case)+'/'],check=True)
accepted=subprocess.run([sys.executable,'-B',S/'accept_case.py',name])
raise SystemExit(0 if r.returncode==accepted.returncode==0 else 1)
