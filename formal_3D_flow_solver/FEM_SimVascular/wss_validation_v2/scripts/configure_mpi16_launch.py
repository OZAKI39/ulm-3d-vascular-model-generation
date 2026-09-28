from pathlib import Path
import shutil,json,time,hashlib
V=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
c=V/'stage3/execution_checks/vessel_baseline_cpu16_subLU';old=c.with_name(c.name+'_no_slots');assert not old.exists();assert 'not enough slots' in (c/'run/solver.log').read_text();c.rename(old);c.mkdir();shutil.copytree(old/'SV_MESH',c/'SV_MESH');(c/'run').mkdir()
for name in ['solver.xml','PETSC_OPTIONS.txt']:shutil.copy2(old/'run'/name,c/'run'/name)
shutil.copy2(old/'policy.json',c/'policy.json');shutil.copy2(old/'input_hashes.json',c/'input_hashes.json');(c/'reports').mkdir();shutil.copy2(old/'reports/preregistration.json',c/'reports/preregistration.json')
for case in [c,V/'stage3/vessel_fine']:
 assert not (case/'run/solver.log').exists();oldpol=json.loads((case/'policy.json').read_text());assert oldpol['MPI_ranks']==16
 archive=case/'reports/before_explicit_mpi16_launch';archive.mkdir()
 for name in ['policy.json','input_hashes.json']:shutil.copy2(case/name,archive/name)
 oldpol['MPI_launch_args']=['--oversubscribe','--bind-to','none','--mca','mpi_yield_when_idle','1'];put(case/'policy.json',oldpol)
 hashes=json.loads((case/'input_hashes.json').read_text());hashes['policy.json']=sha(case/'policy.json');put(case/'input_hashes.json',hashes)
 put(case/'reports/explicit_mpi16_launch.json',dict(unix=time.time(),reason='OpenMPI default physical slots fewer than16; explicit oversubscription for smaller subdomains, not more allocated CPU resources. CPU cgroup quota unchanged. No mesh/PDE/BC/tolerance change.',MPI_launch_args=oldpol['MPI_launch_args'],failed_launch_case=str(old)))
 print(case.name,'explicit MPI16 launch recorded')
