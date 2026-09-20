"""One isolated diagnostic replay; no code patch, no vascular or timing run."""
import shutil,json
from runner_remote import *
w=load('compatibility_winner');sv=load('svmp_gpu_build');old=load('official_fluid_gpu_smoke_execution')
assert old['exit_code']==1 and 'MPI_Comm_rank() function was called after MPI_FINALIZE' in (BASE/old['log']).read_text()
case=BASE/'outputs/official_exit_diagnostic';assert not case.exists();case.mkdir()
shutil.copyfile(BASE/'outputs/official_fluid_gpu_smoke/solver.xml',case/'solver.xml')
# Copy only original official input geometry/boundaries, not generated state.
manifest=json.loads((BASE/'configs/flow_input_manifest.json').read_text())
for f in manifest['files']:
 p=Path(f['destination'])
 if p.parts[0]=='official_fluid_gpu_smoke' and p.name!='solver.xml':
  dst=case/Path(*p.parts[1:]);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(BASE/'outputs'/p,dst)
mpi=json.loads((BASE/'configs/mpi_resolution.json').read_text());wrapper=w['candidate_wrapper']
command=[mpi['working_launcher'],'-n','1',wrapper,'gdb','--batch','-ex','set breakpoint pending on','-ex','set pagination off','-ex','break ompi_mpi_abort','-ex','run','-ex','bt 35','--args',sv['executable'],'solver.xml']
r=run(command,'official_exit_diagnostic',cwd=case,timeout=180,cuda=wrapper,extra_env={'LD_LIBRARY_PATH':sv['runtime_library_path'],'PETSC_OPTIONS':old['PETSC_OPTIONS']})
write('official_exit_diagnosis',dict(status='CAPTURED',run=r,scientific_xml_sha256=digest(case/'solver.xml'),official_input_sha256=old['xml_sha256'],compatibility_repair_applied=False,diagnostic_not_acceptance=True))
