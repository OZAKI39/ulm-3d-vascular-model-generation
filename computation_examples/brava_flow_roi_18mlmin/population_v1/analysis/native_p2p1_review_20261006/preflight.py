"""Reverify actual retained source and production identities, read live resources."""
from common import *
import datetime
import subprocess

def main():
    old = read(PREVIOUS/'PRODUCTION_SOURCE_IDENTITY.json')
    manifest = read(PREVIOUS/'EXACT_DEPLOYED_SOURCE_MANIFEST.json')
    assert all(sha(RETAINED/n) == h for n,h in manifest.items())
    assert sha(old['actual_binary']) == old['actual_binary_sha256']
    assert sha(old['PETSc']['path']) == old['PETSc']['sha256']
    inputs = {str(p.relative_to(BASE)):sha(p) for p in [
        PREVIOUS/'PRODUCTION_SOURCE_IDENTITY.json', PREVIOUS/'EXACT_DEPLOYED_SOURCE_MANIFEST.json',
        PREVIOUS/'STABILIZATION_SOURCE_AUDIT.json', PREVIOUS/'STABILIZATION_TERMS.json',
        BASE/'mesh/SV_MESH/mesh-complete.mesh.vtu',BASE/'mesh/SV_MESH/mesh-surfaces/WALL.vtp',
        BASE/'cases/balanced_pressure_final/frozen_flow/steady_flow.vtu',
        BASE/'production_workflow/FIELD_VALIDATION.json',BASE/'production_workflow/GPU_MANUFACTURED_RUN_IDENTITY.json',
        BASE/'production_workflow/PARTICLE_CFD_INPUT_CONTRACT.json',
        BASE/'cases/balanced_pressure_final/run/solver.xml']}
    # Previous input manifests already include the original population/source/field identities.
    protected = {}
    for rel,root in [('population_v1/analysis/field_reliability_event_review/INPUT_SHA256.json',BASE),
                     ('population_v1/analysis/field_reliability_event_review/OUTPUT_SHA256.json',HERE.parent/'field_reliability_event_review')]:
        entries=read(BASE/rel)
        assert all(sha(root/k)==v for k,v in entries.items()), rel
        inputs[rel]=sha(BASE/rel)
        protected[rel]={'entries':len(entries),'verified':True}
    resources = {k:(Path('/sys/fs/cgroup')/k).read_text().strip() for k in
                 ['cpu.max','cpuset.cpus.effective','memory.max','memory.current','memory.stat','memory.events']}
    for k,cmd in [('GPU',['nvidia-smi','--query-gpu=name,driver_version,memory.total,memory.used,memory.free,utilization.gpu','--format=csv']),
                  ('GPU_processes',['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv']),
                  ('disk',['df','-h',str(BASE)]),('CUDA',['nvidia-smi'])]:
        resources[k]=subprocess.check_output(cmd,text=True)
    git=lambda *args: subprocess.check_output(['git','-c',f'safe.directory={BASE}','-C',str(BASE),*args],text=True).strip()
    resources.update(UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),project_root=git('rev-parse','--show-toplevel'),
                     branch=git('branch','--show-current'),HEAD=git('rev-parse','HEAD'),preexisting_status=git('status','--short'))
    put('RESOURCES.json',resources)
    put('INPUT_SHA256.json',inputs)
    put('SOURCE_IDENTITY.json',dict(status='CONFIRMED',source_commit=old['source_commit'],
        retained_source=str(RETAINED),verified_files=len(manifest),binary=old['actual_binary'],binary_sha256=old['actual_binary_sha256'],
        PETSc=old['PETSc'],build_config=old['actual_compiler_flags'],production_modified=False,protected=protected))
    tests=[]
    for p in (RETAINED/'tests').rglob('*'):
        if p.is_file() and p.suffix in ['.xml','.cpp','.h','.py','.md']:
            t=p.read_text(errors='replace')
            if 'taylor_hood' in t.lower() or 'TET10' in t or 'Tetra10' in t:
                tests.append(dict(path=str(p.relative_to(RETAINED)),sha256=sha(p)))
    put('NATIVE_REGRESSION_INVENTORY.json',dict(matches=tests,executed=False,
        limitation='Retained Taylor-Hood XML fixtures are Stokes/solid; no retained 3D fluid TET10 GPU regression execution claimed.'))
    print('SOURCE_CONFIRMED',len(manifest),'historical manifests verified',protected)

if __name__=='__main__': main()
