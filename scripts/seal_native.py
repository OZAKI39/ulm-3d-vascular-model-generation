"""Bind the newly built production-source binary to an explicit local profile."""
from pathlib import Path
import datetime, json, subprocess, sys
from human_vessel_flow.paths import ROOT, sha256 as sha, write_json as dump

def main():
    deps=ROOT.parent/'new_vessel_flow_solvers';prefix=deps/'brava_svmultiphysics'
    exe=prefix/'bin/svmultiphysics';petsc=deps/'petsc-3.25.5-cuda12.6/lib/libpetsc.so'
    wrapper=deps/'bin/solver-env'
    sources={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'native/svmultiphysics').rglob('*') if p.is_file()}
    identity=dict(time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=sources,
                  upstream_commit='c3f0bb892b765b718f61069ecd9726dbc6d177fd',
                  scientific_changes=False,P2_experimental_patch_applied=False,
                  executable=str(exe),executable_sha256=sha(exe),petsc_library=str(petsc),petsc_library_sha256=sha(petsc),
                  build_script_sha256=sha(ROOT/'scripts/build_native.sh'),new_CFD_runs=0,
                  prior_shared_binary_differences=['BraVa host-buffer synchronization absent in shared binary source','Reuse header was renamed in shared source'])
    dump(ROOT/'native/BUILD_IDENTITY.json',identity)
    config=json.loads((ROOT/'configs/cfd/native.json').read_text())
    config.update(status='BUILT_PRODUCTION_SOURCE_NOT_CFD_TESTED_THIS_REFACTOR',reason='Native loading and source/binary identity verified; new CFD intentionally not run',
                  executable=str(exe),executable_sha256=sha(exe),petsc_library=str(petsc),petsc_library_sha256=sha(petsc),
                  wrapper=str(wrapper),wrapper_sha256=sha(wrapper),
                  runtime_library_path=f'{deps}/petsc-3.25.5-cuda12.6/lib:{deps}/vtk-9.3/lib:/usr/local/cuda-12.6/lib64:/usr/lib/wsl/lib')
    dump(ROOT/'configs/cfd/native.json',config)
    log=ROOT/'validation/logs';log.mkdir(parents=True,exist_ok=True)
    checks={}
    for name,args in [('native_usage',[str(wrapper),str(exe)]),('native_libraries',[str(wrapper),'ldd',str(exe)]),('mpi_version',['/usr/bin/mpiexec','--version']),('cuda_compiler',['/usr/local/cuda-12.6/bin/nvcc','--version'])]:
        r=subprocess.run(args,capture_output=True,text=True,timeout=30)
        (log/(name+'.log')).write_text(r.stdout+r.stderr)
        if name=='native_usage':
            ok=r.returncode==1 and 'requires the solver input XML file name as an argument' in r.stdout+r.stderr
        elif name=='native_libraries':ok=r.returncode==0 and 'not found' not in r.stdout+r.stderr
        else:ok=r.returncode==0
        checks[name]=dict(status='PASS' if ok else 'FAIL',command=args,returncode=r.returncode)
    dump(ROOT/'validation/NATIVE_ENVIRONMENT.json',dict(checks=checks,CFD_runs=0,native_CFD_accuracy='NOT_RUN'))
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
