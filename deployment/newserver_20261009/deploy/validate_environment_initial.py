"""Bounded environment verification; no production or P2 campaign is resumed."""
import hashlib, importlib.metadata, json, os, re, shlex, shutil, subprocess, sys, time, traceback
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path('/workspace/brava_cfd');PROJECT=ROOT/'project';OUT=ROOT/'validation'
report=dict(status='RUNNING',start_utc=datetime.now(timezone.utc).isoformat(),scope='Isolated BraVa CFD environment and existing one-step pipe smoke only',acceptance=dict(cpu_gpu_relative_l2=1e-8,single_run_timeout_s=180,production_runs=0,p2_runs=0),checks={})
def save():
    (ROOT/'reports/ENVIRONMENT_VALIDATION.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
def checked(name,value):
    report['checks'][name]=value;save();print(name,json.dumps(value),flush=True)
def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

try:
    import numpy as np
    import scipy,pyvista as pv,vtk,cupy as cp
    build=json.loads((ROOT/'reports/NATIVE_BUILD.json').read_text())
    assert build['status']=='PASS'
    assert Path(sys.prefix)==ROOT/'envs/cfd'
    assert np.__version__=='2.2.6' and scipy.__version__=='1.15.3'
    checked('python_isolation',dict(status='PASS',executable=sys.executable,prefix=sys.prefix,versions={n:importlib.metadata.version(n) for n in ['numpy','scipy','vtk','pyvista','matplotlib','cupy-cuda12x','psutil']}))
    # Exact integer-valued float64 test; CUDA actually executes a generated reduction.
    x=cp.arange(4096,dtype=cp.float64);observed=float(cp.sum(x*x));expected=float(np.sum(np.arange(4096,dtype=np.float64)**2))
    assert observed==expected
    prop=cp.cuda.runtime.getDeviceProperties(0)
    checked('cupy_cuda_execution',dict(status='PASS',device=prop['name'].decode(),compute_capability=[prop['major'],prop['minor']],driver_api=cp.cuda.runtime.driverGetVersion(),cupy_runtime_api=cp.cuda.runtime.runtimeGetVersion(),cpu_value=expected,gpu_value=observed))
    cp.get_default_memory_pool().free_all_blocks()
    from particle_3d.field import FrozenFEMField
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    velocity=points@np.diag([1.,2.,3.])+np.array([.5,.25,.125]);pressure=points@np.array([2.,3.,4.])+5.
    field=FrozenFEMField(points,np.array([[0,1,2,3]]),velocity,pressure)
    sample=field.sample(np.array([.25,.25,.25]))
    assert sample.inside_lumen and np.allclose(sample.velocity_m_s,velocity.mean(axis=0),rtol=0,atol=1e-14)
    checked('project_p1_sampler',dict(status='PASS',test='Existing canonical sampler at centroid of one synthetic tetrahedron; not a production trajectory',velocity=sample.velocity_m_s.tolist()))
    fields=[PROJECT/'mesh/SV_MESH/mesh-complete.mesh.vtu',PROJECT/'cases/balanced_pressure_final/frozen_flow/steady_flow.vtu',PROJECT/'production_workflow/gmsh_reference/candidate_v2/mesh/volume.vtu',PROJECT/'production_workflow/manufactured_pipe/R1_GMSH_V2/cfd_gpu/run/8-procs/result_071.vtu']
    loaded=[]
    for path in fields:
        data=pv.read(path);assert data.n_points>0 and data.n_cells>0 and np.isfinite(data.points).all()
        names=list(data.point_data.keys())
        for name in names:
            if name.lower() in {'velocity','pressure'}:assert np.isfinite(data.point_data[name]).all()
        loaded.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path),points=data.n_points,cells=data.n_cells,point_arrays=names))
    checked('frozen_mesh_and_flow_read',dict(status='PASS',files=loaded))
    petlog=(ROOT/'logs/petsc_check.log').read_text()
    assert all(s in petlog for s in ['successfully with 1 MPI process','successfully with 2 MPI processes','successfully with CUDA'])
    checked('petsc_native_checks',dict(status='PASS',cpu_single=True,cpu_mpi2=True,cuda=True,PETSc='3.25.5',toolkit='12.6',source='logs/petsc_check.log'))
    fixture=OUT/'fixture';manifest=json.loads((fixture/'MANIFEST.json').read_text())
    for row in manifest['files']:assert sha(fixture/row['path'])==row['sha256']
    checked('small_fixture_identity',dict(status='PASS',files=len(manifest['files']),solver_xml_sha256=sha(fixture/'solver.xml')))
    options=shlex.split((PROJECT/'cases/balanced_pressure_final/run/PETSC_OPTIONS.txt').read_text())
    # Same production science/RA/ASM options in both backends. GPU timing flags only on GPU.
    cpu=list(options)
    for flag,value in [('-mat_type','aij'),('-vec_type','standard')]:cpu[cpu.index(flag)+1]=value
    cpu=[v for v in cpu if v!='-log_view_gpu_time']
    run_root=OUT/'native_smoke';assert not run_root.exists(),'Never overwrite previous validation runs'
    run_root.mkdir()
    identity=[];states={}
    for name,ranks,backend in [('cpu1',1,'CPU'),('gpu1',1,'CUDA'),('gpu2',2,'CUDA')]:
        work=run_root/name;work.mkdir();shutil.copytree(fixture/'mesh',work/'mesh')
        for f in ['solver.xml','lumen_inlet.flow']:shutil.copy2(fixture/f,work/f)
        env=dict(os.environ,PETSC_OPTIONS=shlex.join(cpu if backend=='CPU' else options),CUDA_VISIBLE_DEVICES='0')
        command=[str(ROOT/'bin/brava-mpi'),'-n',str(ranks),str(ROOT/'bin/brava-solver'),'solver.xml']
        start=time.monotonic()
        with (work/'solver.log').open('w') as log:
            p=subprocess.run(command,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
        text=(work/'solver.log').read_text(errors='replace')
        assert p.returncode==0,(name,p.returncode)
        assert not re.search(r'DIVERGED_|PETSC ERROR|CUDA error|CUDA_ERROR',text),name
        assert re.search(r'CONVERGED_\w+',text),name
        matrix=sorted(set(re.findall(r'^\s*type:\s*((?:mpi|seq)?aij(?:cusparse)?)\s*$',text,re.M)))
        vector=sorted(set(re.findall(r'^\s*type:\s*((?:mpi|seq)?cuda|seq|mpi)\s*$',text,re.M)))
        if backend=='CUDA':assert any('cusparse' in t for t in matrix) and any('cuda' in t for t in vector),(name,matrix,vector)
        output=work/f'{ranks}-procs/result_001.vtu';assert output.is_file(),name
        data=pv.read(output);keys={k.lower():k for k in data.point_data}
        node=np.asarray(data.point_data[keys['globalnodeid']]);order=np.argsort(node)
        u=np.asarray(data.point_data[keys['velocity']])[order];pvalues=np.asarray(data.point_data[keys['pressure']])[order]
        assert np.isfinite(u).all() and np.isfinite(pvalues).all()
        states[name]=(node[order],u,pvalues)
        entry=dict(name=name,status='PASS',ranks=ranks,backend=backend,seconds=time.monotonic()-start,command=command,PETSC_OPTIONS=env['PETSC_OPTIONS'],matrix_types=matrix,vector_types=vector,n_points=data.n_points,n_cells=data.n_cells,output_sha256=sha(output))
        identity.append(entry);checked('native_'+name,entry)
    comparisons=[]
    for name in ['gpu1','gpu2']:
        assert np.array_equal(states['cpu1'][0],states[name][0])
        row=dict(comparison='cpu1 vs '+name,status='PASS')
        for i,label in [(1,'velocity'),(2,'pressure')]:
            ref=states['cpu1'][i];delta=states[name][i]-ref
            relative=float(np.linalg.norm(delta)/max(np.linalg.norm(ref),np.finfo(float).tiny))
            row[label+'_relative_l2']=relative;row[label+'_max_abs']=float(np.max(np.abs(delta)))
            assert relative<=1e-8,(name,label,relative)
        comparisons.append(row)
    checked('native_cpu_gpu_parity',dict(status='PASS',tolerance_relative_l2=1e-8,comparisons=comparisons))
    (OUT/'NATIVE_SMOKE_IDENTITY.json').write_text(json.dumps(identity,indent=2)+'\n')
    report.update(status='PASS',production_CFD_runs=0,P2_runs=0,particle_runs=0,environment_one_step_smoke_runs=3)
except Exception:
    report.update(status='FAIL',exception=traceback.format_exc());print(report['exception'],flush=True)
finally:
    report['end_utc']=datetime.now(timezone.utc).isoformat();save()
raise SystemExit(report['status']!='PASS')
