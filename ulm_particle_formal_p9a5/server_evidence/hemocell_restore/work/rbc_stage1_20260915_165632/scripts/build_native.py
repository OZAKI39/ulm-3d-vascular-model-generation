from pathlib import Path
import subprocess,json,time,hashlib,os,traceback
R=Path(__file__).resolve().parents[1]
def save(n,v):(R/n).write_text(json.dumps(v,indent=2)+'\n')
def cmd(name,c):
 save('BUILD_STATE.json',dict(status='RUNNING',phase=name,unix=time.time()))
 with (R/'logs'/f'{name}.log').open('w') as f:p=subprocess.run(c,stdout=f,stderr=subprocess.STDOUT)
 save('provenance/'+name+'_command.json',dict(command=c,returncode=p.returncode));assert p.returncode==0,name
with (R/'BUILD_NATIVE_STARTED.json').open('x') as f:json.dump(dict(unix=time.time(),automatic_retries=0),f)
try:
 os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',HWLOC_COMPONENTS='-opencl')
 cmd('configure',['cmake','-S',str(R/'source'),'-B',str(R/'build'),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_COMPILER=/usr/bin/g++','-DCMAKE_C_COMPILER=/usr/bin/gcc'])
 cmd('compile',['cmake','--build',str(R/'build'),'--parallel','2','--verbose'])
 def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
 prov=dict(status='PASS',build_scope='Complete HemoCell CPU/MPI library and native RBC geometry construction executable; timestep entry point pending geometry gate',hemocell_commit='5a410848bd5c57d5ae1c171112e78eab4a82e650',palabos_commit='05712164d940a42e06afdd705249912fa0c49f14',new_server_build=True,interior_viscosity=False,GPU_IBM=False,compiler=subprocess.check_output(['g++','--version'],text=True),MPI=subprocess.check_output(['mpirun','--version'],text=True),CMake=subprocess.check_output(['cmake','--version'],text=True),HDF5=subprocess.check_output(['h5pcc','-showconfig'],text=True),library_sha256=sha(R/'build/libhemocell_cpu.a'),geometry_binary_sha256=sha(R/'build/rbc_geometry_probe'),driver_sha256=sha(R/'source/rbc_geometry_probe.cpp'),flags=(R/'build/CMakeFiles/hemocell_cpu.dir/flags.make').read_text(),link_command=(R/'build/CMakeFiles/rbc_geometry_probe.dir/link.txt').read_text())
 save('RBC_STAGE1_BUILD_PROVENANCE.json',prov);save('BUILD_STATE.json',dict(status='PASS',phase='NATIVE_CPU_BUILD_COMPLETE',solver_steps=0,unix=time.time()))
except Exception as e:save('BUILD_STATE.json',dict(status='FAIL',error=str(e),traceback=traceback.format_exc(),unix=time.time()));raise
