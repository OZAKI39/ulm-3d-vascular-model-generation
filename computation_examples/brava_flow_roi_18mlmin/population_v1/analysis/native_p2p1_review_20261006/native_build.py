"""Build a unit harness, reusing all production objects except main; not a solver rebuild."""
from common import *
import subprocess
import os
import numpy as np

def main():
    assert read(HERE/'SOURCE_IDENTITY.json')['status']=='CONFIRMED'
    old=read(PREVIOUS/'PRODUCTION_SOURCE_IDENTITY.json')
    rt=HERE/'runtime'; rt.mkdir(exist_ok=True)
    compile_cmd=list(old['compile_command'])
    compile_cmd[compile_cmd.index('-c')+1]=str(HERE/'native_fixture.cpp')
    compile_cmd[compile_cmd.index('-o')+1]=str(rt/'native_fixture.o')
    # Read headers from the exact retained production tree, not public upstream.
    old_source=old['source_checkout']
    compile_cmd=[a.replace(old_source,str(RETAINED)) for a in compile_cmd]
    link_cmd=list(old['link_command'])
    main_index=next(i for i,a in enumerate(link_cmd) if a.endswith('/main.cpp.o'))
    link_cmd[main_index]=str(rt/'native_fixture.o')
    link_cmd[link_cmd.index('-o')+1]=str(rt/'native_fixture')
    objects={}
    for a in link_cmd:
        if a.endswith(('.o','.a')) and a!=str(rt/'native_fixture.o'):
            p=Path(old['link_cwd'])/a; objects[str(p.resolve())]=sha(p)
    meta=dict(compile_command=compile_cmd,link_command=link_cmd,cwd=old['link_cwd'],
              reused_objects_sha256=objects,harness_source_sha256=sha(HERE/'native_fixture.cpp'),solver_rebuilt=False)
    put('NATIVE_FIXTURE_BUILD.json',meta)
    subprocess.run(compile_cmd,check=True,cwd=old['link_cwd'])
    subprocess.run(link_cmd,check=True,cwd=old['link_cwd'])
    rng=np.random.default_rng(read(HERE/'PROTOCOL.json')['rng_seed'])
    vertices=np.array([[1,0,0],[0,1,0],[0,0,1],[0,0,0]],float)
    edges=np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]])
    nodes=np.vstack([vertices,vertices[edges].mean(1)])
    points=np.vstack([nodes,rng.dirichlet(np.ones(4),4096)[:,:3]])
    np.savetxt(rt/'native_points.txt',points,fmt='%.17g')
    env=dict(os.environ)
    libdirs=[str(Path(old['PETSc']['path']).parent),'/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/lib']
    env['LD_LIBRARY_PATH']=':'.join(libdirs+[env.get('LD_LIBRARY_PATH','')])
    result=subprocess.run([str(rt/'native_fixture'),str(rt/'native_points.txt'),str(rt/'native_basis.txt'),str(HERE/'NATIVE_UNIT_FIXTURE.json')],env=env,check=True)
    values=np.loadtxt(rt/'native_basis.txt')
    np.savez_compressed(HERE/'NATIVE_BASIS_SAMPLES.npz',xi=points,N=values[:,:10],dN=values[:,10:40].reshape(-1,10,3),
                        P=values[:,40:44],dP=values[:,44:].reshape(-1,4,3),nodes=nodes,edges=edges)
    unit=read(HERE/'NATIVE_UNIT_FIXTURE.json')
    assert unit['status']=='PASS' and unit['native_pressure_edge_handling_pass']
    meta.update(harness_binary_sha256=sha(rt/'native_fixture'),harness_exit_code=result.returncode,
                native_samples_sha256=sha(HERE/'NATIVE_BASIS_SAMPLES.npz'))
    put('NATIVE_FIXTURE_BUILD.json',meta)
    print(unit)

if __name__=='__main__':main()
