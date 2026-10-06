"""Use the pinned package without modifying/upgrading it; no new mechanics."""
import hashlib,json,os,subprocess,time
from pathlib import Path
import numpy as np
import psutil
from pd_clot.pfpd_candidate.higher_order_rk.oracle_compare import tensor_difference

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1];DEST=OUT/'julia_crosscheck'
EXT=ROOT/'verification/pfpd_candidate_v1/external';PKG=EXT/'Peridynamics.jl'
def dump(name,v):(OUT/name).write_text(json.dumps(v,indent=2)+'\n')
def main():
    commit=subprocess.check_output(['git','-C',str(PKG),'rev-parse','HEAD'],text=True).strip()
    assert commit=='b6e86563845a83b5cb06c805af19a817bde6b106'
    files=[p for p in (PKG/'src').rglob('*') if p.is_file()]+[PKG/'Project.toml',PKG/'Manifest.toml']
    files=[p for p in files if p.exists()];before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    dump('JULIA_PACKAGE_SOURCE_HASHES.json',before)
    overlay=DEST/'depot_overlay';overlay.mkdir(exist_ok=True)
    env=os.environ.copy();env.update(JULIA_DEPOT_PATH=str(overlay)+':'+str(EXT/'julia_depot'),JULIA_PKG_OFFLINE='true',JULIA_NUM_THREADS='1',
        JULIA_NUM_PRECOMPILE_TASKS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    cmd=[str(EXT/'julia-1.10.12/bin/julia'),'--startup-file=no','--compiled-modules=no','--project='+str(EXT/'julia_env'),str(OUT/'julia_pd2_oracle.jl'),str(DEST)]
    start=time.monotonic();peak=0;headroom=1.
    with (DEST/'RUN.log').open('w') as log:
        p=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
        while p.poll() is None:
            try:
                q=psutil.Process(p.pid);peak=max(peak,sum(v.memory_info().rss for v in [q]+q.children(recursive=True)))
            except psutil.Error:pass
            mem=psutil.virtual_memory();headroom=min(headroom,mem.available/mem.total);time.sleep(.1)
    resource=dict(wall_seconds=time.monotonic()-start,peak_process_tree_RSS_bytes=peak,minimum_available_RAM_fraction=headroom,workers=1,exit_code=p.returncode,
        compiled_modules='disabled',old_depot_writes=False,package_installs_or_upgrades=False)
    dump('JULIA_RESOURCE_USAGE.json',resource)
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in before.items())
    if p.returncode:raise RuntimeError('Julia adapter failed; inspect RUN.log before classifying availability')
    a=np.load(DEST/'PYTHON_REPRESENTATIVES.npz');results={};allrows=[]
    for name in ['axial','shear','cross']:
        F=np.loadtxt(DEST/('julia_'+name+'_F.csv'),delimiter=',').reshape(-1,3,3);delta=tensor_difference(a[name],F)
        results[name]=dict(max_F_difference_scaled=delta,PASS=delta<=1e-10)
        for k,id in enumerate(a['ids']):allrows.append(dict(field=name,node=int(id),python_F=a[name][k].tolist(),julia_F=F[k].tolist(),difference_scaled=float(np.linalg.norm(a[name][k]-F[k])/.004)))
    jw=np.loadtxt(DEST/'julia_gradient_weights.csv',delimiter=',');source=a['source'];target=a['target'];g=a['g'];weight_diff=[]
    for row in jw:
        i,j=int(row[0])-1,int(row[1])-1;index=np.flatnonzero((source==i)&(target==j));assert len(index)==1
        weight_diff.append(np.linalg.norm(g[index[0]]-row[2:])*251.25e-6)
    status='PASS' if all(x['PASS'] for x in results.values()) else 'FAIL'
    report=dict(status=status,commit=commit,package_version='0.6.0-DEV',package_path=str(PKG),Julia_version='1.10.12',
        monomial='PD2',package_basis=['x','y','z','x2','xy','xz','y2','yz','z2'],Python_basis=['x','y','z','x2','y2','z2','xy','xz','yz'],
        selectors='Both extract first 3 coordinates and apply 1/delta; both weights include V_j',regularization=dict(epsilon=0.,lambda_=0.,beta=0.),
        source_files_unchanged=len(before),force_calls=0,case_results=results,representative_nodes=a['ids'].tolist(),
        gradient_weights_compared=len(jw),max_gradient_weight_difference_times_delta=float(max(weight_diff)),records=allrows,resources=resource,
        note='Unmodified package uses untruncated SVD reciprocal singular values; primary Python uses direct linear solves. No pseudoinverse repair or kernel change.')
    dump('PERIDYNAMICSJL_CROSSCHECK.json',report)
    assert status=='PASS','Independent package contradicts Python; STOP'
    print(status,results,'weight difference * delta',max(weight_diff),flush=True)

if __name__=='__main__':main()
