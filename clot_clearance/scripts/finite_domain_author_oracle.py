"""One new L=0.75 affine cross-implementation check; no historical writer."""
import os,subprocess,time,tempfile
import numpy as np
import psutil
from scripts.finite_domain_boundary_audit import ROOT,OUT,DELTA,material,save

def main():
    ext=ROOT/'verification/pfpd_candidate_v1/external';pin='b6e86563845a83b5cb06c805af19a817bde6b106'
    assert subprocess.check_output(['git','-C',str(ext/'Peridynamics.jl'),'rev-parse','HEAD'],text=True).strip()==pin
    dest=OUT/'author_L075';dest.mkdir(exist_ok=True);m=material()
    with np.load(OUT/'L075/AFFINE_STATE.npz',allow_pickle=False) as s:
        a={k:s[k] for k in s.files}
    for key,data in {'X':a['X'],'x':a['x'],'volume':a['volume'],'parameters':[DELTA,m.mu,m.bulk]}.items():
        np.savetxt(dest/f'input_{key}.csv',data,delimiter=',')
    depot=tempfile.mkdtemp(prefix='clot_boundary_julia_');env=os.environ.copy()
    env.update(JULIA_DEPOT_PATH=depot+':'+str(ext/'julia_depot'),JULIA_NUM_THREADS='1',JULIA_NUM_PRECOMPILE_TASKS='1',JULIA_PKG_OFFLINE='true',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    cmd=[str(ext/'julia-1.10.12/bin/julia'),'--startup-file=no','--project='+str(ext/'julia_env'),str(ROOT/'scripts/finite_domain_author_oracle.jl'),str(dest)]
    start=time.monotonic();peak=swap=0
    with (dest/'RUN.log').open('w') as log:
        p=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
        while p.poll() is None:
            try:
                q=psutil.Process(p.pid);members=[q]+q.children(recursive=True)
                peak=max(peak,sum(q.memory_info().rss for q in members))
            except psutil.NoSuchProcess:pass
            time.sleep(.2)
    resource=dict(wall_seconds=time.monotonic()-start,tree_peak_RSS_bytes=peak,exit_code=p.returncode,workers=1,compiled_cache_overlay=depot,old_depot_writes=False)
    save(dest/'RESOURCES.json',resource)
    assert p.returncode==0,'See author_L075/RUN.log'
    ids=np.loadtxt(dest/'julia_ids.csv',delimiter=',',dtype=int)-1;jf=np.loadtxt(dest/'julia_force.csv',delimiter=',');actual=np.empty_like(jf);actual[ids]=jf
    expected=a['force'];mask=a['full_support'];energy=float(np.loadtxt(dest/'julia_energy.csv',delimiter=','));W=.00902702025
    err=float(np.linalg.norm(expected-actual)/np.linalg.norm(expected));core=float(np.linalg.norm((expected-actual)[mask])/np.linalg.norm(expected[mask]))
    save(dest/'COMPARISON.json',dict(repository='https://github.com/kaipartmann/Peridynamics.jl',commit=pin,case='L075',particles=len(expected),
         strain=.003,phase='DISABLED',material=vars(m),force_relative_error=err,full_support_core_force_relative_error=core,
         energy_relative_error=abs(energy/(W*.00075**3)-1),passed=err<1e-9 and core<1e-8,integrations_run=0,resources=resource))
    assert err<1e-9 and core<1e-8
    print('Pinned author oracle',err,'core',core)

if __name__=='__main__':main()
