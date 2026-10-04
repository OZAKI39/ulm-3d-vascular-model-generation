"""Pinned author-package comparison with new inputs and writable-depot overlay."""
from scripts.spatial_audit_v2 import *
import subprocess,hashlib

def main():
    old=ROOT/'verification/pfpd_candidate_v1';ext=old/'external';rows=[]
    for n in LEVELS:
        c,f,mat,top=fixture(n);folder=OUT/f'external_reference/C{n}';folder.mkdir(exist_ok=True)
        for name,a in dict(X=f.X,volume=f.volume,parameters=[DELTA,mat.mu,mat.bulk]).items():
            np.savetxt(folder/f'input_{name}.csv',a,delimiter=',')
        x=f.X@np.array([[1.004,.002,0],[0,1.003,.001],[0,0,1.01]]).T
        rng=np.random.default_rng(438);xn=x+L*1e-4*rng.normal(size=x.shape)
        for name,a in [('affine',x),('nonaffine',xn)]:np.savetxt(folder/f'input_x_{name}.csv',a,delimiter=',')
        shutil.copyfile(old/'comparisons/official_oracle.jl',folder/'official_oracle.jl')
        env=os.environ.copy();env['JULIA_DEPOT_PATH']=str(OUT/'external_reference/depot')+':'+str(ext/'julia_depot');env['JULIA_NUM_THREADS']='1';env['JULIA_PKG_OFFLINE']='true'
        cmd=[str(ext/'julia-1.10.12/bin/julia'),'--startup-file=no','--project='+str(ext/'julia_env'),str(folder/'official_oracle.jl')]
        start=time.time()
        with (folder/'JULIA.log').open('w') as log:proc=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
        if proc.returncode:raise RuntimeError(f'Julia C{n}: see {folder}/JULIA.log')
        lookup={tuple(p):i for i,p in enumerate(f.pairs)}
        for name,xx in [('affine',x),('nonaffine',xn)]:
            state=evaluate(f,xx,mat);ids=np.loadtxt(folder/f'julia_{name}_ids.csv',delimiter=',',dtype=int)-1
            jf=np.loadtxt(folder/f'julia_{name}_force.csv',delimiter=',');jF=np.loadtxt(folder/f'julia_{name}_Fbar.csv',delimiter=',').reshape(-1,3,3).transpose(0,2,1)
            je=float(np.loadtxt(folder/f'julia_{name}_energy.csv',delimiter=','));jb=np.loadtxt(folder/f'julia_{name}_Ftilde.csv',delimiter=',')
            bid=[lookup[tuple(p.astype(int)-1)] for p in jb[:,:2]];jft=jb[:,2:].reshape(-1,3,3).transpose(0,2,1)
            errors=dict(energy=abs(state.energy-je)/abs(je),force=float(np.linalg.norm(state.force[ids]-jf)/np.linalg.norm(jf)),
                Fbar=float(abs(state.Fbar[ids]-jF).max()),Ftilde=float(abs(state.Ftilde[bid]-jft).max()))
            rows.append(dict(level=f'C{n}',case=name,errors=errors,passed=max(errors.values())<1e-9,wall_seconds=time.time()-start))
        print(f'Julia C{n} verified={all(r["passed"] for r in rows)}',flush=True)
        save(OUT/'external_reference/COMPARISON.json',dict(status='VERIFIED' if len(rows)==6 and all(r['passed'] for r in rows) else 'IN_PROGRESS',
            commit=subprocess.check_output(['git','-C',str(ext/'Peridynamics.jl'),'rev-parse','HEAD'],text=True).strip(),
            repository='https://github.com/kaipartmann/Peridynamics.jl',julia_version='1.10.12',license='MIT',
            manifest_sha256=hashlib.sha256((ext/'julia_env/Manifest.toml').read_bytes()).hexdigest(),cases=rows))
    shutil.copyfile(old/'resources/PERIDYNAMICS_LICENSE.txt',OUT/'external_reference/PERIDYNAMICS_LICENSE.txt')

if __name__=='__main__':main()
