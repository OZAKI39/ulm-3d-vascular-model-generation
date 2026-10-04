"""Compare public author RKCMaterial(:C1) oracle with independent BAQP candidate."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.spatial import cKDTree
from pd_clot.pfpd_candidate import Family,evaluate
O=Path(__file__).resolve().parents[1]/'verification/pfpd_candidate_v1';C=O/'comparisons'
X=np.loadtxt(C/'input_X.csv',delimiter=',');V=np.loadtxt(C/'input_volume.csv',delimiter=',');delta,mu,bulk=np.loadtxt(C/'input_parameters.csv',delimiter=',')
pairs=cKDTree(X).query_pairs(delta,output_type='ndarray');f=Family.from_arrays(X,V,pairs,delta)
rows=[]
for name in ['affine','nonaffine']:
    x=np.loadtxt(C/f'input_x_{name}.csv',delimiter=',');s=evaluate(f,x)
    ids=np.loadtxt(C/f'julia_{name}_ids.csv',delimiter=',',dtype=int)-1
    jf=np.loadtxt(C/f'julia_{name}_force.csv',delimiter=',');jF=np.loadtxt(C/f'julia_{name}_Fbar.csv',delimiter=',').reshape(-1,3,3).transpose(0,2,1)
    je=float(np.loadtxt(C/f'julia_{name}_energy.csv',delimiter=','))
    jb=np.loadtxt(C/f'julia_{name}_Ftilde.csv',delimiter=',');lookup={tuple(p):i for i,p in enumerate(pairs)}
    bid=[lookup[tuple(p.astype(int)-1)] for p in jb[:,:2]];jft=jb[:,2:].reshape(-1,3,3).transpose(0,2,1)
    r=dict(case=name,relative_energy_error=abs(s.energy-je)/abs(je),relative_force_error=float(np.linalg.norm(s.force[ids]-jf)/np.linalg.norm(jf)),Fbar_max_absolute_error=float(abs(s.Fbar[ids]-jF).max()),Ftilde_max_absolute_error=float(abs(s.Ftilde[bid]-jft).max()),python_energy_J=s.energy,julia_energy_J=je)
    assert max(r[k] for k in ['relative_energy_error','relative_force_error','Fbar_max_absolute_error','Ftilde_max_absolute_error'])<1e-9,r
    rows.append(r)
manifest=O/'external/julia_env/Manifest.toml'
result=dict(status='VERIFIED',repository='https://github.com/kaipartmann/Peridynamics.jl',commit='b6e86563845a83b5cb06c805af19a817bde6b106',julia_version='1.10.12',manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),cases=rows,scope='Intact SVK BAQP only; no public PFPD evolution used',point_energy_warning='Author endpoint-normalized and paper symmetric point energies differ, globally sum identically')
(C/'OFFICIAL_COMPARISON.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
# Preserve the small lockfiles outside ignored external binaries/depot.
for name in ['Project.toml','Manifest.toml']:(C/name).write_bytes((O/'external/julia_env'/name).read_bytes())
