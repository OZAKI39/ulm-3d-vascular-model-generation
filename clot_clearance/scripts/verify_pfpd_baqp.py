"""Intact mechanics gate, run before prescribed degradation studies."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import json,copy
from pathlib import Path
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from pd_clot.geometry import make_cloud
from pd_clot.pfpd_candidate import Family,evaluate,Material
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'verification/pfpd_candidate_v1'

def fixture(cells=(4,3,3),ratio=2.01):
    c=json.loads((ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json').read_text())['clot']
    c['cells']=list(cells);c['horizon_ratio']=ratio
    return make_cloud(c)

def deformation(X,nonaffine=False):
    F=np.array([[1.02,.012,.003],[0,.992,.005],[0,0,1.004]])
    x=X@F.T
    if nonaffine:
        q=(X-X.min(axis=0))/(X.max(axis=0)-X.min(axis=0))
        x=x+1.25e-6*np.column_stack([np.sin(2*q[:,1]+q[:,2]),np.cos(3*q[:,0]+q[:,2]),np.sin(q[:,0]+4*q[:,1])])
    return x

def affine_case(arg):
    name,F=arg;cl=fixture((9,9,9),2.01);f=Family.from_cloud(cl);s=evaluate(f,cl.X@np.asarray(F).T)
    interior=np.all((cl.X-cl.X.min(axis=0)>cl.horizon)&(cl.X.max(axis=0)-cl.X>cl.horizon),axis=1)
    err=np.max(abs(s.Fbar-np.asarray(F)),axis=(1,2));bint=interior[f.pairs].all(axis=1)
    berr=np.max(abs(s.Ftilde-np.asarray(F)),axis=(1,2))
    assert interior.any() and bint.any()
    return dict(name=name,bulk_Fbar_error=float(err[interior].max()),boundary_Fbar_error=float(err[~interior].max()),bulk_Ftilde_error=float(berr[bint].max()),boundary_Ftilde_error=float(berr[~bint].max()))

def fd_case(model):
    cl=fixture();f=Family.from_cloud(cl);x=deformation(cl.X,True);mat=Material(model=model)
    state=evaluate(f,x,mat);rows=[]
    for fraction in [1e-4,3e-5,1e-5,3e-6,1e-6]:
        eps=cl.spacing*fraction;num=np.zeros_like(x)
        for k in range(x.size):
            xp=x.copy();xm=x.copy();xp.flat[k]+=eps;xm.flat[k]-=eps
            num.flat[k]=-(evaluate(f,xp,mat).energy-evaluate(f,xm,mat).energy)/(2*eps)
        error=num-state.force
        rows.append(dict(epsilon_m=eps,relative_L2_error=float(np.linalg.norm(error)/np.linalg.norm(state.force)),max_absolute_error_N=float(abs(error).max())))
    assert any(rows[k]['relative_L2_error']<2e-6 and rows[k+1]['relative_L2_error']<2e-6 for k in range(len(rows)-1))
    scale=np.linalg.norm(state.force,axis=1).sum()
    moment=np.cross(x-x.mean(axis=0),state.force).sum(axis=0)
    balance=float(np.linalg.norm(state.force.sum(axis=0))/scale)
    momentrel=float(np.linalg.norm(moment)/(scale*np.max(np.linalg.norm(x-x.mean(axis=0),axis=1))))
    a=.43;R=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])
    rotated=evaluate(f,x@R.T,mat)
    obj=abs(rotated.energy-state.energy)/state.energy
    equiv=np.linalg.norm(rotated.force-state.force@R.T)/np.linalg.norm(state.force)
    trans=evaluate(f,cl.X+[.001,-.002,.003],mat);rotation=evaluate(f,cl.X@R.T,mat)
    assert balance<5e-12 and momentrel<5e-12 and obj<1e-10 and equiv<1e-10
    assert trans.energy<1e-25 and rotation.energy<1e-25
    return dict(material=model,fd=rows,force_balance_relative=balance,moment_balance_relative=momentrel,objectivity_relative=obj,force_rotation_relative=equiv,translation_energy_J=trans.energy,rigid_rotation_energy_J=rotation.energy,translation_force_max_N=float(abs(trans.force).max()),rotation_force_max_N=float(abs(rotation.force).max()))

def main():
    import psutil
    ram=psutil.virtual_memory();workers=max(1,min(4,os.cpu_count(),int((ram.available-.2*ram.total)/(1024**3))))
    cases=[('tension',np.diag([1.02,1,1])),('compression',np.diag([.98,1,1])),('shear',[[1,.03,0],[0,1,0],[0,0,1]]),('isotropic',np.eye(3)*1.01)]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        affine=list(ex.map(affine_case,cases));mechanics=list(ex.map(fd_case,['svk','neo_hookean']))
    maxerr=max(v for r in affine for k,v in r.items() if k!='name');assert maxerr<5e-12
    result=dict(status='VERIFIED',stage='Intact mechanics only; official comparison separately required',affine=affine,mechanics=mechanics,affine_max_error=maxerr,workers=workers,available_ram_before=ram.available)
    (OUT/'baqp_verification/MECHANICS_GATE.json').write_text(json.dumps(result,indent=2))
    cl=fixture();co=OUT/'comparisons'
    np.savetxt(co/'input_X.csv',cl.X,delimiter=',');np.savetxt(co/'input_volume.csv',cl.volume,delimiter=',')
    np.savetxt(co/'input_parameters.csv',[cl.horizon,1000.,50000.],delimiter=',')
    for name,non in [('affine',False),('nonaffine',True)]:np.savetxt(co/f'input_x_{name}.csv',deformation(cl.X,non),delimiter=',')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
