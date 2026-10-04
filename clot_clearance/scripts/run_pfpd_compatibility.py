"""Read-only prescribed historical patterns; no time integration or damage evolution."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import json,csv
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from pd_clot.regularization.equilibrium import make_fixture
from pd_clot.regularization.damage_energy_audit import audit_state,decompose
from pd_clot.pfpd_candidate import Family,evaluate,Material
from pd_clot.pfpd_candidate.degradation import energetic,kinematic
from pd_clot.pfpd_candidate.normalization import c0,coupon,bulk_planar_ratio
ROOT=Path(__file__).resolve().parents[1];O=ROOT/'verification/pfpd_candidate_v1'

def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False))

def hist_task(p):
    p=Path(p);c=json.loads((ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json').read_text())
    cl=make_fixture(c);family=Family.from_cloud(cl);material=Material(model='neo_hookean',mu=c['material']['shear_modulus_Pa'],bulk=c['material']['bulk_modulus_Pa'])
    with np.load(p/'state.npz') as z:x=z['x'].copy();D0=z['D_before'].copy();D1=z['D_after'].copy()
    assert np.all(D1>=D0) and D1.max()<.95
    for a in [x,D0,D1]:a.setflags(write=False)
    old=audit_state(cl,x,c['material'],c['safety'],D=D0);new=audit_state(cl,x,c['material'],c['safety'],D=D1)
    parts,summary=decompose(cl,old,new)
    prior=json.loads((p/'AUDIT.json').read_text())
    err=max(abs(old.global_energy['U_total_J']-prior['U_before_damage_J']),abs(new.global_energy['U_total_J']-prior['U_after_damage_J']))
    assert err<1e-25
    row=dict(case=p.parent.name,jump=p.name,source=str(p.relative_to(ROOT)),D_max=float(D1.max()),replay_error_J=err,
        A_delta_U_J=-summary['delta_U_total_J'],A_reconstruction_delta_U_J=-summary['delta_U_reconstruction_J'],A_softening_delta_U_J=-summary['delta_U_softening_J'],A_stabilization_delta_U_J=-summary['delta_U_stabilization_J'],
        A_max_relative_delta_F=summary['max_relative_F_change'],A_max_relative_delta_K=summary['max_relative_K_change'],
        A_max_absolute_delta_F=float(parts['absolute_F_change'].max()),A_max_absolute_delta_K=float(parts['absolute_K_change'].max()))
    arrays=dict(x=x,D_before=D0,D_after=D1,A_K_before=old.particle['K'],A_K_after=new.particle['K'],A_F_before=old.particle['F'],A_F_after=new.particle['F'])
    states={}
    for name,sc in [('B',1.),('C',.95),('D',0.)]:
        a=evaluate(family,x,material,g=energetic(D0),h=kinematic(D0,sc));b=evaluate(family,x,material,g=energetic(D1),h=kinematic(D1,sc))
        states[name]=(a,b)
        dk=np.linalg.norm(b.K-a.K,axis=(1,2));df=np.linalg.norm(b.Fbar-a.Fbar,axis=(1,2))
        row.update({name+'_delta_U_J':float((b.bond_energy-a.bond_energy).sum()),name+'_U_before_J':a.energy,name+'_U_after_J':b.energy,
                    name+'_max_relative_delta_K':float(np.max(dk/np.linalg.norm(a.K,axis=(1,2)))),name+'_max_relative_delta_F':float(np.max(df/np.linalg.norm(a.Fbar,axis=(1,2)))),
                    name+'_max_absolute_delta_K':float(dk.max()),name+'_max_absolute_delta_F':float(df.max()),name+'_max_delta_Ftilde':float(np.max(abs(b.Ftilde-a.Ftilde)))})
        pairweight=2*family.wbar*cl.volume[cl.pairs[:,0]]*cl.volume[cl.pairs[:,1]]
        row[name+'_softening_delta_U_J']=float(np.sum(pairweight*(energetic(D1)-energetic(D0))*a.W0))
        row[name+'_reconstruction_delta_U_J']=float(np.sum(pairweight*energetic(D1)*(b.W0-a.W0)))
        arrays.update({name+'_K_before':a.K,name+'_K_after':b.K,name+'_F_before':a.Fbar,name+'_F_after':b.Fbar,name+'_Ftilde_before':a.Ftilde,name+'_Ftilde_after':b.Ftilde,name+'_bond_energy_before':a.bond_energy,name+'_bond_energy_after':b.bond_energy})
    for beforeafter in (0,1):
        for key in ('K','Fbar','Ftilde','bond_energy'):
            assert np.array_equal(getattr(states['B'][beforeafter],key),getattr(states['C'][beforeafter],key)), 'STOP: B/C mismatch'
    for name in ('B','C'):
        assert row[name+'_max_absolute_delta_K']==0 and row[name+'_max_absolute_delta_F']==0 and row[name+'_max_delta_Ftilde']==0
        assert row[name+'_delta_U_J']<=0, 'Separated elastic energy increased'
    dest=O/'historical_compatibility'/p.parent.name/p.name;dest.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(dest/'comparison.npz',**arrays);write(dest/'comparison.json',row)
    return row

def coupon_task(arg):return coupon(*arg)

def main():
    import psutil
    assert json.loads((O/'baqp_verification/MECHANICS_GATE.json').read_text())['status']=='VERIFIED'
    mem=psutil.virtual_memory();workers=max(1,min(4,os.cpu_count(),int((mem.available-.2*mem.total)/(1024**3))))
    paths=sorted((ROOT/'results/streaming_erosion_equilibrium_v1').glob('dt_erosion_*/damage_jump_*/state.npz'))
    assert len(paths)==35
    with ProcessPoolExecutor(max_workers=workers) as ex:
        rows=list(ex.map(hist_task,[str(p.parent) for p in paths]))
        constants=list(ex.map(c0,['constant','linear','quadratic']))
        norm=list(ex.map(coupon_task,[(n,m) for n in [20,32,48] for m in [3.01,5.01,8.01]]))
    write(O/'historical_compatibility/ALL_35_JUMPS.json',rows)
    with (O/'historical_compatibility/ALL_35_JUMPS.csv').open('w',newline='') as f:w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
    write(O/'normalization/NORMALIZATION_CONSTANTS.json',constants)
    bulk=[dict(horizon_ratio=m,ratio=bulk_planar_ratio(m)) for m in [3.01,5.01,8.01,16.01,32.01]]
    write(O/'normalization/GRIFFITH_STUDY.json',dict(coupon=norm,bulk_lattice_refinement=bulk,Gc_demo=1.,boundary_correction='none',meaning='NUMERICAL NORMALIZATION TEST; NOT CLOT CALIBRATION'))
    summary=dict(count=len(rows),workers=workers,available_ram_before=mem.available,max_D=max(r['D_max'] for r in rows),replay_max_error_J=max(r['replay_error_J'] for r in rows),models={k:dict(max_relative_delta_F=max(r[k+'_max_relative_delta_F'] for r in rows),max_relative_delta_K=max(r[k+'_max_relative_delta_K'] for r in rows),max_delta_U_J=max(r[k+'_delta_U_J'] for r in rows),min_delta_U_J=min(r[k+'_delta_U_J'] for r in rows),energy_increase_count=sum(r[k+'_delta_U_J']>0 for r in rows)) for k in 'ABCD'})
    write(O/'historical_compatibility/SUMMARY.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
