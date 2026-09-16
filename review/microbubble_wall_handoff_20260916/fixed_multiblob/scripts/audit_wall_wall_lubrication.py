from pathlib import Path
import sys,os,json,time
import numpy as np
R=Path(__file__).resolve().parents[1];W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['local_work']);os.environ['PECNUT_CODE']=str(W/'runtime/stokesian_dynamics');sys.path.insert(0,str(R/'src'))
from pecnut_fixed_wall_adapter import FixedWall,lattice,posdata,index,drag,dense,generate_grand_resistance_matrix
rows=[]
for spacing in [1.,1.01,1.05]:
 x,a,ids,meta=lattice('HEX',.5,1,1.55,spacing);p=np.array([.13,.07,1.1]);xx=np.vstack([p,x]);aa=np.r_[1.,a];N=len(aa);t=time.perf_counter();full,_,_,_=generate_grand_resistance_matrix(posdata(xx,aa),[],regenerate_Minfinity=True,cutoff_factor=2,mu=1.);full=dense(full);ix=np.concatenate([index(N,i) for i in range(N)]);bulk=np.concatenate([drag(v) for v in aa]);scaled=full[np.ix_(ix,ix)]/np.sqrt(np.outer(bulk,bulk));ev=np.linalg.eigvalsh((scaled+scaled.T)/2);wall=FixedWall(x,a);q=wall.query(p,[1.]);target=scaled[:6,:6];er=np.linalg.norm(q['R_scaled']-target)/np.linalg.norm(target);assert er<=1e-10
 np.savez_compressed(R/f'raw/WALL_WALL_LUBRICATION_S{spacing:g}.npz',coords=xx,radii=aa,full_native_scaled=scaled,target_native=target,target_cached=q['R_scaled'],wall_M_scaled=wall.C)
 dd=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2);np.fill_diagonal(dd,np.inf);gap=float(dd.min()-1)
 rows.append(dict(spacing_ratio=spacing,N_wall=len(a),minimum_wall_wall_gap_over_a=gap,full_R_scaled_condition=float(ev[-1]/ev[0]),full_R_min_eigenvalue=float(ev[0]),wall_M_scaled_condition=float(np.linalg.cond(wall.C)),target_R_scaled_condition=float(np.linalg.cond(target)),native_cached_target_error=float(er),seconds=time.perf_counter()-t,wall_wall_kinematics='exact zero',nearfield_clamp_s_dash=2.00001,farfield_contact_protection_s_dash=2.001,formal_test_spacing=spacing==1.01))
(R/'validation/WALL_WALL_LUBRICATION_DIAGNOSTIC.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
