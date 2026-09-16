from pathlib import Path
import sys,json,os,time,numpy as np
R=Path(__file__).resolve().parents[1];W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['local_work']);os.environ['PECNUT_CODE']=str(W/'runtime/stokesian_dynamics');sys.path.insert(0,str(R/'src'))
from pecnut_fixed_wall_adapter import native_effective,FixedWall,lattice,drag
C=json.loads((R/'contracts/FIXED_MULTIBLOB_WALL_FEASIBILITY_CONTRACT.json').read_text());free=[];raw=[]
for size,a in C['radii_m'].items():
 for factor in [.5,1,2]:
  mu=.001*factor;r,t=native_effective([[0,0,0]],[a],mu=mu);scale=np.sqrt(drag(a,mu)[:6]);s=r/np.outer(scale,scale);raw.append(r);errT=np.max(abs(s[:3,:3]-np.eye(3)));errR=np.max(abs(s[3:,3:]-np.eye(3)));cross=max(np.max(abs(s[:3,3:])),np.max(abs(s[3:,:3])));assert errT<=1e-4 and errR<=1e-4 and cross<=1e-8
  free.append(dict(size=size,radius_m=a,mu=mu,translation_error=float(errT),rotation_error=float(errR),TR_RT_absolute=float(cross)))
np.savez_compressed(R/'raw/FREE_SPACE_FIRST.npz',matrices=np.array(raw));print('FREE_SPACE_PASS',flush=True)
checks=[];matrices=[]
for layout,layers,spacing in [('HEX',1,1.01),('SQUARE',2,1.01),('HEX',1,1.0)]:
 x,a,layer,meta=lattice(layout,.5,layers,1.55,spacing);target=np.array([.12,.07,1.05]);wall=FixedWall(x,a);r=wall.query([target],[1.],True);xx=np.vstack([target,x]);aa=np.r_[1.,a];native,t=native_effective(xx,aa);D=np.sqrt(drag(1)[:6]);native_scaled=native/np.outer(D,D);err=np.linalg.norm(r['R_scaled']-native_scaled)/np.linalg.norm(native_scaled);assert err<=1e-10,(layout,layers,spacing,err)
 matrices.append((r['R_scaled'],native_scaled));checks.append(dict(layout=layout,layers=layers,spacing_ratio=spacing,N_wall=len(a),minimum_wall_wall_gap=float(min(np.linalg.norm(x[i]-x[j])-a[i]-a[j] for i in range(len(a)) for j in range(i))),cache_native_matrix_error=float(err),wall_factor_seconds=wall.factor_seconds,query_seconds=r['total_seconds'],native_seconds=t['total_seconds'],wall_rotations='PRESCRIBED_ZERO',wall_translations='PRESCRIBED_ZERO',wall_strain='PRESCRIBED_ZERO',wall_wall_lubrication_effect_on_target='exactly zero because all fixed wall kinematics zero; native full-R equality verified'))
 np.savez_compressed(R/'raw'/('NATIVE_CACHE_'+str(len(checks))+'_FIRST.npz'),wall_x=x,wall_radii=a,target=target,R_native=native_scaled,R_cached=r['R_scaled'])
 print('CACHE_NATIVE_PASS',checks[-1],flush=True)
# Native SI radius/viscosity scaling on same two-body geometry, not just algebraic rescaling.
scaling=[];ref=None
for size,a in C['radii_m'].items():
 for mu in [.0005,.001,.002]:
  r,_=native_effective(np.array([[0,0,1.1],[0,0,-.5]])*a,np.array([1,.5])*a,mu=mu);d=np.sqrt(drag(a,mu)[:6]);s=r/np.outer(d,d)
  if ref is None:ref=s
  er=np.linalg.norm(s-ref)/np.linalg.norm(ref);assert er<=1e-10;scaling.append({'size':size,'mu':mu,'scaled_error':float(er)} )
result={'status':'PASS','primary':'Pecnut/stokesian-dynamics','free_space':free,'native_vs_cached':checks,'radius_viscosity_scaling':scaling,'wall_wall_regularization':'np.interp clamps s_dash<2.00001 to first upstream table entry; touching diagnostic only; main spacing1.01 is strictly positive clearance','no_pseudoinverse':True,'upstream_math_modified':False};(R/'validation/TOOL_FREE_SPACE_GATE.json').write_text(json.dumps(result,indent=2)+'\n');print('TOOL_CAPABILITY_PASS')
