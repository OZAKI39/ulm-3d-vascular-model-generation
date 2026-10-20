"""Offline prescribed opening of the FULL original endpoint-star model; not Gc."""
from common import *
from native import predict
import csv,time
start=time.perf_counter();ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'HETEROGENEITY_FIELD.npz');model=energy_model(ref);bulk=ref['bond_types']==1;cut=ref['initial_cut'];b=ref['bonds'];sign=np.where(ref['X'][:,2]>0,1.,-1.);cross=(sign[b[:,0]]!=sign[b[:,1]])&~cut
on=S0[ref['bond_types']];base=SC[ref['bond_types']];cases={'uniform':base.copy(),'spatial_actual':field['bond_sc'].copy()}
for name,value in [('minimum',field['bond_sc'][bulk].min()),('median',np.median(field['bond_sc'][bulk])),('maximum',field['bond_sc'][bulk].max())]:
 cases[name]=base.copy();cases[name][bulk]=value
rows=[];summaries={}
for name,off in cases.items():
 k=off*cut;d=cut.astype(float);D=0.;W=0.;lastF=0.;lastu=0.;first=None;initial_slope=None;minimum_drop=0.
 for u in np.linspace(0,180e-6,721):
  x=ref['X'].copy();x[:,2]+=.5*u*sign
  k,dnext=predict(model,ref,x,k,on,off);old=model.evaluate(x,d);new=model.evaluate(x,dnext);drop=old['U']-new['U'];minimum_drop=min(minimum_drop,drop);assert drop>=-1e-20
  D+=drop;d=dnext;f=-.5*np.dot(new['force'][:,2],sign);W+=.5*(lastF+f)*(u-lastu)
  if u>0 and initial_slope is None:initial_slope=f/u
  if first is None and np.any((d>0)&~cut):first=dict(opening_m=u,force_N=f)
  rows.append(dict(case=name,opening_m=u,reaction_N=f,U_J=new['U'],D_J=D,work_trapezoid_J=W,offline_quadrature_residual_J=W-new['U']-D,failed_crossing_bonds=int(np.sum((d==1)&cross)),total_crossing_bonds=int(cross.sum())))
  lastF=f;lastu=u
 summaries[name]=dict(first_softening_grid_sample=first,initial_stiffness_N_m=initial_slope,path_dissipation_J=D,final_stored_U_J=new['U'],complete_crossing_bond_count=int(np.sum((d==1)&cross)),crossing_bonds=int(cross.sum()),minimum_increment_J=minimum_drop,final_quadrature_residual_J=W-new['U']-D)
 assert np.all(d[cross]==1),'Opening path did not separate the prescribed plane'
assert len(set(v['initial_stiffness_N_m'] for v in summaries.values()))==1
with (O/'PRESCRIBED_PATH_RESPONSE.csv').open('w') as f:wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
dump(O/'PRESCRIBED_PATH_ENERGY.json',dict(status='PASS',scope='OFFLINE_FULL_VARIATIONAL_ENERGY_PRESCRIBED_MODE_I_PATH_NOT_NATIVE_DYNAMICS_NOT_GC_EFFECTIVE',path='Opposite rigid translations of original nodes above/below reference z=0, displacement = +/- opening/2; original precut retained, no new initial bond deletion',summaries=summaries,wall_s=time.perf_counter()-start,negative_dissipation=False,bulk_onset=.15,initial_elastic_response_exact=True,interpretation='Only discrete prescribed-path energy; sc/Delta_s are not proportional to calibrated toughness. The three uniform representative sc cards are offline diagnostics, not formal shear cases.'))
print('PATH_AUDIT_PASS',summaries,flush=True)
