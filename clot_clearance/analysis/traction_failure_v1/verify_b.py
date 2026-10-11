"""Independent artifact reconstruction and deterministic offline interface checks."""
from analysis.traction_failure_v1.common import *
from analysis.traction_failure_v1.map_b import apply_nodal_forces,constraints
from scipy.sparse import csr_matrix
s=np.load(PREV/'refined/CLOT_SURFACE_TRACTIONS.npz');d=np.load(B/'CLOT_NODAL_TRACTION_FORCES.npz');m=np.load(B/'CFD_TO_LAMMPS_MAPPING_MATRIX.npz');W=csr_matrix((m['data'],m['indices'],m['indptr']),shape=tuple(m['shape']));X=d['reference_coordinates_m'];x0=d['reference_point_m'];area=s['triangle_area_m2'];ix=m['loaded_node_indices'];L=float(m['length_scale_m']);C6=constraints(X[ix],x0,L);checks={};ind={}
for key,field in [('pressure','pressure_traction_Pa'),('shear','tangential_traction_Pa'),('viscous_normal','viscous_normal_traction_Pa')]:
 f0=W@(s[field]*area[:,None]);f=f0.copy();f[ix]+=(np.repeat(m[key+'_D'],3)*(C6.T@m[key+'_lambda'])).reshape(-1,3);err=np.linalg.norm(f-d[key+'_force_N'])/np.linalg.norm(f);checks[key+'_saved_reconstruction']=err<1e-13
 # Direct exact triangle integral: linear p barycentric first/second moments.
 if key=='pressure':
  bary=np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]]);pnode=np.linalg.solve(bary,s['quadrature_pressure_Pa'].T).T;xyz=s['triangles_xyz_m'];int_xp=area[:,None]/12*(xyz.sum(1)*pnode.sum(1)[:,None]+np.einsum('fij,fi->fj',xyz,pnode));moment=np.cross(int_xp-x0*(area*pnode.mean(1))[:,None],-s['n_clot']).sum(0)
 else:moment=np.cross(s['surface_coordinates_m']-x0,s[field]*area[:,None]).sum(0)
 F=(s[field]*area[:,None]).sum(0);scale=np.linalg.norm(s[field]*area[:,None],axis=1).sum();fe=np.linalg.norm(f.sum(0)-F)/scale;me=np.linalg.norm(np.cross(X-x0,f).sum(0)-moment)/(scale*L);checks[key+'_force_moment']=max(fe,me)<1e-8;ind[key]=dict(force_error_scaled=fe,moment_error_scaled=me,artifact_relative_error=err)
 # Translation of torque center, an independent sign / lever-arm test.
 x1=x0+np.array([.0002,-.0003,.0001]);newM=np.cross(X-x1,f).sum(0);checks[key+'_torque_translation']=np.linalg.norm(newM-(moment-np.cross(x1-x0,F)))/(scale*L)<1e-8
ids=np.array([701,9,88]);forces=np.array([[1.,2.,3.],[-.4,.5,.6]]);density=np.arange(9.).reshape(3,3)*1e9;volume=np.array([2e-13,3e-13,4e-13]);nodal=density*volume[:,None];out=apply_nodal_forces(ids,nodal,[9,701],forces,[88]);expected=nodal.copy();expected[1]+=forces[0];expected[0]+=forces[1];checks['density_conversion_once_then_N']=np.array_equal(out,expected)
try:apply_nodal_forces(ids,nodal,[88],forces[:1],[88]);checks['reject_anchor']=False
except AssertionError:checks['reject_anchor']=True
writer=(ROOT/'pd_clot/lammps_fragmentation/lammps_writer.py').read_text();checks['production_source_order']=writer.index('{force_conversion()}')<writer.index('fix stream surface addforce');checks['raw_density_source']=all(x in writer for x in ['fx*d_vfrac','fy*d_vfrac','fz*d_vfrac'])
regions={};edges=np.quantile(X[ix,0],[0,.25,.5,.75,1])
for name in ['pressure','shear','total']:
 field={'pressure':'pressure_traction_Pa','shear':'tangential_traction_Pa','total':'total_traction_Pa'}[name];ff=s[field]*area[:,None];fc=d[name+'_force_N'];fi=d[name+'_initial_force_N'];rows=[]
 for k in range(4):
  lo=-np.inf if k==0 else edges[k];hi=np.inf if k==3 else edges[k+1];sf=(s['surface_coordinates_m'][:,0]>=lo)&(s['surface_coordinates_m'][:,0]<hi);sn=(X[:,0]>=lo)&(X[:,0]<hi)
  rows.append(dict(x_interval_m=[None if k==0 else lo,None if k==3 else hi],CFD_force_N=ff[sf].sum(0),node_initial_force_N=fi[sn].sum(0),node_corrected_force_N=fc[sn].sum(0),CFD_abs_fraction=np.linalg.norm(ff[sf],axis=1).sum()/np.linalg.norm(ff,axis=1).sum(),node_abs_fraction=np.linalg.norm(fc[sn],axis=1).sum()/np.linalg.norm(fc,axis=1).sum()))
 regions[name]=rows
save(B/'INDEPENDENT_VERIFICATION.json',dict(checks=checks,passed=all(checks.values()),components=ind,tests='Reloaded CSR + saved low-rank correction; analytic linear-pressure triangle moment independent of exported quadrature; torque-center translation; shuffled nonconsecutive IDs; density conversion exactly once, then integrated N; reject anchors.',source_and_nodal_x_regions=regions));assert all(checks.values())
print(json.dumps(ind,indent=2))
