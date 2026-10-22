"""Read-only frozen V2 decomposition; no meshing, CFD or acceptance changes."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import pyvista as pv
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'scripts'))
from poiseuille import reference
from pipe_probe_geometry import strict_cells
from particle_field_inputs import curl,shear
CASE=ROOT/'manufactured_pipe/R1_GMSH_V2/cfd_gpu'
OUT=HERE/'error_diagnosis'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def magnitude(x):return np.linalg.norm(x.reshape(len(x),-1),axis=1)
def stats(x):
 a=magnitude(x);l2=float(np.linalg.norm(a));rms=l2/np.sqrt(len(a))
 return dict(L2=l2,RMS=rms,median=float(np.median(a)),P95=float(np.percentile(a,95)),max=float(a.max()))
def decompose(actual,space,exact):
 total=actual-exact;es=space-exact;ed=actual-space;den=float(np.linalg.norm(exact));ns=float(np.linalg.norm(es));nd=float(np.linalg.norm(ed));nt=float(np.linalg.norm(total))
 cross=float(2*np.sum(es*ed)/den**2)
 point_res=float(np.max(np.abs(total-es-ed)))
 energy_res=abs(nt**2-ns**2-nd**2-2*float(np.sum(es*ed)))
 scale=max(float(np.max(np.abs(actual))),float(np.max(np.abs(space))),float(np.max(np.abs(exact))),np.finfo(float).tiny)
 assert point_res<=32*np.finfo(float).eps*scale
 assert energy_res<=128*np.finfo(float).eps*max(den**2,nt**2,ns**2,nd**2)
 return dict(total_relative_L2=nt/den,space_relative_L2=ns/den,discrete_relative_L2=nd/den,cross_term_normalized=cross,
  alignment_cosine=float(np.sum(es*ed)/(ns*nd)) if min(ns,nd)>64*np.finfo(float).eps*den else None,
  point_identity_max_abs=point_res,energy_identity_abs_residual=energy_res,exact_L2=den,exact_RMS=den/np.sqrt(len(exact)),
  total_absolute=stats(total),space_absolute=stats(es),discrete_absolute=stats(ed))
def classification(d):
 # PASS takes precedence in the overlapping user rules: a passing CFD band remains PASS.
 if d['total_relative_L2']<=.05:return 'PASS'
 if d['space_relative_L2']>.05:
  return 'P1_SPACE_LIMIT_SUPPORTED' if d['discrete_relative_L2']<=.05 else 'MIXED_ERROR'
 return 'DISCRETE_SOLUTION_LIMIT_SUPPORTED'
def load_inputs():
 f=read(ROOT/'FIELD_VALIDATION.json');ident=read(ROOT/'GPU_MANUFACTURED_RUN_IDENTITY.json')
 assert f['status']=='REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC'
 assert f['GPU_execution_identity_sha256']==sha(ROOT/'GPU_MANUFACTURED_RUN_IDENTITY.json')
 for k,h in f['source_sha256'].items():assert sha(ROOT/k)==h,k
 for k in ['final_VTU','final_native_checkpoint']:
  s=ident['steady_saved_state'];assert sha(ROOT/s[k])==s[k+'_sha256']
 with np.load(CASE/'SV_MESH/mesh_arrays.npz') as data:points=data['points_m'];tets=data['tetra']
 with np.load(CASE/'analysis/FIELD_SAMPLES.npz') as data:frozen={k:data[k] for k in data.files}
 plan=read(HERE/'GMSH_REFERENCE_MESH_V2_PLAN.json');proto=read(ROOT/'manufactured_pipe/PROTOCOL.json')
 return f,ident,points,tets,frozen,plan,proto

def main():
 OUT.mkdir(exist_ok=True)
 f,ident,points,tets,v,plan,proto=load_inputs()
 radius,length=plan['radius_m'],plan['length_m'];mean,mu=proto['mean_velocity_m_s'],proto['mu_Pa_s']
 exact=lambda x:reference(x,radius,length,mean,mu)
 # Direct exact values at EVERY node. No wall reset, CFD data or recovered gradients.
 node=exact(points)['u_m_s'];dx=points[tets[:,1:]]-points[tets[:,:1]]
 gp=np.linalg.solve(dx,node[tets[:,1:]]-node[tets[:,:1]]).transpose(0,2,1)
 grid=pv.read(CASE/'SV_MESH/mesh-complete.mesh.vtu')
 def sample(x,c=None):
  if c is None:c,_=strict_cells(grid,points,tets,x)
  assert np.all(c>=0),'No moving, extrapolation or deletion permitted'
  b=np.linalg.solve(dx[c].transpose(0,2,1),(x-points[tets[c,0]])[...,None])[...,0];w=np.c_[1-b.sum(1),b]
  return c,w,np.einsum('ni,nij->nj',w,node[tets[c]]),np.einsum('ni,nij->nj',w,v['nodal_velocity'][tets[c]])
 x=v['points_m'];c,w,up,uc=sample(x,v['cell'])
 sampling_roundoff=float(np.max(np.abs(uc-v['velocity'])))
 assert sampling_roundoff<=64*np.finfo(float).eps*np.max(np.abs(v['nodal_velocity']))
 uc=v['velocity']  # authoritative saved probe values, not recomputed CFD metrics
 ue=exact(x)['u_m_s'];ge=exact(x)['G_s_inv'];gc=v['gradient']
 assert np.array_equal(ue,v['exact_velocity'])
 n=-np.c_[x[:,:2]/np.linalg.norm(x[:,:2],axis=1)[:,None],np.zeros(len(x))]
 fields={'G':(gc,gp[c],ge),'s':(shear(gc,n),shear(gp[c],n),shear(ge,n)),'curl':(curl(gc),curl(gp[c]),curl(ge))}
 rows=[]
 offsets=points[tets[c]]-x[:,None,:];radial=-n
 # Exact quadratic Jensen deficit, separated in the local radial/tangential plane.
 vr=np.einsum('nij,nj->ni',offsets,radial)**2
 tangential=np.c_[-radial[:,1],radial[:,0],np.zeros(len(x))]
 vt=np.einsum('nij,nj->ni',offsets,tangential)**2
 er=-2*mean/radius**2*np.sum(w*vr,axis=1);et=-2*mean/radius**2*np.sum(w*vt,axis=1)
 assert np.max(np.abs(er+et-(up-ue)[:,2]))<1e-15
 for g,old in enumerate(f['groups']):
  m=v['groups']==g;d=decompose(uc[m],up[m],ue[m]);assert abs(d['total_relative_L2']-old['metrics']['u']['relative_L2'])<1e-14
  derivative={k:decompose(a[m],b[m],e[m]) for k,(a,b,e) in fields.items()}
  for k,val in derivative.items():assert abs(val['total_relative_L2']-old['metrics'][k]['relative_L2'])<1e-12
  axial={k:stats(val[m,2:3]) for k,val in [('total',uc-ue),('space',up-ue),('discrete',uc-up)]}
  transverse={k:dict(**stats(val[m,:2]),RMS_um_s=stats(val[m,:2])['RMS']*1e6,RMS_over_Umean=stats(val[m,:2])['RMS']/mean) for k,val in [('CFD',uc),('P1',up)]}
  for name in ['total_absolute','space_absolute','discrete_absolute']:
   d[name]['RMS_mm_s']=d[name]['RMS']*1e3;d[name]['RMS_over_Umean']=d[name]['RMS']/mean
   d[name]['L2_mm_s']=d[name]['L2']*1e3
   d[name]['statistics_mm_s']={key:value*1e3 for key,value in stats({'total_absolute':uc-ue,'space_absolute':up-ue,'discrete_absolute':uc-up}[name][m]).items()}
  d['exact_L2_mm_s']=d['exact_L2']*1e3;d['exact_RMS_mm_s']=d['exact_RMS']*1e3
  for name,val in [('CFD',uc),('P1',up)]:
   transverse[name]['statistics_um_s']={key:value*1e6 for key,value in stats(val[m,:2]).items()}
  den=np.linalg.norm(ue[m])
  rows.append(dict(wall_distance_m=old['wall_distance_m'],role=old['role'],samples=int(m.sum()),velocity=d,axial_absolute_m_s=axial,
   transverse_absolute=transverse,derivatives=derivative,classification=classification(d),
   quadratic_space_identity=dict(radial_relative_L2=float(np.linalg.norm(er[m])/den),tangential_relative_L2=float(np.linalg.norm(et[m])/den)),
   old_wall_zero_interpolant_max_difference_m_s=float(np.max(np.abs(up[m]-v['interpolant_velocity'][m])))))
 # Predeclared 16 rays, theta=(j+1/2)*2pi/16, central z=L/2; 3840 or fewer points.
 depths=np.unique(np.r_[np.geomspace(.1e-6,25e-6,80),np.linspace(25e-6,radius,161)])
 theta=(np.arange(16)+.5)*2*np.pi/16
 xx=np.stack([(radius-depths[:,None])*np.cos(theta),(radius-depths[:,None])*np.sin(theta),np.full((len(depths),16),length/2)],axis=-1).reshape(-1,3)
 cr,wr,ur,hr=sample(xx);ar=exact(xx)['u_m_s']
 profile={}
 for key,a in [('P1',ur),('CFD',hr)]:
  delta=(a-ar)[:,2].reshape(-1,16);ex=ar[:,2].reshape(-1,16)
  rel=np.linalg.norm(delta,axis=1)/np.linalg.norm(ex,axis=1)
  stable={}
  for threshold in [.05,.04,.03]:
   good=rel<=threshold;ids=np.flatnonzero(np.logical_and.accumulate(good[::-1])[::-1])
   sustained=[k for k in range(len(depths)-9) if good[k:k+10].all()]
   stable[str(threshold)]=dict(first_below_and_stays_to_center_m=float(depths[ids[0]]) if len(ids) else None,
    first_ten_consecutive_samples_start_m=float(depths[sustained[0]]) if sustained else None,
    first_ten_consecutive_samples_end_m=float(depths[sustained[0]+9]) if sustained else None)
  profile[key]=dict(relative_axial_L2=rel.tolist(),axial_RMS_m_s=np.sqrt(np.mean(delta**2,axis=1)).tolist(),stable_thresholds=stable)
 # Limited zero-CFD localization analysis, only discrete-limited failing bulk bands.
 localization=[]
 edges=np.stack([points[tets[:,j]]-points[tets[:,i]] for i in range(4) for j in range(i+1,4)],axis=1)
 for g,row in enumerate(rows):
  if row['role']!='BULK' or row['classification']!='DISCRETE_SOLUTION_LIMIT_SUPPORTED':continue
  m=v['groups']==g;cc=c[m];xxp=x[m];error=uc[m]-up[m];mag=magnitude(error)
  centroid=points[tets[cc]].mean(1);dc=radius-np.linalg.norm(centroid[:,:2],axis=1)
  ed=edges[cc];ed/=np.linalg.norm(ed,axis=2)[:,:,None]
  feats=dict(distance_to_inlet_m=xxp[:,2],distance_to_outlet_m=length-xxp[:,2],distance_to_nearest_cap_m=np.minimum(xxp[:,2],length-xxp[:,2]),
    cell_centroid_distance_to_shell_m=np.min(abs(dc[:,None]-np.array([25,50,75])*1e-6),axis=1),
    mean_edge_axial_alignment=np.mean(abs(ed[:,:,2]),axis=1),mean_edge_radial_alignment=np.mean(abs(np.einsum('nij,nj->ni',ed,-n[m])),axis=1))
  correlations={key:float(np.corrcoef(value,mag)[0,1]) if np.ptp(value)>64*np.finfo(float).eps*max(np.max(abs(value)),1e-12) and np.ptp(mag)>1e-16 else None for key,value in feats.items()}
  by_z=[dict(z_over_L=float(z/length),**stats(error[np.isclose(xxp[:,2],z,rtol=0,atol=1e-15)])) for z in np.unique(xxp[:,2])]
  localization.append(dict(wall_distance_m=row['wall_distance_m'],correlation_with_discrete_error_magnitude=correlations,
    notable_associations=[k for k,vv in correlations.items() if vv is not None and abs(vv)>=.5],
    by_axial_plane=by_z,axial_absolute_m_s=stats(error[:,2:3]),transverse_absolute_m_s=stats(error[:,:2]),
    scope='96 frozen points at z/L=0.4,0.5,0.6 only; association is not causation; cannot diagnose end regions directly'))
 qa=read(ROOT/'GMSH_REFERENCE_MESH_V2_RESULT.json');hwall=qa['wall_altitude_m']['median']
 wall_candidates=[dict(wall_distance_m=r['wall_distance_m'],space_relative_L2=r['velocity']['space_relative_L2'],h_wall_estimate_m=hwall*np.sqrt(.04/r['velocity']['space_relative_L2'])) for r in rows if r['role']=='PARTICLE_RELEVANT']
 hrec=min(hwall,min(q['h_wall_estimate_m'] for q in wall_candidates));hrec=np.floor(hrec*1e6)/1e6
 bulk_candidates=[dict(wall_distance_m=r['wall_distance_m'],space_relative_L2=r['velocity']['space_relative_L2'],target_estimate_m=plan['bulk_size_m']*np.sqrt(.04/r['velocity']['space_relative_L2'])) for r in rows if r['role']=='BULK' and r['classification']=='P1_SPACE_LIMIT_SUPPORTED']
 bulk_min=min([v['target_estimate_m'] for v in bulk_candidates],default=plan['bulk_size_m'])
 near=next(r for r in rows if r['role']=='PARTICLE_RELEVANT')
 radial_term=near['quadratic_space_identity']['radial_relative_L2']
 tangential_term=near['quadratic_space_identity']['tangential_relative_L2']
 chord_limit=plan['actual_wall_chord_m']*np.sqrt(max(0.,.04-radial_term)/tangential_term)
 angular=int(np.ceil((2*np.pi*radius/chord_limit)/4)*4)
 proposed_chord=2*radius*np.sin(np.pi/angular)
 proposal=dict(status='PROPOSAL_ONLY',mesh_generated=False,CFD_runs=0,design_model='h_new=h_old*sqrt(0.04/E_space); local design estimate, NOT a verified convergence rate',
  recommended_wall_altitude_m=float(hrec),wall_estimates=wall_candidates,old_wall_altitude_m=hwall,
  old_bulk_target_m=plan['bulk_size_m'],bulk_estimates=bulk_candidates,recommended_bulk_target_range_m=[float(np.ceil(.9*bulk_min*1e6/10)*10e-6),float(np.floor(bulk_min*1e6/10)*10e-6)],
  wall_altitude_alone_sufficient=False,
  tangential_design_estimate=dict(current_chord_m=plan['actual_wall_chord_m'],current_intervals=plan['circumferential_intervals'],
    allowable_chord_m=float(chord_limit),proposed_intervals_multiple_of_four=angular,proposed_chord_m=float(proposed_chord),
    proposed_sagitta_m=float(radius*(1-np.cos(np.pi/angular))),
    rule='Retain measured radial-error norm without assuming it improves; require radial_norm + tangential_norm*(chord_new/chord_old)^2 <=0.04. Conditional design model only.',
    same_cylinder_geometry=True,requires_conforming_cap_rim_if_ever_authorized=True),
  bulk_scope='Only representation-supported failing central core bands; do not globally refine the 29.8/59.7 um discrete-limited shell/transition regions',
  caveats=['Fixed-point near-wall errors need not scale as h squared when only normal altitude changes.',
   'Tangential quadratic interpolation error and finite wall facets remain; wall recommendation is necessary estimate, not a guaranteed PASS.',
   'Passing bands are not failures and do not drive global refinement. No solver, BC, cap or acceptance change authorized.'])
 put(ROOT/'PROPOSED_V3_MESH_POLICY.json',proposal)
 inputs={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'FIELD_VALIDATION.json',ROOT/'GPU_MANUFACTURED_RUN_IDENTITY.json',CASE/'analysis/FIELD_SAMPLES.npz',CASE/'SV_MESH/mesh_arrays.npz',HERE/'candidate_v2/QA_SAMPLES.npz',ROOT/'PARTICLE_CFD_INPUT_CONTRACT.json',Path(__file__)]}
 result=dict(status='MIXED_FIELD_ERROR',DIAGNOSTIC_ONLY=True,CFD_runs=0,mesh_changes=False,original_field_gate_unchanged=True,input_sha256=inputs,
  point_count=960,units=dict(velocity='m/s unless explicitly suffixed',gradient_s_curl='1/s',Umean_m_s=mean),CFD_reconstruction_max_roundoff_m_s=sampling_roundoff,exact_P1_nodal_wall_values_not_overwritten=True,classification_precedence='PASS first, then A/B/C for failing bands',
  vector_identity='e_total=e_space+e_discrete; verified with 32 eps field-scale tolerance; energy with 128 eps energy-scale tolerance',
  cross_term_units='Dimensionless squared relative norm; not percent; multiply by 10000 for percent-squared',
  alignment_zero_policy='N/A if either error norm <=64 eps times reference norm',groups=rows,
  radial_profile=dict(z_over_L=.5,ray_count=16,point_count=len(xx),theta_rad=theta.tolist(),wall_distances_m=depths.tolist(),fields=profile,
   stability_definition='Report first qualifying point that stays below threshold to center, plus explicitly limited first ten-consecutive-sample window'),
  discrete_bulk_localization=localization,localization_definition='Pearson |r|>=0.5 is a descriptive association flag only; constant features=N/A',
  proposal_sha256=sha(ROOT/'PROPOSED_V3_MESH_POLICY.json'))
 np.savez_compressed(OUT/'DECOMPOSITION_SAMPLES.npz',points_m=x,cell=c,groups=v['groups'],u_exact=ue,u_P1=up,u_CFD=uc,G_P1=gp[c],
  quadratic_radial_error_m_s=er,quadratic_tangential_error_m_s=et,radial_points_m=xx,radial_cell=cr,radial_u_exact=ar,radial_u_P1=ur,radial_u_CFD=hr)
 result['sample_sha256']=sha(OUT/'DECOMPOSITION_SAMPLES.npz')
 put(ROOT/'V2_FIELD_ERROR_DECOMPOSITION.json',result)
 print(json.dumps(dict(status=result['status'],groups=[dict(d_um=r['wall_distance_m']*1e6,classification=r['classification'],total=r['velocity']['total_relative_L2'],space=r['velocity']['space_relative_L2'],discrete=r['velocity']['discrete_relative_L2'],cos=r['velocity']['alignment_cosine']) for r in rows],proposal=proposal),indent=2))
if __name__=='__main__':main()
