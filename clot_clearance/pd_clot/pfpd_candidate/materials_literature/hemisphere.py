"""The two authorized intact Ogden cases on the unchanged accepted hemisphere."""
import time
from pathlib import Path
import numpy as np
import ufl
from dolfinx import fem
from ..curved_clot_audit.native_fem import Native
from ..curved_clot_audit.common import read,dump,digest,weighted_stats
from ..intact_physics.fem_backend import build_forms
from ..fem_fracture_candidate.native import Coupled,scalar
from .ogden import batch

def run_case(r):
    start=time.perf_counter();out=Path(r['output_dir']);C=read(r['task_configuration']);old=Path(r['baseline_task']);n=Native('fine',degree=2,restore=False);reference=read(old/'requests/LM.json');reference['load']['p_Pa']=3. if r['request_id']=='MAT_LM' else 0.;reference['load']['tau_x_Pa']=3. if r['request_id']=='MAT_LM' else 0.;loadforms=build_forms(n,reference,10);load=loadforms.external
    boundary=fem.Function(n.V);bc=fem.dirichletbc(boundary,np.flatnonzero(n.fixed).astype(np.int32));fixed=np.repeat(n.fixed,3);s=Coupled(n.domain,n.V,n.u,[bc],fixed,quadrature=C['hemisphere']['material_quadrature'],external=load,prefix=r['request_id']+'_',options=dict(ksp_type='preonly',pc_type='lu',pc_factor_mat_solver_type='mumps') if 'recovery' in r else None);s.d.x.array[:]=0;s.old.x.array[:]=0
    levels=[0.] if r['request_id']=='MAT_LM' else C['hemisphere']['wall_strain_steps'];history=[];previous=0.
    if 'recovery' in r:
        saved=np.load(r['recovery']['checkpoint']);assert np.array_equal(saved['dof_coords'],n.coords) and np.array_equal(saved['dofmap'],n.V.dofmap.list);n.u.x.array[:]=saved['u_coefficients'];previous=r['recovery']['last_trusted_epsilon'];levels=[e for e in levels if e>previous]
    for step,epsilon in enumerate(levels):
        dest=out/f'steps/{step:03d}';dest.mkdir(parents=True,exist_ok=False);boundary.x.array[:]=0;boundary.x.array.reshape(-1,3)[:,0]=epsilon*n.coords[:,0];prescribed=boundary.x.array.copy();
        if 'recovery' in r and previous>0:n.u.x.array[:]*=epsilon/previous
        else:n.u.x.array.reshape(-1,3)[:,0]+=(epsilon-previous)*n.coords[:,0]
        n.u.x.array[fixed]=prescribed[fixed];iterations=[]
        def monitor(snes,it,norm):
            iterations.append(dict(iteration=it,residual_N=float(norm)));dump(dest/'NEWTON.json',iterations);print(r['request_id'],step,it,norm,flush=True)
        s.mechanics.solver.setMonitor(monitor);(dest/'EQUILIBRIUM_CALLED').write_text(str(time.time()));s.mechanics.solve();s.mechanics.solver.cancelMonitor();s.update();quality,reaction,residual=s.mechanical_quality(3*(62.5e-6)**2,3*np.pi*n.R**2,n.R,prescribed)
        gates=dict(free=quality['free_residual_relative']<=1e-8,force=quality['global_force_relative']<=1e-6,moment=quality['current_moment_relative']<=1e-6,BC=quality['BC_error_over_L']<=1e-10,J=quality['J_min']>0,finite=quality['finite'],SNES=quality['mechanical_reason']>0)
        n.checkpoint(dest);meta=read(dest/'CHECKPOINT_IDENTITY.json');meta.update(material_id=r['material_id'],method_id=r['method_id'],configuration_sha256=digest(r['task_configuration']),source_sha256=digest(__file__),schema_version='2.0',material_source_sha256=digest(Path(__file__).parent/'ogden.py'),coupled_backend_sha256=digest(Path(__file__).parent.parent/'fem_fracture_candidate/native.py'));dump(dest/'CHECKPOINT_IDENTITY.json',meta)
        row=dict(step=step,epsilon_w=epsilon,quality=quality,energies=s.energies(),PASS=all(gates.values()),failed_gates=[k for k,v in gates.items() if not v],iterations=iterations);dump(dest/'summary.json',row);history.append(row);dump(out/'HISTORY.json',history)
        if not row['PASS']:return dict(success=False,solve_status='NUMERICAL_QUALITY_FAILED',failed_gates=row['failed_gates'],response=None,diagnostic_response=row)
        previous=epsilon
    n.update_gradient();points=np.load(old/'shared/sample_points.npz');mapping=np.load(old/'shared/certified_mapping.npz');assert np.array_equal(n.original_cells,mapping['original_cell_ids']) and np.array_equal(n.vertex_ids,mapping['cell_vertex_ids']);X=points['X'];ids=mapping['sample_local_cell'];u=n.sample(X,ids);F=n.sample(X,ids,n.gf).reshape(-1,3,3);z=batch(F);sigma=z['P']@F.swapaxes(-1,-2)/z['J'][:,None,None]
    np.savez_compressed(out/'common_fields.npz',X=X,u=u,F=F,P=z['P'],sigma=sigma,J=z['J'],W0=z['W0'],W_iso=z['W_iso'],W_vol=z['W_vol'],sample_ids=points['sample_ids'],weights=points['weights'],ROI=points['ROI'])
    stats={}
    for region,mask in [('full',np.ones(len(X),bool)),('ROI',points['ROI'])]:
        stats[region]={k:weighted_stats(z[k],points['weights'],X,mask,points['sample_ids']) for k in ['W0','W_iso','W_vol','J']}
    vertices=fem.Expression(n.u,n.vertex_ref).eval(n.domain,n.entities);un=np.zeros_like(n.X);weight=np.zeros(len(n.X));np.add.at(un,n.vertex_ids.ravel(),(vertices*n.cv[:,None,None]).reshape(-1,3));np.add.at(weight,n.vertex_ids.ravel(),np.repeat(n.cv,4));un/=weight[:,None];apex=int(np.argmax(n.X[:,2]));uRMS=np.sqrt(scalar(ufl.dot(n.u,n.u)*s.dx)/n.cv.sum())
    np.savez_compressed(out/'display.npz',X=n.X,cells=n.m['cells'],triangles=n.m['triangles'],u=un,dof_coords=n.coords,u_dofs=n.u.x.array.reshape(-1,3),fixed_dofs=n.fixed,reaction_N=reaction.reshape(-1,3),external_N=s.full_extern.reshape(-1,3),prescribed_u=prescribed.reshape(-1,3))
    profiles={}
    for key in ['axis','transverse']:profiles[key+'_X']=points[key+'_X'];profiles[key+'_u']=n.sample(points[key+'_X'],mapping[key+'_local_cell'])
    np.savez_compressed(out/'profiles.npz',**profiles)
    response=dict(apex_position_m=n.X[apex],apex_displacement_m=un[apex],native_u_RMS_m=float(uRMS),common_volume=stats,energies=s.energies(),quality=quality,wall_on_clot_resultant_N=reaction.reshape(-1,3).sum(0),Psi_plus=dict(status='NOT_DEFINED_FOR_THIS_MATERIAL',value=None),Ftilde=dict(status='NOT_APPLICABLE',value=None),damage=0.,wall_strain=epsilon,hemisphere_mesh_sha256=r['geometry']['sha256'],sampling_id=reference['response_spec']['sampling_id'],measure_id=reference['response_spec']['measure_id'],wall_path_work_J=None,wall_path_work_status='not inferred from final force dot displacement; intact case path retained separately')
    artifacts={k:dict(path=str(out/file),sha256=digest(out/file)) for k,file in [('fields','common_fields.npz'),('display','display.npz'),('profiles','profiles.npz'),('native_final',str(dest.relative_to(out)/'native_checkpoint.npz'))]}
    return dict(success=True,solve_status='CONVERGED_INTACT',response=response,artifacts=artifacts,new_scientific_cases=1,load_steps=len(levels),material_calls=s.material_calls,material_seconds=s.material_seconds,interface_wall_seconds=time.perf_counter()-start,reference_scope='ONE_EXISTING_FINE_MESH_P2; imported literature material, not recalibrated; no P2/P3 or mesh convergence extension')
