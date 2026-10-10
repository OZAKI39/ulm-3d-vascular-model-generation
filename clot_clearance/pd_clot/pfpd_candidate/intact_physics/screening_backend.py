"""Intact-only screening using the existing native FEM/material assembly.

The reused mechanical adapter holds d=0 and never constructs a phase solver.
There is no fracture solve, surrogate, superposition, or new mesh here.
"""
import hashlib
import resource
import time
from pathlib import Path
import numpy as np
import ufl
from dolfinx import fem
from petsc4py import PETSc
from ..curved_clot_audit.native_fem import Native
from ..curved_clot_audit.common import read,dump,digest,weighted_stats
from ..materials_literature.ogden import batch
from ..fem_fracture_candidate.native import Coupled,scalar
from .fem_backend import build_forms
from .policy import risk_classification

def prepare(r):
    PETSc.Log.begin();t=time.perf_counter();n=Native('fine',degree=2,restore=False)
    profile=dict(mesh_state_loading_s=time.perf_counter()-t);t=time.perf_counter()
    C=read(r['task_configuration']);forms=build_forms(n,r,10)
    bc_fun=fem.Function(n.V);bc=fem.dirichletbc(bc_fun,np.flatnonzero(n.fixed).astype(np.int32))
    fixed=np.repeat(n.fixed,3)
    s=Coupled(n.domain,n.V,n.u,[bc],fixed,quadrature=C['material_quadrature'],
              external=forms.external,prefix='screen_'+r['request_id']+'_',options=C['solver']['options'])
    s.d.x.array[:]=0;s.old.x.array[:]=0
    profile['UFL_JIT_material_setup_s']=time.perf_counter()-t
    return dict(schema_version='3.0',request=r,n=n,s=s,forms=forms,bc=bc_fun,fixed=fixed,
                external_before=s.full_extern.copy(),profile=profile,configuration=C)

def solve_or_load(p):
    r=p['request'];n=p['n'];s=p['s'];out=Path(r['output_dir']);history=[]
    epsilon=r['wall']['epsilon_w'];p['bc'].x.array.reshape(-1,3)[:,0]=epsilon*n.coords[:,0]
    prescribed=p['bc'].x.array.copy();n.u.x.array[:]=prescribed
    def monitor(snes,it,norm):
        history.append(dict(iteration=it,residual_N=float(norm)))
        dump(out/'NEWTON.json',history);print(r['request_id'],it,norm,flush=True)
    s.mechanics.solver.setMonitor(monitor)
    stage=PETSc.Log.Stage('screen_solve');stage.push();t=time.perf_counter()
    (out/'EQUILIBRIUM_CALLED').open('x').write(str(time.time()))
    s.mechanics.solve();p['profile']['nonlinear_solve_inclusive_s']=time.perf_counter()-t
    stage.pop();s.mechanics.solver.cancelMonitor()
    events={}
    for name in ['SNESSolve','SNESFunctionEval','SNESJacobianEval','KSPSolve','PCSetUp','MatLUFactorNum']:
        events[name]=PETSc.Log.Event(name).getPerfInfo(stage.id)
    p['profile']['PETSc_solve_stage_events']=events
    p['profile']['constitutive_in_solve_s']=s.material_seconds
    p['profile']['linear_solve_inclusive_s']=events['KSPSolve']['time']
    p['profile']['linear_setup_s']=events['PCSetUp']['time']
    p['profile']['residual_jacobian_inclusive_s']=events['SNESFunctionEval']['time']+events['SNESJacobianEval']['time']
    p['profile']['timing_semantics']='PETSc event times are inclusive and nested; do not sum all event times. Nonlinear solve includes material, assembly and linear solve.'
    t=time.perf_counter();s.update()
    q,reaction,res=s.mechanical_quality(3*(62.5e-6)**2,3*np.pi*n.R**2,n.R,prescribed)
    after=n.vector(p['forms'].external);unchanged=np.array_equal(after,p['external_before'])
    energies=s.energies()
    finite=q['finite'] and all(np.isfinite(v) for v in energies.values()) and np.isfinite(s.P.x.array).all()
    gates=dict(free_residual=q['free_residual_relative']<=1e-8,global_force=q['global_force_relative']<=1e-6,
               current_moment=q['current_moment_relative']<=1e-6,prescribed_DOF=q['BC_error_over_L']<=1e-10,
               native_quadrature_J=q['J_min']>0,finite=bool(finite),SNES=q['mechanical_reason']>0,
               load_identity=unchanged,intact_d_zero=bool(np.all(s.d.x.array==0)))
    p['profile']['numerical_quality_checks_s']=time.perf_counter()-t
    t=time.perf_counter();n.checkpoint(out)
    meta=read(out/'CHECKPOINT_IDENTITY.json');meta.update(schema_version='3.0',material_id=r['material_id'],
         method_id=r['method_id'],configuration_sha256=digest(r['task_configuration']),source_sha256=digest(__file__),
         material_source_sha256=digest(Path(__file__).parent.parent/'materials_literature/ogden.py'))
    dump(out/'CHECKPOINT_IDENTITY.json',meta)
    load_id=dict(before_sha256=hashlib.sha256(p['external_before'].tobytes()).hexdigest(),
                 after_sha256=hashlib.sha256(after.tobytes()).hexdigest(),bitwise_equal=unchanged)
    np.savez_compressed(out/'native_diagnostics.npz',reaction_N=reaction,external_N=after,prescribed_u=prescribed,residual_N=res)
    p['profile']['IO_s']=time.perf_counter()-t
    solution=dict(success=all(gates.values()),quality=q,gates=gates,energies=energies,reaction=reaction,
                  prescribed=prescribed,external_identity=load_id,newton=history)
    dump(out/'NUMERICAL_GATES.json',dict(quality=q,gates=gates,external_load_identity=load_id))
    return solution

def evaluate(p,sol):
    r=p['request'];n=p['n'];s=p['s'];out=Path(r['output_dir']);profile=p['profile']
    if not sol['success']:
        dump(out/'PROFILE.json',profile)
        return dict(success=False,quality=False,solve_status='NUMERICAL_QUALITY_FAILED',response=None,
                    diagnostic_response={k:v for k,v in sol.items() if k not in ['reaction','prescribed']},
                    failed_gates=[k for k,v in sol['gates'].items() if not v],new_equilibrium_solve_calls=1)
    t=time.perf_counter();n.update_gradient()
    points=np.load(Path(r['shared_inputs'])/'sample_points.npz');mapping=np.load(Path(r['shared_inputs'])/'certified_mapping.npz')
    assert np.array_equal(n.original_cells,mapping['original_cell_ids']) and np.array_equal(n.vertex_ids,mapping['cell_vertex_ids'])
    X=points['X'];ids=mapping['sample_local_cell'];u=n.sample(X,ids);F=n.sample(X,ids,n.gf).reshape(-1,3,3)
    z=batch(F);sigma=z['P']@F.swapaxes(-1,-2)/z['J'][:,None,None]
    stretch=np.sqrt(np.linalg.eigvalsh(F.swapaxes(-1,-2)@F)[:,-1])
    source=dict(W0=z['W0'],W_iso=z['W_iso'],W_vol=z['W_vol'],J=z['J'],lambda_max=stretch)
    stats={region:{k:weighted_stats(v,points['weights'],X,mask,points['sample_ids']) for k,v in source.items()}
           for region,mask in [('full',np.ones(len(X),bool)),('ROI',points['ROI'])]}
    vertices=fem.Expression(n.u,n.vertex_ref).eval(n.domain,n.entities)
    un=np.zeros_like(n.X);weight=np.zeros(len(n.X));np.add.at(un,n.vertex_ids.ravel(),(vertices*n.cv[:,None,None]).reshape(-1,3))
    np.add.at(weight,n.vertex_ids.ravel(),np.repeat(n.cv,4));un/=weight[:,None]
    # Discontinuous element-centre values for an honest field display, no smoothing of hotspots.
    centers=n.xn.mean(1);Fc=n.sample(centers,np.arange(n.nc,dtype=np.int32),n.gf).reshape(-1,3,3)
    zc=batch(Fc);lc=np.sqrt(np.linalg.eigvalsh(Fc.swapaxes(-1,-2)@Fc)[:,-1])
    profiles={}
    for k in ['axis','transverse']:
        profiles[k+'_X']=points[k+'_X'];profiles[k+'_u']=n.sample(points[k+'_X'],mapping[k+'_local_cell'])
    apex=int(np.argmax(n.X[:,2]));uRMS=np.sqrt(scalar(ufl.dot(n.u,n.u)*s.dx)/n.cv.sum())
    R=sol['reaction'].reshape(-1,3);positions=n.coords+n.u.x.array.reshape(-1,3)
    wall_work_conjugate=float(np.sum(R[:,0]*n.coords[:,0]))
    response=dict(apex_position_m=n.X[apex],apex_displacement_m=un[apex],native_u_RMS_m=float(uRMS),
                  energy_J=sol['energies']['U_el_J'],energies={k:v for k,v in sol['energies'].items() if k in ['U_el_J','intact_W_iso_J','intact_W_vol_J']},
                  wall_on_clot_resultant_N=R.sum(0),wall_on_clot_current_moment_Nm=np.cross(positions,R).sum(0),
                  wall_stretch_generalized_reaction_J=wall_work_conjugate,
                  wall_stretch_generalized_reaction_definition='sum R_x X; conjugate to dimensionless epsilon_w; not the net resultant',
                  common_volume=stats,numerical_quality=sol['quality'],external_load_identity=sol['external_identity'],
                  RISK_PROXY_W0=dict(quantity='intact total strain-energy density per reference volume',units='Pa',
                                    baseline_case='G_P',classification=risk_classification(stats['ROI']['W0']['P99'])),
                  RISK_PROXY_LAMBDA_MAX=dict(quantity='maximum singular value of F',units='1',
                                            classification=risk_classification(stats['ROI']['lambda_max']['P99'])),
                  Psi_plus=dict(status='NOT_DEFINED_FOR_THIS_MATERIAL',value=None),
                  wall_strain=r['wall']['epsilon_w'],sampling_id=r['response_spec']['sampling_id'],
                  measure_id=r['response_spec']['measure_id'],finite_sample_tail='P95/P99/max estimates; no strict tail bound',
                  fracture='NOT_SUPPORTED_FOR_FRACTURE',patient_specific_validation='NOT_ESTABLISHED')
    profile['postprocessing_s']=time.perf_counter()-t;t=time.perf_counter()
    np.savez_compressed(out/'common_fields.npz',X=X,u=u,F=F,P=z['P'],sigma=sigma,**source,
                        sample_ids=points['sample_ids'],weights=points['weights'],ROI=points['ROI'])
    np.savez_compressed(out/'display.npz',X=n.X,cells=n.m['cells'],triangles=n.m['triangles'],u=un,
                        cell_centers=centers,cell_W0=zc['W0'],cell_lambda_max=lc,cell_ids=n.original_cells)
    np.savez_compressed(out/'profiles.npz',**profiles)
    artifacts={name:dict(path=str(out/name),sha256=digest(out/name)) for name in
               ['common_fields.npz','display.npz','profiles.npz','native_checkpoint.npz','native_diagnostics.npz','NUMERICAL_GATES.json']}
    # Verify the saved mechanical identity without a new solve.
    saved=np.load(out/'native_checkpoint.npz')
    readback=bool(np.array_equal(saved['u_coefficients'],n.u.x.array) and np.array_equal(saved['dofmap'],n.V.dofmap.list))
    assert readback
    profile['IO_s']+=time.perf_counter()-t
    profile.update(peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                   nonlinear_iterations=s.mechanics.solver.getIterationNumber(),material_calls=s.material_calls,
                   saved_readback_identity=readback,MPI_ranks=1,BLAS_OMP_threads=1)
    dump(out/'PROFILE.json',profile)
    viewer=PETSc.Viewer().createASCII(str(out/'PETSC_LOG.txt'));PETSc.Log.view(viewer);viewer.destroy()
    return dict(success=True,quality=True,solve_status='CONVERGED_INTACT_SCREENING',response=response,
                artifacts=artifacts,new_equilibrium_solve_calls=1,profile=profile,
                use_approval='NUMERICAL_IMPLEMENTATION_PASSED; APPLICATION_ACCURACY_NOT_ESTABLISHED')
