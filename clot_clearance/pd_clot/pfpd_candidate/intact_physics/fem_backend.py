"""Native FEM adapter reusing the accepted mesh/checkpoint/containment helper.

Only PETSc SNES solves equilibrium. Material energy, residual and tangent share UFL.
"""
import os,time,resource,hashlib
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from .schema import RequestError,identity_hash
from .material import ufl_energy,fields
from ..curved_clot_audit.native_fem import Native,fem,ufl,PETSc,basix
from ..curved_clot_audit.common import read,dump,digest,bijective_map,weighted_stats,ROOT,OLD
from dolfinx.fem.petsc import NonlinearProblem,assemble_matrix

def phi_numpy(X,request):
    if request['load_mode']=='GLOBAL':return np.ones(X.shape[:-1])
    return np.exp(-np.sum((X-np.asarray(request['load']['center_m']))**2,axis=-1)/(2*request['load']['ell_m']**2))
def wall_values(coords,request):
    a=np.zeros_like(coords)
    if request['wall_mode']=='PRESCRIBED_INPLANE_STRETCH':a[:,0]=request['wall']['epsilon_w']*coords[:,0]
    return a
def build_forms(n,r,q):
    X=ufl.SpatialCoordinate(n.domain);N=ufl.FacetNormal(n.domain)
    phi=1. if r['load_mode']=='GLOBAL' else ufl.exp(-sum((X[i]-r['load']['center_m'][i])**2 for i in range(3))/(2*r['load']['ell_m']**2))
    dx=ufl.Measure('dx',domain=n.domain,metadata={'quadrature_degree':q});ds=ufl.Measure('ds',domain=n.domain,subdomain_data=n.mt,metadata={'quadrature_degree':q})
    vv=ufl.TestFunction(n.V);F=ufl.Identity(3)+ufl.grad(n.u);psi=ufl_energy(F);energy=psi*dx;internal=ufl.derivative(energy,n.u,vv)
    ex=ufl.as_vector((1.,0.,0.));ey=ufl.as_vector((0.,1.,0.));unitT=[-phi*N,phi*(ex-N[0]*N),phi*(ey-N[1]*N)]
    unitforms=[ufl.dot(t,vv)*ds(3) for t in unitT]
    coeff=np.array([r['load'][k] for k in ['p_Pa','tau_x_Pa','tau_y_Pa']]);constants=[fem.Constant(n.domain,PETSc.ScalarType(x)) for x in coeff]
    T=sum(c*t for c,t in zip(constants,unitT));external=ufl.dot(T,vv)*ds(3)
    residual=internal-external;tangent=ufl.derivative(residual,n.u,ufl.TrialFunction(n.V))
    return SimpleNamespace(energy=energy,internal=internal,external=external,residual=residual,tangent=tangent,units=unitforms,T=T,dx=dx,ds=ds,coefficients=coeff,constants=constants,F=F)
def surface_quadrature(m,r,q):
    ref,w=basix.make_quadrature(basix.CellType.triangle,q);bary=np.column_stack((1-ref.sum(1),ref));ids=np.flatnonzero(~m['wall_faces']);verts=m['X'][m['triangles'][ids]];X=np.einsum('qv,fvi->fqi',bary,verts);weights=2*m['areas'][ids,None]*w[None,:];N=np.broadcast_to(m['normals'][ids,None,:],X.shape)
    phi=phi_numpy(X,r);mode=np.stack([-N,np.eye(3)[0]-N[...,0,None]*N,np.eye(3)[1]-N[...,1,None]*N],axis=0)*phi[None,...,None]
    force=np.einsum('fq,mfqi->mi',weights,mode);moment=np.einsum('fq,mfqi->mi',weights,np.cross(X[None],mode));coef=np.array([r['load'][k] for k in ['p_Pa','tau_x_Pa','tau_y_Pa']]);T=np.einsum('m,mfqi->fqi',coef,mode)
    return dict(X=X,weights=weights,N=N,phi=phi,traction=T,mode_traction=mode,face_ids=ids,quadrature_id=np.arange(len(w)),unit_force=force,unit_reference_moment=moment,total_force=coef@force,total_reference_moment=coef@moment,area=float(weights.sum()),phi_area=float((weights*phi).sum()),phi_min=float(phi.min()),phi_max=float(phi.max()),phi_mean=float(np.sum(weights*phi)/weights.sum()))

def prepare(r):
    t=time.perf_counter();C=read(r['task_configuration']);degree=3 if r['method_id']=='FEM_NATIVE_P3_CHECK' else 2;n=Native('fine',degree=degree,restore=False);out=Path(r['output_dir']);shared=Path(r['shared_inputs'])
    verified=read(shared/'PREPARED_IDENTITY.json')
    if verified['mesh_sha256']!=r['geometry']['sha256'] or verified['configuration_sha256']!=digest(r['task_configuration']):raise RequestError('INPUT_PROVENANCE_MISMATCH','Preflight geometry/config identity changed')
    mapping=np.load(shared/'certified_mapping.npz');assert np.array_equal(n.original_cells,mapping['original_cell_ids']);assert np.array_equal(n.vertex_ids,mapping['cell_vertex_ids'])
    mesh=np.load(C['mesh_file']);assert np.array_equal(n.X,mesh['X']);assert np.array_equal(n.m['fixed'],mesh['fixed'])
    forms=build_forms(n,r,C['quadrature']['primary']);prescribed=wall_values(n.coords,r);wallfun=fem.Function(n.V);wallfun.x.array[:]=prescribed.ravel();wallfun.x.scatter_forward()
    bc=fem.dirichletbc(wallfun,np.flatnonzero(n.fixed).astype(np.int32));mode=np.load(shared/f'load_modes_P{degree}_{r["load_mode"]}.npz')
    assert np.array_equal(n.coords,mode['dof_coords']);assert np.array_equal(n.fixed,mode['fixed'])
    ext=(forms.coefficients@mode['unit_modes']).reshape(-1,3);actual=n.vector(forms.external).reshape(-1,3)
    if np.max(abs(actual-ext))/C['scales']['force_N']>1e-10:raise RequestError('LOAD_MODE_MISMATCH','Mode combination differs from native assembly')
    return SimpleNamespace(r=r,C=C,n=n,out=out,forms=forms,prescribed=prescribed,bc=bc,ext=actual,ext_hash=hashlib.sha256(actual.tobytes()).hexdigest(),mapping=mapping,prepare_seconds=time.perf_counter()-t,stage_times={})

def near_nullspace(problem,n):
    xyz=n.coords/n.R;b=np.zeros((len(xyz),3,6));b[:,0,0]=1;b[:,1,1]=1;b[:,2,2]=1;b[:,1,3]=-xyz[:,2];b[:,2,3]=xyz[:,1];b[:,0,4]=xyz[:,2];b[:,2,4]=-xyz[:,0];b[:,0,5]=-xyz[:,1];b[:,1,5]=xyz[:,0];b[n.fixed]=0
    Q,_=np.linalg.qr(b.reshape(-1,6));vectors=[]
    for j in range(6):v=problem.A.createVecRight();v.array[:]=Q[:,j];vectors.append(v)
    ns=PETSc.NullSpace().create(vectors=vectors,comm=problem.A.comm);problem.A.setNearNullSpace(ns);return ns,vectors
def checkpoint(p):
    n=p.n;n.checkpoint(p.out);meta=read(p.out/'CHECKPOINT_IDENTITY.json');meta.update(configuration_sha256=digest(p.r['task_configuration']),source_sha256=digest(__file__),request_identity=identity_hash(p.r),request_file_sha256=digest(p.out/'request.json'),method_id=p.r['method_id'],material_id=p.r['material_id'],wall_mode=p.r['wall_mode'],load_mode=p.r['load_mode'],helper='curved_clot_audit.Native.checkpoint; full native coefficients/mesh/tags/owners')
    dump(p.out/'CHECKPOINT_IDENTITY.json',meta)
    saved=np.load(p.out/'native_checkpoint.npz');assert np.array_equal(saved['u_coefficients'],n.u.x.array);assert np.array_equal(saved['dof_coords'],n.coords)
def solve_or_load(p):
    t=time.perf_counter();n=p.n;r=p.r;C=p.C
    if r['action']=='LOAD_SAVED_BASELINE':
        expected=OLD/'reference/fine/native_dofs.npz'
        if Path(r['source']).resolve()!=expected.resolve() or r['load_mode']!='GLOBAL' or r['wall_mode']!='FIXED_DISK' or [r['load'][k] for k in ['p_Pa','tau_x_Pa','tau_y_Pa']]!=[3.,3.,0.] or n.degree!=2:raise RequestError('STATE_IDENTITY_MISMATCH','Only original fine global M P2 is the supported native baseline')
        z=np.load(expected);ids,err=bijective_map(z['X'],n.coords,n.R*1e-11);n.u.x.array[:]=z['u'][ids].ravel();n.u.x.scatter_forward();info=dict(status='SAVED_BASELINE_LOADED',snes_reason=None,new_equilibrium_solve_calls=0,coefficient_bitwise=True,source_sha256=digest(expected),source_path=str(expected));checkpoint(p)
    else:
        if n.degree==3:
            if not Path(r['source']).is_file():raise RequestError('P2_INITIAL_STATE_MISSING','LM P3 requires the saved LM P2 checkpoint; no alternate initial case')
            expected=Path(r['shared_inputs']).parent/'cases/LM/native_checkpoint.npz'
            source_result=read(expected.parent/'result.json')
            if Path(r['source']).resolve()!=expected.resolve() or not source_result['success'] or source_result['artifacts']['native_checkpoint']['sha256']!=digest(expected):raise RequestError('P2_INITIAL_STATE_IDENTITY_MISMATCH','Require the successful, unchanged corresponding LM P2 checkpoint')
            initial=np.load(r['source']);space=fem.functionspace(n.domain,('Lagrange',2,(3,)));u2=fem.Function(space);ids,_=bijective_map(initial['dof_coords'],space.tabulate_dof_coordinates(),n.R*1e-11);u2.x.array[:]=initial['u_coefficients'].reshape(-1,3)[ids].ravel();n.u.interpolate(u2)
        else:n.u.x.array[:]=p.prescribed.ravel() # zero or affine wall-compatible field, not a solved target
        n.u.x.array.reshape(-1,3)[n.fixed]=p.prescribed[n.fixed];n.u.x.scatter_forward()
        np.savez_compressed(p.out/'initial_state.npz',u_coefficients=n.u.x.array,dof_coords=n.coords)
        problem=NonlinearProblem(p.forms.residual,n.u,bcs=[p.bc],J=p.forms.tangent,petsc_options_prefix='intact_'+r['request_id']+'_',petsc_options=C['solver']['options'])
        tm=time.perf_counter();assemble_matrix(problem.A,problem.J,bcs=[p.bc]);problem.A.assemble();p.stage_times['initial_matrix_assembly_s']=time.perf_counter()-tm
        nullspace,vectors=near_nullspace(problem,n);history=[]
        def monitor(snes,it,norm):
            row=dict(iteration=it,residual_L2_N=float(norm),wall_seconds=time.perf_counter()-t,linear_iterations=int(snes.getLinearSolveIterations()));history.append(row);dump(p.out/'RESIDUAL_HISTORY.json',history)
            np.savez_compressed(p.out/'last_iterate.npz',u_coefficients=snes.getSolution().getArray(readonly=True).copy(),dof_coordinates=n.coords)
            if time.perf_counter()-t>C['solver']['max_wall_seconds']:snes.setConvergedReason(PETSc.SNES.ConvergedReason.DIVERGED_MAX_IT)
            print(r['request_id'],'Newton',it,'residual',norm,flush=True)
        problem.solver.setMonitor(monitor);(p.out/'EQUILIBRIUM_CALLED').open('x').write(time.ctime()+'\n');tm=time.perf_counter();problem.solve();n.u.x.scatter_forward();p.stage_times['SNES_solve_s']=time.perf_counter()-tm
        info=dict(status='SOLVED',snes_reason=int(problem.solver.getConvergedReason()),new_equilibrium_solve_calls=1,iterations=int(problem.solver.getIterationNumber()),linear_iterations=int(problem.solver.getLinearSolveIterations()),source_sha256=digest(p.out/'initial_state.npz'),source_path=str(p.out/'initial_state.npz'),history=history)
        checkpoint(p) # Save raw coefficients including boundary roundoff; no post-solve overwrite.
    p.stage_times['solve_or_load_s']=time.perf_counter()-t
    return info

def quadrature_extrema(p):
    n=p.n;points,w=basix.make_quadrature(basix.CellType.tetrahedron,p.C['quadrature']['primary']);Fexp=fem.Expression(n.F,points);Xexp=fem.Expression(ufl.SpatialCoordinate(n.domain),points);extrema={};finite=True;minimum=np.inf;count=0
    for start in range(0,n.nc,p.C['resources']['quad_cell_block']):
        ids=n.entities[start:start+p.C['resources']['quad_cell_block']];F=Fexp.eval(n.domain,ids).reshape(-1,len(points),3,3);X=Xexp.eval(n.domain,ids).reshape(-1,len(points),3);v=fields(F);finite &= all(np.isfinite(a).all() for a in v.values()) and np.isfinite(F).all();minimum=min(minimum,float(v['J'].min()));count+=int(np.prod(F.shape[:2]))
        scalars=dict(F_minus_I=v['F_minus_I'],J=v['J'],psi0=v['psi0'],Psi_plus=v['Psi_plus'],P_norm_Pa=np.linalg.norm(v['P'],axis=(-2,-1)),sigma_norm_Pa=np.linalg.norm(v['sigma'],axis=(-2,-1)))
        for key,a in scalars.items():
            for sense in ['min','max']:
                flat=int(np.argmin(a) if sense=='min' else np.argmax(a));cell,q=np.unravel_index(flat,a.shape);value=float(a[cell,q]);label=key+'_'+sense
                if label not in extrema or (value<extrema[label]['value'] if sense=='min' else value>extrema[label]['value']):extrema[label]=dict(value=value,position_m=X[cell,q],original_cell_id=int(n.original_cells[ids[cell]]),quadrature_id=int(q))
    return dict(finite=bool(finite),min_detF=minimum,points=count,points_per_cell=len(points),quadrature_degree=p.C['quadrature']['primary'],extrema=extrema,measure='native reference tetrahedral quadrature points; all cells; no display sampling')

def evaluate(p,solution):
    t=time.perf_counter();n=p.n;r=p.r;C=p.C;n.update_gradient();iv=n.vector(p.forms.internal).reshape(-1,3);ext=n.vector(p.forms.external).reshape(-1,3);res=iv-ext;reaction=np.where(n.fixed[:,None],res,0.);udof=n.u.x.array.reshape(-1,3);current=n.coords+udof
    scale=C['scales'];bal=dict(free_residual_L2_N=float(np.linalg.norm(res[~n.fixed])),free_residual_relative=float(np.linalg.norm(res[~n.fixed])/max(np.linalg.norm(ext[~n.fixed]),scale['free_force_floor_N'])),force_floor_N=scale['free_force_floor_N'],wall_on_clot_resultant_N=reaction.sum(0),clot_on_wall_resultant_N=-reaction.sum(0),wall_current_moment_Nm=np.cross(current,reaction).sum(0),external_resultant_N=ext.sum(0),external_reference_moment_Nm=np.cross(n.coords,ext).sum(0),external_current_moment_Nm=np.cross(current,ext).sum(0),total_force_imbalance=float(np.linalg.norm((ext+reaction).sum(0))/scale['force_N']),total_current_moment_imbalance=float(np.linalg.norm(np.cross(current,ext+reaction).sum(0))/scale['moment_Nm']),prescribed_DOF_error_m=float(np.max(abs(udof[n.fixed]-p.prescribed[n.fixed]))),fixed_DOF_rows=int(n.fixed.sum()),known_external_on_fixed_N=ext[n.fixed].sum(0),reaction_definition='wall-on-clot R_D = grad(U)-f_ext at all constrained native vector DOFs; opposite sign for clot-on-wall; moments use native DOF positions plus native displacement')
    ext_unchanged=hashlib.sha256(ext.tobytes()).hexdigest()==p.ext_hash
    energy=float(fem.assemble_scalar(fem.form(p.forms.energy)));urms=float(np.sqrt(fem.assemble_scalar(fem.form(ufl.dot(n.u,n.u)*p.forms.dx))/n.cv.sum()))
    # Vertex display only: actual native P2/P3 displacement, no recovered gradient claim.
    ut=fem.Expression(n.u,n.vertex_ref).eval(n.domain,n.entities);un=np.zeros_like(n.X);weight=np.zeros(len(n.X));np.add.at(un,n.vertex_ids.ravel(),(ut*n.cv[:,None,None]).reshape(-1,3));np.add.at(weight,n.vertex_ids.ravel(),np.repeat(n.cv,4));un/=weight[:,None]
    apex=int(np.argmax(n.X[:,2]));assert np.linalg.norm(n.X[apex]-np.array([0,0,C['R_m']]))<1e-10*n.R
    mapping=p.mapping;common=np.load(Path(r['shared_inputs'])/'sample_points.npz');X=common['X'];ids=mapping['sample_local_cell'];canon=mapping['sample_canonical_cell'];assert np.all(ids>=0);assert np.array_equal(n.original_cells[ids],canon)
    uu=n.sample(X,ids);F=n.sample(X,ids,n.gf).reshape(-1,3,3);v=fields(F);weights=common['weights'];roi=common['ROI'];sid=common['sample_ids']
    statistics={}
    for region,mask in [('full',np.ones(len(X),bool)),('ROI',roi)]:
        statistics[region]=dict(u_RMS_m=float(np.sqrt(np.average(np.sum(uu[mask]**2,axis=1),weights=weights[mask]))))
        for key in ['psi0','Psi_plus','F_minus_I','J_minus_1']:statistics[region][key]=weighted_stats(v[key],weights,X,mask,sid)
    np.savez_compressed(p.out/'common_volume_fields.npz',X=X,sample_ids=sid,weights=weights,ROI=roi,u=uu,F=F,P=v['P'],sigma=v['sigma'],psi0=v['psi0'],Psi_plus=v['Psi_plus'],J_minus_1=v['J_minus_1'],F_minus_I=v['F_minus_I'],canonical_cell_id=canon)
    profiles={}
    for name in ['axis','transverse']:
        px=common[name+'_X'];pi=mapping[name+'_local_cell'];profiles[name+'_X']=px;profiles[name+'_u']=n.sample(px,pi);pf=n.sample(px,pi,n.gf).reshape(-1,3,3);profiles[name+'_F']=pf;profiles[name+'_Psi_plus']=fields(pf)['Psi_plus']
    np.savez_compressed(p.out/'profiles.npz',**profiles)
    np.savez_compressed(p.out/'display_and_reactions.npz',X=n.X,cells=n.m['cells'],triangles=n.m['triangles'],wall_faces=n.m['wall_faces'],u=un,volume=n.m['volume'],fixed=n.m['fixed'],dof_coords=n.coords,u_dofs=udof,fixed_dofs=n.fixed,prescribed_u_dofs=p.prescribed,reaction_N=reaction,external_N=ext,internal_N=-iv)
    qp=quadrature_extrema(p);dump(p.out/'NATIVE_QUADRATURE_CHECK.json',qp)
    quality=bool(bal['free_residual_relative']<=C['gates']['free_residual'] and bal['total_force_imbalance']<=C['gates']['global_force_and_moment'] and bal['total_current_moment_imbalance']<=C['gates']['global_force_and_moment'] and bal['prescribed_DOF_error_m']<=C['gates']['prescribed_displacement_over_R']*n.R and qp['finite'] and qp['min_detF']>0 and ext_unchanged and (solution['snes_reason'] is None or solution['snes_reason']>0))
    gate_checks=dict(free_residual=bal['free_residual_relative']<=C['gates']['free_residual'],force_balance=bal['total_force_imbalance']<=C['gates']['global_force_and_moment'],current_moment=bal['total_current_moment_imbalance']<=C['gates']['global_force_and_moment'],prescribed_DOF=bal['prescribed_DOF_error_m']<=C['gates']['prescribed_displacement_over_R']*n.R,finite_quadrature=qp['finite'],positive_J=qp['min_detF']>0,external_unchanged=ext_unchanged,SNES_converged=solution['snes_reason'] is None or solution['snes_reason']>0)
    response=dict(apex_position_m=n.X[apex],apex_node_id=apex,apex_displacement_m=un[apex],native_integral_u_RMS_m=urms,nodal_volume_u_RMS_m=float(np.sqrt(np.average(np.sum(un**2,axis=1),weights=n.m['volume']))),energy_J=energy,common_volume=statistics,balance=bal,native_quadrature=qp,external_unchanged=ext_unchanged,known_surface_external_virtual_work_at_final_u_J=float(np.sum(ext*udof)),wall_loading_path_work_J=None,wall_work_status='NOT_INFERRED_FROM_FINAL_FORCE_DOT_U; prescribed-wall path not integrated',Ftilde_status='NOT_APPLICABLE')
    p.stage_times['response_evaluation_s']=time.perf_counter()-t;dump(p.out/'STAGE_TIMES.json',dict(prepare_s=p.prepare_seconds,**p.stage_times));dump(p.out/'SOLVER_SUMMARY.json',solution)
    artifacts={name:dict(path=str(p.out/file),sha256=digest(p.out/file)) for name,file in [('native_checkpoint','native_checkpoint.npz'),('native_identity','CHECKPOINT_IDENTITY.json'),('common_fields','common_volume_fields.npz'),('profiles','profiles.npz'),('display','display_and_reactions.npz')]}
    return dict(success=quality,failed_gates=[k for k,v in gate_checks.items() if not v],solve_status=solution['status'] if quality else 'NUMERICAL_QUALITY_FAILED',new_equilibrium_solve_calls=solution['new_equilibrium_solve_calls'],data_source_identity=dict(native_checkpoint=artifacts['native_checkpoint'],mesh_sha256=r['geometry']['sha256'],sampling_sha256=C['sampling']['points_sha256'],configuration_sha256=digest(r['task_configuration']),input_source=solution['source_path'],input_source_sha256=solution['source_sha256']),response=response,artifacts=artifacts,quantity_definition=dict(u_RMS='native volume integral over entire actual fine mesh; common-volume estimates separately named',common_measure=C['sampling']['measure_id'],sample_count=len(X),original_sample_count=262144,selection='fine own-domain certified points only',Psi_plus='native Green-strain positive spectral diagnostic; finite sampling estimate, not damage',profiles='61 original positions per profile; native vectors; Frobenius L2 is not per-point error',Ftilde='NOT_APPLICABLE'),reference_scope='SAVED_GLOBAL_M_LIMITED_PRIOR_EVIDENCE' if r['action']=='LOAD_SAVED_BASELINE' else ('LM_SAME_GEOMETRY_P2_P3_CHECK_PENDING' if r['request_id'].startswith('LM') else 'EXPLORATORY_SINGLE_DISCRETIZATION_'+r['request_id']),Ftilde=dict(status='NOT_APPLICABLE',value=None))
