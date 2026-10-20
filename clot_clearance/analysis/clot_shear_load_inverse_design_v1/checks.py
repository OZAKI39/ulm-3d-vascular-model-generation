"""Only new field, surface, transfer, and external callback path checks."""
from common import *
from current_surface import CurrentSurface
from shear_field_library import ShearField
from native import native,predict
import time

def audit(tool,x,v,d,t):
    z=tool.evaluate(x,v,d,t);f=z['force'];ff=z['face_force'];n=z['normal'];active=z['active']
    force=np.linalg.norm(f.sum(0)-ff.sum(0))
    moment=np.linalg.norm(np.cross(x,f).sum(0)-np.cross(z['centroid'],ff).sum(0))
    power=abs(np.sum(f*v)-np.sum(ff*z['face_velocity']))
    tangent=np.max(abs(np.sum(z['traction']*n,axis=1)),initial=0)
    budget=np.asarray(tool.W@np.linalg.norm(ff,axis=1));ratio=np.linalg.norm(f,axis=1)/np.maximum(budget,1e-100)
    assert force<1e-17 and moment<1e-20 and power<1e-18 and tangent<1e-10
    assert ratio.max()<=1+1e-12 and z['drag_power']<=1e-20 and not f[tool.n:].any()
    lumped=np.asarray(tool.W@(z['area']*active));nz=lumped>0
    vertex=tool.field(x[tool.loaded_ids],np.tile([0.,0.,1.],(len(tool.loaded_ids),1)),t)*BETA
    nodal=f[tool.loaded_ids]/np.maximum(lumped[tool.loaded_ids,None],1e-100)
    mag=np.linalg.norm(vertex,axis=1);good=(mag>.05*tool.field.theta['nominal_shear_Pa'])&(lumped[tool.loaded_ids]>0)
    angles=np.degrees(np.arccos(np.clip(np.sum(vertex[good]*nodal[good],axis=1)/(mag[good]*np.maximum(np.linalg.norm(nodal[good],axis=1),1e-100)),-1,1)))
    return dict(force_error_N=float(force),current_moment_error_Nm=float(moment),virtual_power_error_W=float(power),
                tangent_error_Pa=float(tangent),max_force_over_positive_local_budget=float(ratio.max()),
                max_node_vs_analytic_initial_plane_direction_deg=float(angles.max(initial=0)),
                direction_scope='Initial flat zero-velocity states only; curvature and velocity affect the force law.',
                min_weight=1/3,max_weight=1/3,negative_weights=0,active_faces=int(active.sum()))

def main():
    start=time.perf_counter();ref=load_npz(O/'provenance/reference.npz');X=ref['X'];d=ref['initial_cut'].astype(float)
    results={};rng=np.random.default_rng(20261017)
    for mode in CASES:
        tool=CurrentSurface(ref,mode);z=tool.evaluate(X,np.zeros_like(X),d,HOLD)
        results[mode]=dict(flat=audit(tool,X,np.zeros_like(X),d,HOLD))
        angle=.35;rotation=np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]])
        center=X[tool.tri].mean((0,1));deformed=(X-center)@np.diag([1.1,.9,1.])@rotation.T+center
        v=rng.normal(size=X.shape)*.02;v[tool.n:]=0
        zr=tool.evaluate(deformed,v,d,HOLD)
        assert np.max(abs(zr['area']/tool.reference_area-.99))<1e-12
        assert np.max(abs(zr['normal']-np.array([0.,0.,1.])@rotation.T))<1e-12
        results[mode]['deformed']=audit(tool,deformed,v,d,HOLD)
        # Retained surfaces moving out of the drive retain passive drag.
        far=X+np.array([.002,0,0]);zf=tool.evaluate(far,v,d,HOLD)
        assert not zf['drive_force'].any() and np.linalg.norm(zf['drag_force'])>0 and zf['drag_power']<0
        results[mode]['outside_drive_passive_drag_W']=zf['drag_power']
        assert float(tool.field.envelope(0))==0 and float(tool.field.envelope(RISE))==1 and float(tool.field.envelope(OFF))==0
        # Original-ID array permutation, before transferring to native ordering.
        perm=rng.permutation(len(X));back=np.empty_like(X);back[perm]=X[perm]
        assert np.array_equal(tool.evaluate(back,np.zeros_like(X),d,HOLD)['force'],z['force'])
    uniform=CurrentSurface(ref,'uniform');v=np.tile([25.,0,0],(len(X),1));v[int(ref['clot_count']):]=0
    assert np.max(abs(uniform.evaluate(X,v,d,HOLD)['traction']))<1e-10
    assert np.max(abs(uniform.evaluate(X,np.zeros_like(X),d,HOLD)['traction']-np.array([1000.,0,0])))<1e-10
    # Independent finite-difference streamfunction and planar divergence.
    field=ShearField('counter_rotating');points=np.array([[-.2e-3,-.1e-3,375e-6],[.05e-3,.08e-3,375e-6],[.2e-3,-.05e-3,375e-6]])
    normals=np.tile([0.,0,1.],(3,1));eps=1e-8
    def psi(p):
        a=p[:,0];b=p[:,1];sigma=field.theta['width_m'];s=field.theta['separation_m']
        return sigma*np.exp(.5)/2*(np.exp(-((a+s/2)**2+b*b)/(2*sigma*sigma))-np.exp(-((a-s/2)**2+b*b)/(2*sigma*sigma)))
    gx=(psi(points+[eps,0,0])-psi(points-[eps,0,0]))/(2*eps)
    gy=(psi(points+[0,eps,0])-psi(points-[0,eps,0]))/(2*eps)
    expected=np.column_stack([-gy,gx,np.zeros(3)]);err=np.max(abs(expected-field.spatial(points,normals)))
    ux=(field(points+[eps,0,0],normals,HOLD)-field(points-[eps,0,0],normals,HOLD))/(2*eps)
    uy=(field(points+[0,eps,0],normals,HOLD)-field(points-[0,eps,0],normals,HOLD))/(2*eps)
    div=np.max(abs(ux[:,0]+uy[:,1]))/(field.U/field.theta['width_m'])
    assert err<1e-7 and div<1e-7
    # Reuse an actual prior accepted damaged state; no invented main-run cuts.
    old=load_npz(PREV/'cases/inverse/keys/final.npz');zo=uniform.evaluate(old['x'],old['velocity'],old['d'],HOLD)
    assert np.array_equal(zo['active'],(old['d'][uniform.edge_ids]<1).all(1)&zo['geometric'])
    assert not zo['face_force'][~zo['active']].any()
    ports=[];version=None
    for mode in CASES:
        ns,C=native(O/'checks'/mode,mode);s=ns.dynamics();version=ns.l.version()
        ns.l.command(f'reset_timestep 0 time {HOLD:.17g}')
        # A small material velocity makes the new slip path nontrivial.
        vel=np.zeros_like(X);vel[:ns.n]=[.02,.003,0]
        ids=ns.l.numpy.extract_atom('id')[:ns.N]-1;ns.l.numpy.extract_atom('v')[:ns.N]=vel[ids]
        ns.l.command('run 0');p=ns.state();oracle=ns.tool.evaluate(p['x'],p['velocity'],p['d'],p['time'])
        uc,fc,_=C(ref).evaluate(p['x'],p['d']);port=np.max(abs(p['force']-p['material_force']-fc-oracle['force']))
        assert port<1e-15
        h=MIN_H;state=ns.step(h)
        predicted=p['x'].copy();predicted[:ns.n]+=h*(p['velocity'][:ns.n]+h*p['force'][:ns.n]/(2*RHO*V))
        vv=p['velocity'][:ns.n]+h*(p['force'][:ns.n]+state['force'][:ns.n])/(2*RHO*V)
        assert np.max(abs(predicted-state['x']))<2e-18 and np.max(abs(vv-state['velocity'][:ns.n]))<1e-12
        ve=np.max(abs(state['force_stage_velocity']-state['velocity']));assert ve<1e-11
        w4=ns.tool.work(p,state);w8=ns.tool.work(p,state,8);assert np.max(abs(w4-w8))<1e-18 and w4[1]<=0
        assert np.array_equal(p['d'],state['d'])
        restart=O/'checks'/mode/'one_step.restart';ns.l.command(f'write_restart {restart}');ns.close()
        resumed,C=native(O/'checks'/(mode+'_restart_port'),mode,restart=restart);rs=resumed.dynamics()
        ro=resumed.tool.evaluate(rs['x'],rs['velocity'],rs['d'],rs['time'])
        re=np.max(abs(ro['force']-state['external_force']));assert re<1e-15
        assert np.array_equal(rs['x'],state['x']) and np.array_equal(rs['d'],state['d']) and np.array_equal(rs['kappa'],state['kappa'])
        resumed.close()
        ports.append(dict(mode=mode,force_port_error_N=float(port),end_velocity_consistency_m_s=float(ve),
                          restart_force_error_N=float(re),drift_work_4_vs_8_J=float(np.max(abs(w4-w8))),
                          native_check_steps=1,scope='One new slip callback step and same-state restart force reevaluation, not a repeated full restart matrix or main fracture case.'))
    dump(O/'SURFACE_TRACTION_MAPPING_AUDIT.json',dict(status='PASS',cases=results,flat_streamfunction_FD_error=float(err),flat_divergence_relative=float(div),old_damaged_state_retained_faces=int(zo['active'].sum()),
         V3_status='NOT_CHANGED_NOT_ADMITTED; no claim to repair its old affine interpolation',
         local_limit='Positive centroid weights preserve force/current moment/power but triangle quadrature and 2/3–4/3 initial nodal area variation remain. Direction differs naturally near a vortex zero; no local mesh convergence claim.'))
    dump(O/'KEY_CHECKS.json',dict(status='PASS',native_version=version,native_port=ports,new_native_check_steps=2,
         callback='Implicit linear final-half-kick slip solve, original-ID ordering, current committed edge damage',
         inherited_validation='Existing material identity, MPI kernel, contact and full restart tests reused; no kernel changes',
         library_sha256=sha(LIB),wall_s=time.perf_counter()-start))
    print('KEY_CHECKS_PASS',version,ports)

if __name__=='__main__':main()
