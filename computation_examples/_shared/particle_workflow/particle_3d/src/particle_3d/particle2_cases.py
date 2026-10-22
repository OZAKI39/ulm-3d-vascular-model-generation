"""Permanent synthetic/real Particle-2 validation cases; no production choices."""
from dataclasses import asdict
import numpy as np
from .rbc import RBCGeometry, RBCState
from .rbc_orientation import (short_axis, quaternion_from_short_axis, rotation_matrix,
                              angular_velocity, advance_orientation, jeffery_axis_derivative)
from .rbc_integrator import advance_single_rbc, equilibrium_rbc
from .particle1_cases import AffineValidationField
from .integrator import euler_position

EPS = np.finfo(np.float64).eps
NORM_ATOL = 512*EPS
IDENTITY_RELATIVE_BUDGET = 1024*EPS
GAMMA_VALIDATION_S_INV = 20.  # SYNTHETIC_VALIDATION_PARAMETER, not a vessel shear rate
SHEAR_ANGLE_STEPS = (1/512, 1/1024, 1/2048)  # gamma*dt, VALIDATION_ONLY
PERIOD_RELATIVE_BOUND_FACTOR = 4.  # first-order budget 4*gamma*dt, fixed before runs


def initial_orientation():
    return quaternion_from_short_axis([1.,2.,3.])


def state_row(state, sample, t, dt, step, event="ACTIVE"):
    row = dict(particle_id=state.particle_id, rbc_id=state.geometry.provenance.rbc_id,
               time_s=float(t), validation_dt_s=float(dt), timestep_role="VALIDATION_ONLY", step=step,
               r=float(state.geometry.r), jeffery_lambda=float(state.geometry.jeffery_lambda),
               tetra_id=int(sample.tetra_id), inside_lumen=bool(sample.inside_lumen),
               pressure_pa=float(sample.pressure_pa), boundary_event=event,
               quaternion_norm_error=float(abs(np.linalg.norm(state.quaternion_wxyz)-1)),
               p_norm_error=float(abs(np.linalg.norm(state.short_axis)-1)),
               omega_norm_s_inv=float(np.linalg.norm(state.angular_velocity_s_inv)))
    for i, a in enumerate("xyz"):
        row.update({f"{a}_m":float(state.position_m[i]), f"V_{a}_m_s":float(state.velocity_m_s[i]),
                    f"u_{a}_m_s":float(sample.velocity_m_s[i]), f"p_{a}":float(state.short_axis[i]),
                    f"Omega_{a}_s_inv":float(state.angular_velocity_s_inv[i]),
                    f"vorticity_{a}_s_inv":float(sample.vorticity_s_inv[i])})
    row.update({f"q_{a}":float(state.quaternion_wxyz[i]) for i,a in enumerate("wxyz")})
    row.update({f"G_{i}{j}_s_inv":float(sample.velocity_gradient_s_inv[i,j]) for i in range(3) for j in range(3)})
    return row


def static_case(geometries):
    """Constant translation and zero gradient; multiple distribution geometries."""
    velocity=np.array([2e-4,-1e-4,.5e-4])
    field=AffineValidationField(velocity,np.zeros((3,3)))
    dt=.0025; count=80  # VALIDATION_ONLY
    rows=[]
    for geometry in geometries:
        q0=initial_orientation()
        state=RBCState(geometry.provenance.rbc_id,[1e-5,2e-5,3e-5],q0,velocity,[0.,0.,0.],geometry)
        x0=state.position_m.copy(); r0=rotation_matrix(q0)
        for step in range(count+1):
            row=state_row(state,field.sample(state.position_m),step*dt,dt,step)
            row["rotation_matrix_error"]=float(np.max(np.abs(rotation_matrix(state.quaternion_wxyz)-r0)))
            row["position_error_m"]=float(np.max(np.abs(state.position_m-(x0+step*dt*velocity))))
            rows.append(row)
            if step<count: state=advance_single_rbc(state,field,dt)
    return rows,dict(passed=max(r["rotation_matrix_error"] for r in rows)<=NORM_ATOL,
        max_rotation_matrix_error=max(r["rotation_matrix_error"] for r in rows),
        max_position_error_m=max(r["position_error_m"] for r in rows), validation_dt_s=dt,
        geometry_count=len(geometries),norm_atol=NORM_ATOL)


def rigid_rotation_case(geometries):
    omega=np.array([0.,0.,3.])
    g=np.array([[0.,-3.,0.],[3.,0.,0.],[0.,0.,0.]])
    field=AffineValidationField(np.zeros(3),g)
    dt=.002; count=500  # VALIDATION_ONLY
    q0=initial_orientation(); p0=short_axis(q0)
    rows=[]
    for geometry in geometries:
        state=RBCState(geometry.provenance.rbc_id,[1e-5,2e-5,3e-5],q0,[0.,0.,0.],omega,geometry)
        state=equilibrium_rbc(state,field.sample(state.position_m))
        for step in range(count+1):
            t=step*dt; angle=3*t
            exact=np.array([np.cos(angle)*p0[0]-np.sin(angle)*p0[1], np.sin(angle)*p0[0]+np.cos(angle)*p0[1],p0[2]])
            row=state_row(state,field.sample(state.position_m),t,dt,step)
            row.update({f"analytic_p_{a}":float(exact[i]) for i,a in enumerate("xyz")})
            row["axis_vector_error"]=float(np.max(np.abs(state.short_axis-exact)))
            row["omega_error_s_inv"]=float(np.max(np.abs(state.angular_velocity_s_inv-omega)))
            rows.append(row)
            if step<count: state=advance_single_rbc(state,field,dt)
    bound=64*(count+1)*EPS
    return rows,dict(passed=max(r["axis_vector_error"] for r in rows)<=bound,
                    max_axis_vector_error=max(r["axis_vector_error"] for r in rows),
                    max_omega_error_s_inv=max(r["omega_error_s_inv"] for r in rows),
                    axis_error_bound=bound,validation_dt_s=dt,geometry_count=len(geometries))


def simple_shear_case(geometries, angle_step, progress=None):
    """Independent batch of single-RBC replays, identical solver kernel.

Start p=(1,0,0) in the shear plane. Detect unwrapped angle=-pi by interpolating
the first bracketing pair. This is the shape-axis period, not the 2*pi directed
vector period. Analytical curve follows exp((W+lambda*E)t)p0, normalized.
    """
    gamma=GAMMA_VALIDATION_S_INV; dt=angle_step/gamma
    g=np.array([[0.,gamma,0.],[0.,0.,0.],[0.,0.,0.]])
    ratios=np.array([s.r for s in geometries]); lambdas=np.array([s.jeffery_lambda for s in geometries])
    analytic=np.pi/gamma*(ratios+1/ratios)
    q=np.tile(quaternion_from_short_axis([1.,0.,0.]),(len(geometries),1))
    cumulative=np.zeros(len(geometries)); previous_angle=np.zeros(len(geometries))
    measured=np.full(len(geometries),np.nan)
    max_q=np.zeros(len(geometries)); max_p=np.zeros(len(geometries)); max_identity=np.zeros(len(geometries)); max_axis=np.zeros(len(geometries))
    count=int(np.ceil(1.05*np.max(analytic)/dt))
    stride=max(1,count//600)  # only CSV plotting decimation; metrics use EVERY step
    rows=[]
    for step in range(count+1):
        t=step*dt
        p=short_axis(q)
        omega=angular_velocity(p,lambdas,g)
        identity=np.max(np.abs(np.cross(omega,p)-jeffery_axis_derivative(p,lambdas,g)),axis=1)
        max_identity=np.maximum(max_identity,identity)
        max_q=np.maximum(max_q,np.abs(np.linalg.norm(q,axis=1)-1))
        max_p=np.maximum(max_p,np.abs(np.linalg.norm(p,axis=1)-1))
        frequency=gamma*ratios/(ratios*ratios+1)
        exact=np.column_stack((np.cos(frequency*t),-np.sin(frequency*t)/ratios,np.zeros(len(ratios))))
        exact/=np.linalg.norm(exact,axis=1)[:,None]
        max_axis=np.maximum(max_axis,np.linalg.norm(p-exact,axis=1))
        if step%stride==0 or step==count:
            for j,geometry in enumerate(geometries):
                rows.append(dict(rbc_id=geometry.provenance.rbc_id,r=float(ratios[j]),jeffery_lambda=float(lambdas[j]),
                    time_s=t,validation_dt_s=dt,timestep_role="VALIDATION_ONLY",unwrapped_angle_rad=float(cumulative[j]),
                    p_x=float(p[j,0]),p_y=float(p[j,1]),p_z=float(p[j,2]),
                    analytic_p_x=float(exact[j,0]),analytic_p_y=float(exact[j,1]),analytic_p_z=0.,
                    quaternion_norm_error=float(abs(np.linalg.norm(q[j])-1)),p_norm_error=float(abs(np.linalg.norm(p[j])-1))))
        if step==count: break
        q=advance_orientation(q,omega,dt)
        new_p=short_axis(q); angle=np.arctan2(new_p[:,1],new_p[:,0])
        increment=np.arctan2(np.sin(angle-previous_angle),np.cos(angle-previous_angle))
        updated=cumulative+increment
        hit=np.isnan(measured)&(updated<=-np.pi)
        measured[hit]=(step+(-np.pi-cumulative[hit])/increment[hit])*dt
        cumulative=updated; previous_angle=angle
        if progress and step and step%10000==0: progress(step,count)
    metrics=[]
    for j,geometry in enumerate(geometries):
        prov=asdict(geometry.provenance)
        relative=abs(measured[j]-analytic[j])/analytic[j]
        metrics.append(dict(rbc_id=prov["rbc_id"],D_um=prov["D_um"],V_fL=prov["V_fL"],
            a_m=float(geometry.a_m),b_m=float(geometry.b_m),c_m=float(geometry.c_m),
            r=float(ratios[j]),jeffery_lambda=float(lambdas[j]),gamma_dot_s_inv=gamma,
            gamma_role="SYNTHETIC_VALIDATION_PARAMETER",validation_dt_s=dt,timestep_role="VALIDATION_ONLY",
            measured_period_s=float(measured[j]),analytic_period_s=float(analytic[j]),relative_error=float(relative),
            max_quaternion_norm_error=float(max_q[j]),max_p_norm_error=float(max_p[j]),
            max_jeffery_identity_error=float(max_identity[j]),max_axis_vector_error=float(max_axis[j]),
            period_relative_error_bound=PERIOD_RELATIVE_BOUND_FACTOR*angle_step,
            period_semantics="shape axis p~-p; directed vector period is 2*T_axis",
            steps=count,plot_stride=stride,
            passed=bool(np.isfinite(measured[j]) and relative<=PERIOD_RELATIVE_BOUND_FACTOR*angle_step and
                        max_q[j]<=NORM_ATOL and max_p[j]<=NORM_ATOL and max_identity[j]<=IDENTITY_RELATIVE_BUDGET*gamma)))
    return rows,metrics


def real_rbc_trajectory(field,classifier,initial,geometry,dt,progress=None):
    state=RBCState(geometry.provenance.rbc_id,initial["initial_position_m"],initial_orientation(),[0.,0.,0.],[0.,0.,0.],geometry)
    sample=field.sample(state.position_m); state=equilibrium_rbc(state,sample)
    rows=[state_row(state,sample,0.,dt,0)]
    event_record=None; segments=0; max_step=0.
    for step in range(1,int(np.ceil(initial["horizon_s"]/dt))+1):
        predicted=euler_position(state.position_m,state.velocity_m_s,dt)
        event=classifier.first_event(state.position_m,predicted); segments+=1
        max_step=max(max_step,float(np.linalg.norm(predicted-state.position_m)))
        if event is not None:
            partial_dt=event.segment_fraction*dt
            q=state.quaternion_wxyz if partial_dt==0 else advance_orientation(state.quaternion_wxyz,state.angular_velocity_s_inv,partial_dt)
            t=((step-1)+event.segment_fraction)*dt
            terminal=RBCState(state.particle_id,event.position_m,q,state.velocity_m_s,state.angular_velocity_s_inv,geometry)
            sample=field.sample(event.position_m); terminal=equilibrium_rbc(terminal,sample)
            rows.append(state_row(terminal,sample,t,dt,step,event.role))
            event_record=dict(role=event.role,time_s=float(t),position_m=event.position_m.tolist(),
                              segment_index=step-1,segment_fraction=event.segment_fraction,
                              segment_start_m=state.position_m.tolist(),unmodified_trial_endpoint_m=predicted.tolist(),
                              triangle_id=event.triangle_id,role_triangle_id=event.role_triangle_id)
            break
        state=advance_single_rbc(state,field,dt)
        sample=field.sample(state.position_m)
        rows.append(state_row(state,sample,step*dt,dt,step))
        if progress and step%500==0: progress(step)
    finite=all(np.isfinite([v for v in row.values() if isinstance(v,(int,float))]).all() for row in rows)
    qerr=max(r["quaternion_norm_error"] for r in rows); perr=max(r["p_norm_error"] for r in rows)
    verr=max(abs(row[f"V_{a}_m_s"]-row[f"u_{a}_m_s"]) for row in rows for a in "xyz")
    role=event_record["role"] if event_record else "HORIZON_REACHED"
    summary=dict(rbc_id=geometry.provenance.rbc_id,geometry=asdict(geometry),r=float(geometry.r),jeffery_lambda=float(geometry.jeffery_lambda),
        validation_dt_s=dt,timestep_role="VALIDATION_ONLY",initial_orientation_role="VALIDATION_INITIAL_ORIENTATION_ONLY",
        initial_axis=short_axis(initial_orientation()).tolist(),row_count=len(rows),checked_segment_count=segments,
        center_and_orientation_finite=finite,all_centers_inside=all(r["inside_lumen"] for r in rows),
        max_velocity_relation_error_m_s=verr,max_quaternion_norm_error=qerr,max_p_norm_error=perr,
        max_trial_step_length_m=max_step,exit_boundary=role,event=event_record,
        wall_crossing=role=="WALL",inlet_crossing=role=="INLET",
        finite_size_wall_clearance="NOT_VALIDATED_PARTICLE3",production_particle_timestep_frozen=False,
        production_orientation_distribution_frozen=False,visual_step_jump_review="PENDING_USER_REVIEW",
        passed=bool(finite and verr==0 and qerr<=NORM_ATOL and perr<=NORM_ATOL and
                    role in ["OUTLET_01","OUTLET_02","OUTLET_03"] and segments==len(rows)-1))
    return rows,summary
