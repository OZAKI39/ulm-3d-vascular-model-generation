"""Moving two-scale reference; local pressure/flow updated on each reference path."""
import time
import h5py
import numpy as np
from mevo_mb.bubbles.contact import project_contacts,constrain_pair_velocities
from mevo_mb.coupling.bifurcation import BranchClassifier
from validation.history_branch_robustness.closure import mean_slip
from validation.history_branch_robustness.periodic import slow_history_step


def run_reference(out,geometry,flow,field,closure,positions,global_ids,config,*,label,dt,history):
    radius = 1.78e-6
    radii = np.full(len(global_ids),radius)
    duration = config["particles"]["duration_s"]
    steps = int(np.ceil(duration/dt))
    increments = np.zeros((steps,len(global_ids),3))
    slip,previous = np.zeros_like(positions),np.zeros_like(positions)
    x = np.array(positions,copy=True)
    active = np.ones(len(x),bool)
    classifier = BranchClassifier(geometry,global_ids)
    frames = [x.copy()]
    flags = [active.copy()]
    velocities = [np.zeros_like(x)]
    start = time.perf_counter()
    tolerance = config["particles"]["contact_tolerance_m"]
    budget = config["particles"]["contact_max_iterations"]
    drag = 6*np.pi*config["medium"]["viscosity_Pa_s"]*radius
    maximum_projection,wall_events,contact_events,peak_re = 0.,0,0,0.
    min_pressure,max_pressure = np.inf,0.

    def project(q,r):
        return project_contacts(q,r,geometry,tolerance,budget,"ckdtree")

    def constrain(q,r,v):
        v = geometry.constrain_velocity(q,r,v)
        v = constrain_pair_velocities(q,r,v,tolerance,budget,"ckdtree")
        return geometry.constrain_velocity(q,r,v)

    for step in range(1,steps+1):
        ids = np.flatnonzero(active)
        if not len(ids):
            break
        t = (step-1)*dt
        old = x[ids].copy()
        u0 = flow.velocity(old,t)
        old_velocity = constrain(old,radii[ids],u0+slip[ids])
        predicted,_,_,_ = project(old+dt*old_velocity,radii[ids])
        pressure = field.interpolate_complex_pressure(predicted)
        gradient = field.interpolate_complex_gradient(predicted)
        min_pressure = min(min_pressure,float(abs(pressure).min()))
        max_pressure = max(max_pressure,float(abs(pressure).max()))
        force = np.zeros_like(x)
        force[ids] = drag*mean_slip(closure,label,pressure,gradient)
        candidate = slow_history_step(force,slip,previous,increments,step,dt,radius,history=history,
                                     density=config["medium"]["density_kg_m3"],viscosity=config["medium"]["viscosity_Pa_s"])
        u1 = flow.velocity(predicted,t+dt)
        end_velocity = constrain(predicted,radii[ids],u1+candidate[ids])
        proposed = old+.5*dt*(old_velocity+end_velocity)
        new,wall,contacts,_ = project(proposed,radii[ids])
        maximum_projection = max(maximum_projection,float(np.linalg.norm(new-proposed,axis=1).max()))
        wall_events += int(wall.sum())
        contact_events += len(contacts)
        candidate[ids] = end_velocity-u1
        candidate[~active] = slip[~active]
        peak_re = max(peak_re,float((2*1000*radius*np.linalg.norm(candidate[ids],axis=1)/.001).max()))
        increments[step-1] = candidate-slip
        previous,slip = slip,candidate
        x[ids] = new
        events = classifier.observe(global_ids[ids],old,new,old_velocity,end_velocity,t,t+dt)
        for event in events:
            index = int(np.flatnonzero(global_ids==event["particle_id"])[0])
            active[index] = False
            x[index] = event["position_m"]
        frames.append(x.copy())
        flags.append(active.copy())
        full_velocity = np.full_like(x,np.nan)
        full_velocity[ids] = end_velocity
        velocities.append(full_velocity)
    with h5py.File(out / (label+"_dt"+str(dt)+".h5"),"x") as h:
        h["positions_m"] = np.array(frames)
        h["active"] = np.array(flags)
        h["total_velocity_m_s"] = np.array(velocities)
        h["time_s"] = np.arange(len(frames))*dt
        h["particle_id"] = global_ids
        h.attrs["reference_only"] = True
        h.attrs["complete"] = True
        h.attrs["model"] = "Adiabatic periodic variable-radius closure plus equilibrium-radius slow-slip added mass/Abel kernel"
    return dict(label=label,dt_s=dt,history_enabled=history,branch=classifier.summary(),
                final_positions_m=x.tolist(),global_ids=global_ids.tolist(),runtime_s=time.perf_counter()-start,
                wall_projection_events=wall_events,contact_events=contact_events,maximum_projection_m=maximum_projection,
                maximum_slow_slip_Re=peak_re,pressure_range_Pa=[min_pressure,max_pressure],
                scope="Validation-only two-scale moving reference; geometric walls/contacts, no secondary; fast/slow closure and omitted other-particle contacts limit attribution")
