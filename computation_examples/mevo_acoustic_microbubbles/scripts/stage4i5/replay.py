"""Auditable short checkpoint replays. Never resumes or reruns the long gold cases."""
import argparse
from collections import defaultdict
from contextlib import contextmanager
import copy
import cProfile
import hashlib
import json
import time
import traceback

import h5py
import numpy as np

from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.coupling.bifurcation import BranchClassifier
from mevo_mb.coupling.integrators import GlobalSynchronousIntegrator
from mevo_mb.flow.fields import VTKFlowField
from scripts.stage4.geometry_certificate import install as install_certificate
from scripts.stage4i5.common import BASE, GOLD, new_run, finish, write, sha


class Timers:
    """Exclusive nested wall timers. Instrumentation overhead stays explicit."""
    def __init__(self):
        self.seconds = defaultdict(float)
        self.calls = defaultdict(int)
        self.stack = []
        self.undo = []

    @contextmanager
    def measure(self, name):
        frame = [time.perf_counter(), 0.]
        self.stack.append(frame)
        try:
            yield
        finally:
            elapsed = time.perf_counter()-frame[0]
            self.stack.pop()
            self.seconds[name] += elapsed-frame[1]
            self.calls[name] += 1
            if self.stack:
                self.stack[-1][1] += elapsed

    def wrap(self, obj, key, name):
        original = getattr(obj, key)

        def wrapped(*args, **kwargs):
            with self.measure(name):
                return original(*args, **kwargs)

        self.undo.append((obj, key, original))
        setattr(obj, key, wrapped)

    def restore(self):
        for obj, key, original in reversed(self.undo):
            setattr(obj, key, original)


class DistanceProxy:
    def __init__(self, original, timers):
        self.original, self.timers = original, timers

    def EvaluateFunction(self, point):
        with self.timers.measure("wall_closest_distance"):
            return self.original.EvaluateFunction(point)

    def EvaluateGradient(self, point, normal):
        with self.timers.measure("wall_normal_including_repeated_closest_point"):
            return self.original.EvaluateGradient(point, normal)


def load_case(case):
    config = copy.deepcopy(json.loads((GOLD/"cases"/case/"RUN_IDENTITY.json").read_text())["config"])
    def frozen(path):
        path = str(path)
        if path.startswith("/workspace/mevo_acoustic_microbubbles/"):
            path = path[len("/workspace/mevo_acoustic_microbubbles/"):]
        return str(GOLD/"inputs"/path)
    config["geometry"]["mesh_surface"] = frozen(config["geometry"]["mesh_surface"])
    config["bubble_response"]["path"] = frozen(config["bubble_response"]["path"])
    config["particles"]["pair_benchmark_path"] = str(BASE/"protocol/PAIR_BACKEND_STATE.json")
    field = ComplexAcousticField.load(frozen(config["effective_acoustic_field"]["path"]))
    f = config["flow"]
    flow = VTKFlowField(frozen(f["snapshot"]), velocity_array=f["velocity_array"],
        coordinate_scale_to_m=f["coordinate_scale_to_m"], velocity_scale_to_m_s=f["velocity_scale_to_m_s"],
        spatial_method=f["spatial_method"])
    return config, field, flow


def instrumentation(tracker, timers):
    import mevo_mb.bubbles.contact as contact
    import mevo_mb.bubbles.tracker as tracking
    import scipy.optimize
    timers.wrap(contact, "contact_candidates", "contact_candidate_search")
    timers.wrap(tracking, "project_contacts", "contact_position_projection")
    timers.wrap(tracking, "constrain_pair_velocities", "contact_velocity_assembly")
    timers.wrap(scipy.optimize, "nnls", "NNLS")
    timers.wrap(tracker.response, "harmonic_response", "response_lookup")
    timers.wrap(tracker.pair, "forces", "pair_force")
    timers.wrap(tracker.flow, "velocity", "flow_interpolation")
    timers.wrap(tracker.field, "interpolate_complex_pressure", "acoustic_interpolation")
    timers.wrap(tracker.field, "interpolate_complex_gradient", "acoustic_interpolation")
    timers.wrap(tracker.monitor, "observe", "physics_validity")
    if getattr(tracker,"_i5_wall_backend",None) is not None:
        timers.wrap(tracker._i5_wall_backend,"query","wall_combined_closest_and_normal")
    tracker.cylinder._mesh_distance = DistanceProxy(tracker.cylinder._mesh_distance, timers)
    timers.wrap(tracker.cylinder, "level_and_gradient", "wall_cache_and_batch")
    # Install the original gold certificate AFTER the distance proxy.
    certificate = install_certificate(tracker.cylinder, per_particle=True)
    timers.wrap(tracker.cylinder, "project", "wall_projection")
    timers.wrap(tracker.cylinder, "constrain_velocity", "wall_velocity_constraint")
    return certificate


def replay(segment, backend="P0", *, profile=False, max_steps=None, component_timers=True):
    if profile and not component_timers:
        raise ValueError("Minimal-instrumentation timing must not use cProfile")
    kind="profiles" if profile else "p1_vtk_contact" if component_timers else "timing_confirmation"
    out = new_run(kind, dict(segment=segment, backend=backend,
        cprofile=profile, component_timers=component_timers, max_steps=max_steps,
        protocol_sha256=sha(BASE/"protocol/PROTOCOL.json")))
    replay.last_output=out
    initialization = time.perf_counter()
    config, field, flow = load_case(segment["case"])
    tracker = GlobalSynchronousIntegrator(config, field, flow)
    timers = Timers()
    backend_info = {}
    if backend != "P0":
        from scripts.stage4i5.variants import install_variant
        backend_info = install_variant(tracker, backend)
    certificate = (instrumentation(tracker, timers) if component_timers
                   else install_certificate(tracker.cylinder, per_particle=True))
    saved = np.load(BASE/"protocol"/(segment["id"]+".npz"))
    ids = saved["particle_id"].copy()
    radii = saved["radii_m"].copy()
    active = ids[saved["active"]].copy()
    x, t = saved["positions_m"][saved["active"]].copy(), float(saved["time_s"])
    frozen = saved["positions_m"].copy()
    nominal = config["particles"]["dt_s"]
    end = segment["end_time_s"]
    classifier = BranchClassifier(tracker.cylinder, ids, config.get("classification_tolerance_m"))
    original_branch = json.loads((GOLD/"cases"/segment["case"]/"BRANCH_METRICS.json").read_text())
    for event in original_branch["crossings"]:
        if event.get("crossing_time_s", np.inf) <= t:
            classifier.assignment[event["particle_id"]] = event["branch"]
            classifier.events.append(event)
    state = tracker.evaluate(x, radii[active], t)
    initial_errors = {key: float(np.max(np.abs(value-saved[key][saved["active"]]))) for key,value in state.items()}
    initial_equal = {key: bool(np.array_equal(value, saved[key][saved["active"]])) for key,value in state.items()}
    write(out/"INITIAL_RECONSTRUCTION.json", dict(bitwise=initial_equal, max_abs_error=initial_errors,
        state_choice="Recomputed original evaluator; no hidden force replacement with checkpoint values"))
    initialization_s = time.perf_counter()-initialization
    timers.seconds.clear()
    timers.calls.clear()
    rows, position_rows, velocity_rows, force_rows, mask_rows, wall_rows = [], [], [], [], [], []
    contact_hashes, crossing_rows = [], []
    full_vel = np.full((len(ids),3), np.nan)
    full_force = np.full((len(ids),2,3), np.nan)
    begin = time.perf_counter()
    cpu_begin = time.process_time()
    last_progress = begin
    profiler = cProfile.Profile() if profile else None
    if profiler:
        profiler.enable()
    try:
        macro = int(np.floor(t/nominal+1e-10))+1
        while len(active) and t < end-1e-15 and (max_steps is None or len(rows)<max_steps):
            target = min(macro*nominal,end)
            if t >= target-1e-15:
                macro += 1
                continue
            step_start = time.perf_counter()
            with timers.measure("timestep_computation"):
                speed = np.linalg.norm(state["total_velocity_m_s"], axis=1)
                candidates = config["particles"]["max_displacement_per_radius"]*radii[active]/np.maximum(speed,1e-30)
                limiter = int(np.argmin(candidates))
                particle_dt = float(candidates[limiter])
                dt = min(target-t,particle_dt)
                reason = 0 if target-t <= particle_dt else 1
                if dt <= 1e-15:
                    raise RuntimeError("Translation step collapsed")
            with timers.measure("diagnostic_recording"):
                frozen[active] = x
                mask = np.zeros(len(ids),bool)
                mask[active] = True
                position_rows.append(frozen.copy())
                mask_rows.append(mask)
                full_vel.fill(np.nan)
                full_force.fill(np.nan)
                full_vel[active] = state["total_velocity_m_s"]
                full_force[active,0] = state["primary_force_N"]
                full_force[active,1] = state["secondary_force_N"]
                velocity_rows.append(full_vel.copy())
                force_rows.append(full_force.copy())
                limiter_id = int(active[limiter])
                regime = tracker.monitor.last["pair_regime"].copy()
            old, old_state = x.copy(), state
            x, wall, contacts, correction, iterations = tracker.step(x,radii[active],t,dt,state)
            end_time = min(t+dt,target)
            state = tracker.evaluate(x,radii[active],end_time)
            with timers.measure("branch_classification"):
                crossings = classifier.observe(active,old,x,old_state["total_velocity_m_s"],
                    state["total_velocity_m_s"],t,end_time)
            tracker.monitor.record_accepted_interval(dt)
            elapsed = time.perf_counter()-step_start
            with timers.measure("diagnostic_recording"):
                global_pairs = active[np.asarray(contacts,int).reshape(-1,2)]
                contact_hashes.append(hashlib.sha256(global_pairs.astype("<i8").tobytes()).hexdigest())
                all_wall = np.zeros(len(ids),bool)
                all_wall[active] = wall
                wall_rows.append(all_wall)
                crossing_rows.extend(crossings)
                rows.append([t,end_time,dt,particle_dt,reason,limiter_id,len(active),len(contacts),int(wall.sum()),
                    correction,iterations,elapsed,*[regime[key+"_fraction"] for key in ("far_field","transition","near_field","contact")],
                    bool(tracker.monitor.last["invalid_reason"])])
                frozen[active] = x
                leaving = {event["particle_id"] for event in crossings}
                for event in crossings:
                    frozen[event["particle_id"]] = event["position_m"]
                keep = np.array([int(pid) not in leaving for pid in active],bool)
                active, x = active[keep], x[keep]
                state = {key:value[keep] for key,value in state.items()}
                t = end_time
            if leaving and len(active):
                state = tracker.evaluate(x,radii[active],t)
            if time.perf_counter()-last_progress > 15:
                write(out/"PROGRESS.json", dict(segment=segment["id"],backend=backend,time_s=t,
                    end_time_s=end,accepted_steps=len(rows),elapsed_s=time.perf_counter()-begin))
                print(f"{segment['id']} {backend} t={t:.9f}/{end:.9f} steps={len(rows)}",flush=True)
                last_progress=time.perf_counter()
    except Exception as error:
        write(out/"EXECUTION_FAILURE.json",dict(status="FAIL",error=repr(error),
            traceback=traceback.format_exc(),last_accepted_time_s=t,accepted_steps=len(rows),
            backend=backend,segment=segment,production_eligible=False))
        np.savez_compressed(out/"last_state.npz",positions_m=x,radii_m=radii[active],particle_id=active,time_s=t)
        raise
    finally:
        if profiler:
            profiler.disable()
            profiler.dump_stats(out/"profile.pstats")
        timers.restore()
        if hasattr(tracker,"_i5_cleanup"):
            tracker._i5_cleanup()
    elapsed_s = time.perf_counter()-begin
    cpu_s = time.process_time()-cpu_begin
    frozen[active]=x
    array=np.asarray(rows)
    output_begin=time.perf_counter()
    with h5py.File(out/"replay.h5","x") as h:
        for key,value in dict(rows=array, positions_m=position_rows, active=mask_rows,
            total_velocity_m_s=velocity_rows, forces_N=force_rows, wall_hits=wall_rows,
            final_positions_m=frozen, final_active_ids=active, particle_id=ids, radii_m=radii,
            contact_sha256=np.asarray(contact_hashes,dtype="S64")).items():
            h.create_dataset(key,data=value,compression="gzip",compression_opts=1)
        h.attrs["row_columns"]=json.dumps(["start_s","end_s","dt_s","particle_dt_s","limiter_reason_0_macro_1_particle",
            "limiter_global_id","active_count","contact_count","wall_hit_count","projection_m","iterations","step_wall_s",
            "far_field_fraction","transition_fraction","near_field_fraction","contact_fraction","any_invalid_flag"])
        h.attrs["complete"] = bool(t >= end-1e-15)
    output_s = time.perf_counter()-output_begin
    write(out/"crossings.json",crossing_rows)
    write(out/"branch.json",classifier.summary())
    with h5py.File(GOLD/"cases"/segment["case"]/"particles/trajectories.h5") as h:
        reference=h["positions_m"][segment["end_frame"]]
        reference_active=h["particle_id"][:][h["active"][segment["end_frame"]]]
        gold_error=float(np.max(np.linalg.norm(frozen-reference,axis=1))) if t>=end-1e-15 else None
    times=dict(timers.seconds)
    accounted=sum(times.values())
    for attr,key in [("_i5_neighbor_cache","neighbor_cache"),("_i5_wall_cache","exact_wall_cache")]:
        value=getattr(tracker,attr,None)
        if value is not None:
            backend_info[key]={name:getattr(value,name) for name in ("calls","hits","misses","rebuilds","reuses","maximum_candidates") if hasattr(value,name)}
    wall_backend=getattr(tracker,"_i5_wall_backend",None)
    if wall_backend is not None:
        backend_info["wall_backend"]={key:getattr(wall_backend,key) for key in
            ("build_seconds","calls","points_queried","query_seconds","launch_and_synchronize_seconds","transfer_seconds","host_to_device_bytes","device_to_host_bytes","zero_copy_buffer_alias_checks") if hasattr(wall_backend,key)}
    return finish(out,dict(status="COMPLETE" if t>=end-1e-15 else "DIAGNOSTIC_PREFIX",segment=segment,
        backend=backend,backend_info=backend_info,initialization_s=initialization_s,
        initial_bitwise=initial_equal,initial_max_abs_errors=initial_errors,
        elapsed_s=elapsed_s,cpu_s=cpu_s,solver_wall_s=elapsed_s-times.get("diagnostic_recording",0),
        output_HDF_s=output_s,VTK_output_s=0.,VTK_output_scope="Original Stage4 scientific loop emits HDF/events, no per-step VTK output",
        timers_exclusive_s=times,unclassified_python_and_timer_overhead_s=elapsed_s-accounted,
        calls=dict(timers.calls),accepted_steps=len(rows),dt_min_s=float(array[:,2].min()),
        dt_mean_s=float(array[:,2].mean()),dt_max_s=float(array[:,2].max()),
        start_time_s=segment["start_time_s"],end_time_s=t,contact_pair_events=int(array[:,7].sum()),
        wall_particle_events=int(array[:,8].sum()),certificate=certificate,
        pair_backend_decisions=tracker.pair.decisions,gold_endpoint_max_position_error_m=gold_error,
        gold_endpoint_active_ids_equal=bool(np.array_equal(active,reference_active)) if t>=end-1e-15 else None,
        physics_validity=tracker.monitor.summary(),cprofile=profile,component_timers=component_timers,
        timing_scope=("Explicit nested component timers, including per-point VTK calls" if component_timers
            else "No component wrappers or per-point VTK timers; common accepted-step diagnostics, dt and classifier timers retained"),
        scope="Saved-state replay; initial and endpoint reconstruction errors explicit. Not a new full gold trajectory."))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend",default="P0")
    parser.add_argument("--profile",action="store_true")
    parser.add_argument("--segments",nargs="*")
    parser.add_argument("--max-steps",type=int)
    parser.add_argument("--continue-on-failure",action="store_true")
    parser.add_argument("--no-component-timers",action="store_true")
    args=parser.parse_args()
    protocol=json.loads((BASE/"protocol/PROTOCOL.json").read_text())
    for segment in protocol["segments"]:
        if not args.segments or segment["id"] in args.segments:
            try:
                replay(segment,args.backend,profile=args.profile,max_steps=args.max_steps,
                    component_timers=not args.no_component_timers)
            except Exception as error:
                last=getattr(replay,"last_output",None)
                if last is not None and not (last/"EXECUTION_FAILURE.json").exists():
                    write(last/"EXECUTION_FAILURE.json",dict(status="FAIL",phase="initialization_or_output",
                        backend=args.backend,segment=segment,error=repr(error),traceback=traceback.format_exc(),
                        production_eligible=False))
                if not args.continue_on_failure:
                    raise
                traceback.print_exc()


if __name__=="__main__":
    main()
