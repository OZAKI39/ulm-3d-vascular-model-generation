"""Y-case glue around the existing svMultiPhysics executable and audited monitors."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
import yaml

from analysis.dual_track_v1.common import *
from analysis.dual_track_v1.resources import sample
import signal
write=save
def finish(out,obj):save(out/"CFD_RESULT.json",obj)


def main(level):
    from pd_clot.lammps_fragmentation.y_bifurcation_sources import YSources
    sources=YSources();c=sources.config
    w=yaml.safe_load(Path('/workspace/y_clot/deps/y_workflow/configs/stage4/workflow.yaml').read_text())
    out=B/level;mesh=json.loads((out/'MESH_QUALITY_AUDIT.json').read_text());assert mesh['status']=='PASS'
    mesh['mesh']=str((out/'SV_MESH/mesh-complete.mesh.vtu').relative_to(ROOT));mesh['geometry_hash']=sha(B/'clot_present_lumen.vtp')
    freeze=json.loads((B/'DESIGN_FREEZE.json').read_text())
    sys.path.insert(0, w["solver_vendor"])
    from vascular_validation.postprocess import SolutionMeasurements
    from vascular_validation.validation import triangle_flux
    from flow_solver_support.steady_monitor import SteadyStopMonitor
    from flow_solver_support.transient_checks import checkpoint_audit
    from flow_solver_support.solver_checks import linear_gate, nonlinear_gate
    from flow_solver_support.flow_parser import parse_solver_log

    class YMeasurements(SolutionMeasurements):
        def measure(self, velocity, pressure):
            flows, pressures, backflow = {}, {}, {}
            for face_id, role in [(2,"INLET"),(3,"OUTLET_A"),(4,"OUTLET_B")]:
                tri = self.boundary[self.tags == face_id]
                flows[role] = triangle_flux(self.points, tri, velocity)
                xyz = self.points[tri]
                area_vectors = .5*np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
                areas = np.linalg.norm(area_vectors,axis=1)
                pressures[role] = float(np.dot(areas,pressure[tri].mean(axis=1))/areas.sum())
                fluxes = np.einsum("ij,ij->i",velocity[tri].mean(axis=1),area_vectors)
                backflow[role] = float(areas[fluxes < 0].sum()/areas.sum())
            incoming = -flows["INLET"]
            outlets = {n: flows[n] for n in ("OUTLET_A","OUTLET_B")}
            total = sum(outlets.values())
            wall = np.linalg.norm(velocity[self.wall_nodes],axis=1)
            return dict(Q_target_m3_s=self.Q,Q_in_m3_s=incoming,Q_out_total_m3_s=total,
                        outlet_flows_m3_s=outlets,signed_outward_boundary_flows_m3_s=flows,
                        epsilon_Q=abs(incoming-self.Q)/self.Q,epsilon_mass=abs(total-incoming)/self.Q,
                        outlet_fractions={n:v/total for n,v in outlets.items()} if total else None,
                        area_average_pressure_pa=pressures,negative_outward_flux_area_fraction=backflow,
                        velocity_L2=self.velocity_l2(velocity),velocity_max_m_s=float(np.linalg.norm(velocity,axis=1).max()),
                        pressure_range_pa=[float(pressure.min()),float(pressure.max())],
                        velocity_finite=bool(np.isfinite(velocity).all()),pressure_finite=bool(np.isfinite(pressure).all()),
                        wall_velocity_max_m_s=float(wall.max()),wall_noslip_pass=bool(wall.max() <= 1e-10*self.Umean))

    assert sha(w["solver_executable"]) == w["solver_sha256"]
    assert json.loads(Path(w["solver_backend_validation"]).read_text())["passed"]
    build = json.loads(Path(w["solver_build_manifest"]).read_text())
    stack = build["PETSc_build"]
    library = Path(stack["prefix"]) / "lib/libpetsc.so"
    assert sha(library) == build["PETSc_library_sha256"]
    source_mesh = (ROOT / mesh["mesh"]).parent
    assert (out / "SV_MESH").resolve()==source_mesh.resolve()
    run = out / "run"
    run.mkdir()
    flow = c["flow"]
    dt, ranks = flow["steady_relaxation_dt_s"], flow["mpi_ranks"]
    q=freeze["flow"]["Q_target_m3_s"]
    assert abs(q/(flow["inlet_mean_velocity_m_s"]*mesh["faces"]["INLET"]["area_m2"])-1)<1e-10
    tree = ET.parse(w["solver_xml_template"])
    root = tree.getroot()
    for tag, value in {"Number_of_time_steps":flow["maximum_steps"],"Time_step_size":dt,
                       "Increment_in_saving_VTK_files":flow["save_interval_steps"],
                       "Increment_in_saving_restart_files":flow["save_interval_steps"],
                       "Density":c["medium"]["density_kg_m3"],"Viscosity/Value":c["medium"]["viscosity_Pa_s"]}.items():
        root.find(".//"+tag).text = str(value)
    block = root.find("Add_mesh")
    for face in block.findall("Add_face"):
        block.remove(face)
    for name in ("WALL","INLET","OUTLET_A","OUTLET_B","CLOT"):
        face = ET.SubElement(block,"Add_face",name=name)
        ET.SubElement(face,"Face_file_path").text = "../SV_MESH/mesh-surfaces/"+name+".vtp"
    equation = root.find("Add_equation")
    for bc in equation.findall("Add_BC"):
        equation.remove(bc)
    conditions = [("INLET",dict(Type="Dir",Time_dependence="Steady",Profile="Parabolic",Impose_flux="true",Zero_out_perimeter="true",Value=-q)),
                  ("WALL",dict(Type="Dir",Time_dependence="Steady",Value=0)),
                  ("CLOT",dict(Type="Dir",Time_dependence="Steady",Value=0)),
                  *[(name,dict(Type="Neu",Time_dependence="Steady",Value=p)) for name,p in zip(("OUTLET_A","OUTLET_B"),flow["outlet_pressure_Pa"],strict=True)]]
    for name, values in conditions:
        bc = ET.SubElement(equation,"Add_BC",name=name)
        for tag,value in values.items():
            ET.SubElement(bc,tag).text = str(value)
    ET.indent(tree)
    tree.write(run / "solver.xml",encoding="utf-8",xml_declaration=True)
    shutil.copy2(w["petsc_options_template"],run / "PETSC_OPTIONS.txt")
    policy = dict(save_interval_steps=flow["save_interval_steps"],steady_last_intervals=flow["steady_consecutive_intervals"],
                  velocity_change_limit=flow["steady_velocity_relative_tolerance"],flow_change_limit=1e-6,
                  mass_limit=flow["mass_relative_tolerance"])
    write(out / "INPUTS.json",dict(policy=policy,Q_target_m3_s=q,inlet_profile="Parabolic, flux normalized",
          coordinate_unit="m",velocity_unit="m/s",time_unit="s",pressure_unit="Pa",rho_unit="kg/m3",mu_unit="Pa s",
          mesh_sha256=mesh["mesh_sha256"],solver_sha256=sha(w["solver_executable"]),
          solver_xml_sha256=sha(run / "solver.xml"),petsc_library_sha256=sha(library),
          backend_validation_sha256=sha(w["solver_backend_validation"])))
    env = {k:os.environ[k] for k in ("HOME","USER","LOGNAME","LANG") if k in os.environ}
    env.update(PATH="/usr/bin:/bin:/usr/local/cuda/bin",LC_ALL="C",OMP_NUM_THREADS="1",CUDA_VISIBLE_DEVICES="0",
               LD_LIBRARY_PATH=build["runtime_library_path"],PETSC_OPTIONS=(run / "PETSC_OPTIONS.txt").read_text().strip())
    command = [stack["candidate_wrapper"],stack["launcher"],"-n",str(ranks),stack["candidate_wrapper"],w["solver_executable"],"solver.xml"]
    from analysis.dual_track_v1.measure_b import Measurements
    measure = Measurements(source_mesh / "mesh_arrays.npz",q,flow["inlet_mean_velocity_m_s"])
    measure.wall_nodes=np.unique(measure.boundary[np.isin(measure.tags,[1,5])])
    monitor = SteadyStopMonitor(measure,policy)
    resource_records=[];hard_stop=None;last_resource=0
    seen, retries = set(), []
    started = time.perf_counter()
    result_dir = run / (str(ranks)+"-procs")
    stop = None

    def observe(done=False):
        nonlocal stop
        for path in sorted(result_dir.glob("result_*.vtu"),key=lambda p:int(p.stem.rsplit("_",1)[1])):
            step = int(path.stem.rsplit("_",1)[1])
            if step == 0 or step % policy["save_interval_steps"] or step in seen:
                continue
            cp = path.with_name("stFile_%03d.bin" % step)
            if not cp.exists() or (not done and time.time()-max(cp.stat().st_mtime,path.stat().st_mtime)<2):
                continue
            try:
                checkpoint_audit(cp,step,dt,ranks=ranks)
                u,p = measure.read(path)
            except Exception as error:
                retries.append(dict(step=step,error=repr(error)))
                break
            seen.add(step)
            state = dict(measure.measure(u,p),step=step,time_s=step*dt)
            if monitor.observe(state,u) and stop is None:
                stop = dict(step=step,reason="Five consecutive steady intervals and measured finite/mass/wall gates")
                (run / "STOP_SIM").write_text("0\n")
            write(out / "PROGRESS.json",dict(states=monitor.states,intervals=monitor.intervals,stop=stop,
                                              elapsed_s=time.perf_counter()-started,read_retries=retries))

    with (run / "solver.log").open("w") as log:
        proc = subprocess.Popen(command,cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while proc.poll() is None:
            observe()
            if time.perf_counter()-last_resource>15:
                r=sample(proc.pid);resource_records.append(r);last_resource=time.perf_counter();write(out/'RESOURCES.json',resource_records)
                if r['RSS_bytes']>freeze['budgets']['process_tree_RSS_GiB']*1024**3:hard_stop='MEMORY_BUDGET'
            if time.perf_counter()-started>freeze['budgets']['per_CFD_wall_s']:hard_stop='WALL_BUDGET'
            if hard_stop:
                os.killpg(proc.pid,signal.SIGTERM)
                try:proc.wait(timeout=20)
                except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL)
                break
            time.sleep(3)
        code = proc.wait()
    write(out/"NATIVE_EXIT_STATUS.json",dict(exit_code=code,source="Python Popen.wait of own MPI launcher"))
    observe(True)
    history = parse_solver_log((run / "solver.log").read_text(errors="replace"),dt)
    write(out / "SOLVER_HISTORY.json",history)
    failures = [] if hard_stop is None else [dict(gate="resource",error=hard_stop)]
    for name, gate in [("linear",lambda:linear_gate(dict(exit_code=code,history=history))),
                       ("nonlinear",lambda:nonlinear_gate(history,max_iterations=24))]:
        try:
            gate()
        except Exception as error:
            failures.append(dict(gate=name,error=repr(error)))
    outputs = list(result_dir.glob("result_*.vtu"))
    if not outputs:
        finish(out,dict(status="FAIL",exit_code=code,health_failures=failures,error="No native VTU output"))
        raise RuntimeError("CFD produced no solution; native log preserved")
    final = max(outputs,key=lambda p:int(p.stem.rsplit("_",1)[1]))
    u,p = measure.read(final)
    metrics = measure.measure(u,p)
    checks = dict(solver=code == 0 and not failures,steady=stop is not None,
                  mass=max(metrics["epsilon_Q"],metrics["epsilon_mass"]) <= policy["mass_limit"],
                  finite=metrics["velocity_finite"] and metrics["pressure_finite"],wall=metrics["wall_noslip_pass"],
                  forward_outlets=all(v>0 for v in metrics["outlet_flows_m3_s"].values()))
    write(out/"CFD_MASS_BALANCE.json",dict(source="actual_vtu_surface_integration",**metrics))
    finish(out,dict(status="PASS_WITH_SCOPE" if all(checks.values()) else "FAIL",checks=checks,metrics=metrics,
                    exit_code=code,health_failures=failures,stop=stop,runtime_s=time.perf_counter()-started,
                    snapshot=str(final.relative_to(ROOT)),snapshot_sha256=sha(final),velocity_array="Velocity",
                    coordinate_scale_to_m=1.,velocity_scale_to_m_s=1.,time_unit="s",mesh_sha256=mesh["mesh_sha256"],
                    geometry_hash=mesh["geometry_hash"],scope="Steady SI water flow, numerical coupling case; no CFD mesh-independence or local-physics accuracy claim"))
    if not all(checks.values()):
        raise RuntimeError("Steady CFD acceptance failed; preserve diagnostics before dependent work")


if __name__ == "__main__":
    main(sys.argv[1])
