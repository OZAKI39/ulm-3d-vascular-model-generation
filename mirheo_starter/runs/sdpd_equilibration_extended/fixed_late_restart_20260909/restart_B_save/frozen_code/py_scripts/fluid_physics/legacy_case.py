"""Freeze the accepted run, not an amalgam of current and historical parameters.

Mathematical acceptance definitions reimplemented from the user's local
utils/cfd_flow/{validated_contract,steady_state}.py; no pipeline imported.
Source hashes and locations are emitted alongside the case. No root license was
found in that user-owned project; no blanket license is asserted for its code.
"""
from dataclasses import asdict
import csv
import math
from pathlib import Path
import re
import struct
import numpy as np
import yaml
from py_scripts.vessel_geometry.export import load_package
from .common import require_file, read_json, sha256_file


def close(a, b, label, rtol=1e-11, atol=0.0):
    if not math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=atol):
        raise ValueError(f"SOURCE_CONFLICT {label}: {a} != {b}")


def temporal_residual(previous, current, floor):
    """Temporal convergence, NOT deviation from prescribed pressure/velocity."""
    return abs(current-previous)/max(abs(current), floor)


def time_mean(t, values, required_duration):
    t, values = np.asarray(t, float), np.asarray(values, float)
    if len(t) < 2 or not np.isfinite(t).all() or not np.isfinite(values).all() or np.any(np.diff(t) <= 0):
        raise ValueError("invalid physical-time samples")
    if t[-1]-t[0] < required_duration*(1-1e-12):
        raise ValueError("WINDOW_INSUFFICIENT")
    return float(np.trapz(values, t)/(t[-1]-t[0]))


def flow_closure(q_in, q_out):
    return {"mass_residual": abs(q_in-sum(q_out))/max(abs(q_in), np.finfo(float).tiny),
            "significant_backflow": [v < 0 and abs(v) > 0.05*abs(q_in) for v in q_out]}


def acceptance_mapping(config, code_path):
    a = config["steady_acceptance"]
    rows = [
        ("R_mass_short", "abs(mean Qin - sum mean Qout)/max(abs(mean Qin), tiny)", a["mass_residual"], "1", "TIME_AVERAGED_CONSERVATION", "steady_state.py:79"),
        ("R_mass_long", "same formula over long window", a["mass_residual"], "1", "TIME_AVERAGED_CONSERVATION", "steady_state.py:88"),
        ("physical_volume_closure", "abs(Qin-sum Qout)/max(abs(Qin),tiny), instantaneous physical aperture velocity flux checkpoint", a["physical_volume_closure"], "1", "CONSERVATION", "physical_port_flux.py:1135,1262"),
        ("R_velocity", "abs(U(c)-U(b))/max(abs(U(c)),1e-12 m/s)", a["velocity_residual"], "1", "TEMPORAL_RESIDUAL", "steady_state.py:89"),
        ("R_pressure", "max over inlet gauge and three inlet-outlet drops: abs(p(c)-p(b))/max(abs(p(c)),1 Pa)", a["pressure_residual"], "1", "TEMPORAL_RESIDUAL_NOT_TARGET_ERROR", "steady_state.py:93"),
        ("R_inlet", "abs(mean Qin_short-Qtarget)/Qtarget", a["inlet_residual"], "1", "TARGET_DEVIATION", "steady_state.py:104"),
        ("flow_fraction_drift", "max_j(max of three fj checkpoints - min), fj = Qj/sum Qout (NOT Qj/Qin)", a["flow_fraction_drift"], "absolute fraction", "TEMPORAL_DRIFT", "physical_port_flux.py:1248; steady_state.py:105"),
        ("no_significant_mean_backflow", "fail if mean Qj < 0 and abs(mean Qj) > 0.05 abs(mean Qin); both windows", 0.05, "fraction", "SIGNED_TIME_AVERAGE", "steady_state.py:73"),
        ("Q_density_consistency", "abs(Q_rho_u_over_rho0-Q_velocity)/abs(Q_rho_u_over_rho0); denominator_kind=abs_Q_rho_u_over_rho0 in accepted physical_flux_history.json; zero denominator has no accepted example", 0.01, "1", "SCALING_DIAGNOSTIC_ONLY", "steady_state.py:111; validated_contract.py:42"),
        ("rho_range", "0.9 <= rho/rho_reference <= 1.1", [a["rho_min"], a["rho_max"]], "1", "ADAPTED_PHYSICAL_DENSITY_SANITY", "validated_contract.py:43"),
    ]
    backend = {"minimum_pdf_positive": ">0", "tau": "1", "maximum_lattice_speed": "<0.05", "controller_target": "<=1e-8", "controller_controlled_flux": "<=1e-8", "full_timestep_identity": "<=1e-8"}
    return {"status": "MAPPED", "short_window_s": a["short_window_s"], "long_window_s": a["long_window_s"],
            "mean_rule": "Trapezoidal integration in physical time; old equally spaced iterations multiplied by fixed dt.",
            "scope": "Future vessel validation on the frozen finite port polygons. Periodic material tests do not execute vessel gates.",
            "criteria": [dict(name=n, formula=f, threshold=v, unit=u, category=cat, source=str(code_path/p.split(':')[0])+(':'+p.split(':',1)[1] if ':' in p else ''), status="SOURCE_VALUE") for n,f,v,u,cat,p in rows],
            "backend_specific": {k: {"source_value":v, "status":"NOT_APPLICABLE_TO_DPD", "source":str(code_path/"validated_contract.py")} for k,v in backend.items()},
            "instantaneous_reverse_particle_velocity": "ALLOWED; never clipped or removed to pass backflow gate",
            "new_material_test_thresholds": "PROPOSED in config snapshot; not inherited or user approved"}


def freeze_case(c):
    root = Path(c["legacy_project"]).resolve()
    run = root/c["accepted_run"]
    files = {"current_config": root/c["legacy_config"], "report":run/"qc/reference_scaled_base_final.json",
             "runtime_contract":run/"qc/reference_scaled_base_runtime_contract.json", "generated_lua":run/"fresh_initial_segment.lua",
             "planes":root/c["plane_contract"], "transform":root/c["rigid_transform"], "solver_wall":root/c["solver_wall"],
             "extensions":root/c["surface_run"]/"qc/extension_geometry_qc.json", "radius":root/c["surface_run"]/"qc/radius_fidelity.json",
             "historical_extension_correction":root/c["surface_run"]/"bc/extension_pressure_correction_vmtk_boundarynormal_crossseam.csv",
             "historical_bc":root/c["surface_run"]/"bc/boundary_conditions_vmtk_boundarynormal_crossseam.json",
             "surface_config":root/c["surface_run"]/"input/cfd_surface_prepare.yaml",
             "physical_flux_history":run/"qc/reference_scaled_base_physical_flux_history.json",
             "acceptance_code":root/"utils/cfd_flow/steady_state.py", "contract_code":root/"utils/cfd_flow/validated_contract.py",
             "plane_code":root/"utils/cfd_flow/physical_port_flux.py"}
    surface_config=yaml.safe_load(require_file(files['surface_config']).read_text())
    files['original_ports']=root/surface_config['paths']['cfd_preprocess_run']/"roi/port_classification.csv"
    hashes = {str(p):sha256_file(require_file(p)) for p in files.values()}
    current = yaml.safe_load(files["current_config"].read_text())
    report = read_json(files["report"])
    rc = report["runtime_contract"]
    separate = read_json(files["runtime_contract"])
    # The standalone contract is required to be exactly the accepted report's contract.
    if separate.get("contract") != rc:
        raise ValueError("SOURCE_CONFLICT accepted runtime contracts differ")
    lua = files["generated_lua"].read_text()
    def number(key):
        m = re.search(r"^"+re.escape(key)+r"\s*=\s*([-+0-9.eE]+)\s*$", lua, re.M)
        if not m:
            raise ValueError("missing generated constant "+key)
        return float(m[1])
    pairs = {"dx":"dx_m", "dt":"dt_s", "rho0_phy":"rho0_kg_m3", "nu_phy":"nu_m2_s", "bulk_viscosity_phy":"bulk_nu_m2_s", "pressure_reference_phy":"pressure_reference_pa"}
    for a,b in pairs.items():
        close(number(a),rc[b],a)
    flow_match = re.search(r"mass_flowrate=([-+0-9.eE]+)",lua)
    if not flow_match or "kind='adaptive_flux_pressure'" not in lua:
        raise ValueError("SOURCE_CONFLICT generated inlet definition")
    close(float(flow_match[1]),rc["target_mass_flow_kg_s"],"generated inlet mass flow")
    for port, value in rc["outlet_absolute_pressure_pa"].items():
        m = re.search(r"function "+port+r"_pressure\(x,y,z,t\) return ([-+0-9.eE]+) end",lua)
        if not m:
            raise ValueError("SOURCE_CONFLICT outlet not constant "+port)
        close(float(m[1]),value,port)
        close(value-rc["pressure_reference_pa"],rc["outlet_gauge_pressure_pa"][port],port+" gauge",rtol=1e-9)
    for key, old in (("density_kg_m3","rho0_kg_m3"),("kinematic_viscosity_m2_s","nu_m2_s"),("bulk_viscosity_m2_s","bulk_nu_m2_s")):
        close(current["physics"][key],rc[old],key)
    bc = current["boundary_conditions"]
    for key in ("target_mass_flow_kg_s","target_volume_flow_m3_s"):
        close(bc[key],rc[key],key)
    if bc["outlet_gauge_pressures_pa"] != rc["outlet_gauge_pressure_pa"]:
        raise ValueError("SOURCE_CONFLICT current vs accepted pressures")
    close(rc["target_mass_flow_kg_s"],rc["rho0_kg_m3"]*rc["target_volume_flow_m3_s"],"rho Q")
    close(rc["pressure_reference_pa"],rc["rho0_kg_m3"]*(rc["dx_m"]/rc["dt_s"])**2/3,"LBM pressure offset")
    history=read_json(files['physical_flux_history'])
    if history['plane_contract_sha256']!=rc['physical_plane_contract_sha256']:
        raise ValueError('SOURCE_CONFLICT accepted flux history planes')
    # Recompute the stored diagnostics; do not trust historical PASS strings.
    for sample in history['samples']:
        qout=sum(sample['ports'][n]['Q_velocity_m3_s'] for n in rc['outlet_gauge_pressure_pa'])
        qin=sample['ports']['inlet']['Q_velocity_m3_s']
        close(abs(qin-qout)/abs(qin),sample['physical_volume_closure'],'accepted physical closure',rtol=1e-9)
        for name,p in sample['ports'].items():
            qv,qm=p['Q_velocity_m3_s'],p['Q_rho_u_over_rho0_m3_s']
            close(abs(qm-qv)/abs(qm),p['R_Q_density_consistency'],'accepted density flux denominator',rtol=1e-9)
            if name!='inlet':close(qv/qout,sample['flow_fractions'][name],'accepted fraction denominator',rtol=1e-9)
    planes = read_json(files["planes"])
    if sha256_file(files["planes"]) != rc["physical_plane_file_sha256"] or planes["contract_sha256"] != rc["physical_plane_contract_sha256"]:
        raise ValueError("SOURCE_CONFLICT physical measurement planes")
    if planes["source_geometry_sha256"] != sha256_file(files["solver_wall"]):
        raise ValueError("SOURCE_CONFLICT physical planes wall")
    package = Path(c["geometry_package"])
    if sha256_file(require_file(package/"migration_manifest.json")) != c["geometry_manifest_sha256"]:
        raise ValueError("SOURCE_CONFLICT stage-one package")
    g = load_package(package)
    hashes.update({str(p):sha256_file(p) for p in package.iterdir() if p.is_file()})
    tr = read_json(files["transform"])
    mat, inv = np.array(tr["forward_homogeneous_transform_4x4"]), np.array(tr["inverse_homogeneous_transform_4x4"])
    if not np.allclose(mat@inv, np.eye(4),atol=1e-12) or not np.allclose(mat[:3,:3].T@mat[:3,:3],np.eye(3),atol=1e-12) or not np.isclose(np.linalg.det(mat[:3,:3]),1):
        raise ValueError("SOURCE_CONFLICT not a rigid transform")
    raw = files["solver_wall"].read_bytes(); count = struct.unpack_from("<I",raw,80)[0]
    dtype = np.dtype([("n","<f4",(3,)),("p","<f4",(3,3)),("a","<u2")])
    wall = np.frombuffer(raw,dtype,count=count,offset=84)["p"]
    expected = g.points_m[g.triangles[g.entity_ids==1]]@mat[:3,:3].T+mat[:3,3]
    if wall.shape != expected.shape:
        raise ValueError("SOURCE_CONFLICT rotated wall topology")
    wall_error = float(np.max(np.abs(wall-expected)))
    if wall_error > 5e-11:
        raise ValueError("SOURCE_CONFLICT rotated wall coordinates")
    ext = {r["port_id"]:r for r in read_json(files["extensions"])["boundaries"]}
    with files['original_ports'].open(encoding='utf-8-sig') as f:original_ports={r['port_id']:r for r in csv.DictReader(f)}
    with files["historical_extension_correction"].open(encoding='utf-8-sig') as f:corrections = list(csv.DictReader(f))
    historic_bc=read_json(files['historical_bc'])
    for name,row in zip(('outlet_01','outlet_02','outlet_03'),historic_bc['outlets']):
        close(row['P_solver_boundary_pa'],rc['outlet_gauge_pressure_pa'][name],'historical corrected distal pressure')
    ports = []
    for patch in g.patches:
        if patch.kind == "wall": continue
        name = next((k for k,v in planes["ports"].items() if v["source_port_id"] == patch.port_id),None)
        if name is None: raise ValueError("PORT_IDENTITY_CONFLICT")
        pp = planes["ports"][name]
        if pp["source_boundary_origin"] != patch.boundary_origin: raise ValueError("PORT_ORIGIN_CONFLICT")
        original=original_ports[patch.port_id]
        original_reference={'center_m':[float(original[k])*1e-6 for k in ('x_um','y_um','z_um')],
                            'outward_normal':[float(original[k]) for k in ('outward_normal_x','outward_normal_y','outward_normal_z')],
                            'source':str(files['original_ports'])+' # '+patch.port_id,
                            'category':'SOURCE_VALUE_CONVERTED_um_TO_m','scope':'original centerline-based physical cut reference before extension; not current distal cap or a claim of unchanged original seam tessellation',
                            'nominal_diameter_m':float(original['diameter_um'])*1e-6,
                            'historical_planned_extension_end_m':[float(original[k])*1e-6 for k in ('extension_end_x_um','extension_end_y_um','extension_end_z_um')]}
        transformed = {}
        for label, plane in pp["planes"].items():
            row = dict(plane)
            row["origin_m"] = (inv[:3,:3]@np.array(plane["origin_m"])+inv[:3,3]).tolist()
            for key in ("unit_normal","basis_u","basis_v"):
                row[key] = (inv[:3,:3]@np.array(plane[key])).tolist()
            row["coordinate_system"] = "anatomical_fMOST_m"
            row["source_category"] = "DERIVED_RIGID_INVERSE_OF_SOURCE_VALUE"
            transformed[label] = row
        if np.dot(transformed["central"]["unit_normal"],patch.outward_normal) < .99999:
            raise ValueError("PORT_NORMAL_CONFLICT")
        ports.append({"name":name,"patch":{k:v for k,v in asdict(patch).items() if k!='face_ids'},
                      "face_count":len(patch.face_ids),"measurement_planes":transformed,
                      "measurement_source":str(files["planes"])+"#/ports/"+name,
                      "legacy_plane_metadata":{k:v for k,v in pp.items() if k!='planes'},
                      "existing_extension":ext[patch.port_id],
                      "original_physical_cut_reference":original_reference,
                      "flux_sign":"positive INTO vessel" if patch.kind=="inlet" else "positive OUT of vessel"})
    def src(value, unit, key, category="SOURCE_VALUE", original=None, converted=False, file="report"):
        return dict(value=value, unit=unit, file=str(files[file]), field=key, category=category,
                    original_value=value if original is None else original, converted=converted, assumption=False)
    rho,nu = rc["rho0_kg_m3"],rc["nu_m2_s"]
    trace = {
        "density_kg_m3":src(rho,"kg/m^3","runtime_contract.rho0_kg_m3"),
        "kinematic_viscosity_m2_s":src(nu,"m^2/s","runtime_contract.nu_m2_s"),
        "dynamic_viscosity_pa_s":src(rho*nu,"Pa s","rho0_kg_m3 * nu_m2_s","DERIVED",[rho,nu],True),
        "bulk_viscosity_m2_s":src(rc["bulk_nu_m2_s"],"m^2/s","runtime_contract.bulk_nu_m2_s"),
        "target_volume_flow_m3_s":src(rc["target_volume_flow_m3_s"],"m^3/s","runtime_contract.target_volume_flow_m3_s"),
        "target_mass_flow_kg_s":src(rc["target_mass_flow_kg_s"],"kg/s","runtime_contract.target_mass_flow_kg_s"),
        "outlet_gauge_pressures_pa":src(rc["outlet_gauge_pressure_pa"],"Pa gauge","runtime_contract.outlet_gauge_pressure_pa"),
        "pressure_numerical_offset_pa":src(rc["pressure_reference_pa"],"Pa (numerical)","runtime_contract.pressure_reference_pa"),
        "temperature_K":dict(value=c["temperature_K"],unit="K",file="user request attachment 214edccc-c8eb-4c37-9be1-8abbb5ecbb3e/pasted-text.txt",field="section 1.3",category="USER_ASSUMPTION",original_value=25,original_unit="degC",converted=True,assumption=True),
    }
    case = {"schema_version":1,"case_type":"legacy_continuum_fluid_comparison","data_origin":"user_confirmed_test_values",
            "status":"SOURCE_MATCHED","geometry_package":str(package),"geometry_manifest_sha256":c["geometry_manifest_sha256"],
            "stage1_user_acceptance":"ACCEPTED_BY_USER_THIS_REQUEST; historical reports left unchanged",
            "density_kg_m3":rho,"kinematic_viscosity_m2_s":nu,"dynamic_viscosity_pa_s":rho*nu,
            "dynamic_viscosity_formula":"mu = rho * nu", "bulk_viscosity_m2_s":rc["bulk_nu_m2_s"],
            "bulk_dynamic_viscosity_pa_s":rho*rc["bulk_nu_m2_s"],"bulk_viscosity_matching":"UNVERIFIED: ordinary DPD cannot independently prescribe shear and bulk viscosity; affects compressibility/transients",
            "temperature_K":c["temperature_K"],"temperature_source":c["temperature_source"],
            "target_volume_flow_m3_s":rc["target_volume_flow_m3_s"],"target_mass_flow_kg_s":rc["target_mass_flow_kg_s"],
            "outlet_gauge_pressures_pa":rc["outlet_gauge_pressure_pa"],"pressure_reference":{"legacy_numerical_offset_pa":rc["pressure_reference_pa"],"formula":"rho0*(1/3)*(dx/dt)^2", "legacy_named_absolute_pa":rc["outlet_absolute_pressure_pa"],"thermodynamic_absolute_pressure_pa":None,"new_gauge_reference_pa":0.0,"new_reference_category":"DESIGN_CHOICE: gauge zero mapped to measured DPD equilibrium, not atmospheric pressure"},
            "temporal_target":"constant boundary tests; statistically stationary mean flow", "inlet_profile":"NOT_SPECIFIED: accepted adaptive_flux_pressure constrains total mass flow; do not inherit older parabolic profile or older Q",
            "walls":{"fixed":True,"deforming":False,"impermeable":True,"no_slip":True,"source":"user requirement; old wall_libb_continuous_q","dpd_verified":False},
            "ports":ports,"trace":trace,"accepted_report":str(files["report"]),"accepted_iteration":report["accepted_restart"]["iteration"],
            "accepted_runtime_contract":rc,"source_matching":{"current_physics_equals_accepted":True,"generated_lua_equals_accepted":True,"plane_and_wall_hash_match":True,"rotated_wall_triangle_max_error_m":wall_error,"tolerance_m":5e-11,"geometry_modified":False,"rigid_transform_source":str(files["transform"]),"old_binary_sha256":rc["binary_sha256"],"current_config_binary_sha256":current["backend"]["musubi_expected_sha256"] if "backend" in current else None},
            "existing_extensions":"Retain existing distal caps and extensions. Do not add their length again.",
            "historical_extension_pressure_correction":{"source":str(files["historical_extension_correction"]),"rows":corrections,"scope":"Historical old Q; preserve accepted distal gauge pressures, do not recompute with current target or apply correction twice."},
            "geometry_scales":{"old_reference_dx_m":rc["dx_m"],"source_centerline_radius_min_um":read_json(files["radius"])["source_qc"]["radius_min_um"],"source_centerline_radius_scope":"109 source nodes; 1 um radius floor, not a proven minimum of the complete lumen", "measured_port_min_hydraulic_diameter_m":min(p["measurement_planes"]["central"]["local_hydraulic_diameter_m"] for p in ports),"global_min_lumen_diameter_m":None},
            "future_constraints":{"animal":"mouse","RBC_source":"same mouse","RBC_diameter_m":None,"RBC_membrane_parameters":None,"bubble_diameter_m":[2e-6,4e-6],"purpose":"transport, target arrival, adhesion; no ultrasound oscillation/collapse","simulated_this_stage":False},
            "geometry_volume_m3":g.checks["original_topology"]["signed_volume_m3"],"wall_area_m2":next(p.area_m2 for p in g.patches if p.kind=="wall")}
    # Find current binary declaration without relying on an assumed section name.
    def find_key(obj,key):
        if isinstance(obj,dict):
            if key in obj:return obj[key]
            for v in obj.values():
                got=find_key(v,key)
                if got is not None:return got
        return None
    case["source_matching"]["current_config_binary_sha256"]=find_key(current,"musubi_expected_sha256")
    case["source_matching"]["binary_note"]="Current config references another backend build; accepted generated snapshot/report are authoritative; no old solver executed."
    case['pressure_reference']['historical_gauge_origin']=historic_bc['pressure_reference']
    case['pressure_reference']['historical_gauge_origin_source']=str(files['historical_bc'])+' # pressure_reference'
    case['pressure_reference']['historical_gauge_origin_note']='Global structural leaves zero gauge in the originating 1D test model; accepted distal corrected pressures match exactly. Not measured atmospheric/absolute pressure.'
    case['source_matching']['accepted_flux_diagnostics_recomputed']=True
    case['source_matching']['historical_profile_and_Q_not_inherited']={'file':str(files['historical_bc']),'old_Q_m3_s':historic_bc['inlet']['flow_rate_m3_s'],'old_profile':historic_bc['inlet']['profile'],'accepted_inlet':'adaptive_flux_pressure with accepted total Q; no profile specified'}
    acceptance=acceptance_mapping(current,root/"utils/cfd_flow")
    acceptance['accepted_flux_history_source']=str(files['physical_flux_history'])
    acceptance['density_flux_and_fraction_denominators_verified_from_actual_samples']=True
    return case, acceptance, hashes, g
