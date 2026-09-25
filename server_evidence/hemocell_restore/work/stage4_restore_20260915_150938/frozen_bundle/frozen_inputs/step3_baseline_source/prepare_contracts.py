#!/usr/bin/python3
"""Read current WSL provenance and frozen Step2 geometry; never run old solvers."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import shutil
import yaml

def write(path,obj):
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def prepare(run):
    old=Path('/home/lzy/projects/ulm_3D_vascular')
    step1=Path('/home/lzy/projects/compre_output/step1/20260912_215759/geometry_contract')
    step2=Path('/home/lzy/projects/compre_output/step2/20260912_225418')
    cfgpath=old/'configs/cfd_flow.yaml';cfg=yaml.safe_load(cfgpath.read_text())
    cfg3path=old/'configs/cfd_preprocess.yaml';cfg3=yaml.safe_load(cfg3path.read_text())
    constpath=old/'utils/cfd_flow/validated_contract.py';constants={}
    for node in ast.parse(constpath.read_text()).body:
        if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
            try:constants[node.targets[0].id]=ast.literal_eval(node.value)
            except (ValueError,TypeError):pass
    base=old/Path(cfg['paths']['accepted_base_qc']).parent.parent
    runtimepath=base/'qc/reference_scaled_base_runtime_contract.json'
    runtime=json.loads(runtimepath.read_text())['contract']
    promoted=old/'outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637'
    promotedpath=promoted/'input/production_numerical_contract.json';promotion=json.loads(promotedpath.read_text())
    s3path=old/'outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/roi/boundary_conditions.json'
    s3=json.loads(s3path.read_text())
    s4path=old/cfg['paths']['source_surface_run']/'bc/boundary_conditions_vmtk_boundarynormal_crossseam.json'
    s4=json.loads(s4path.read_text())
    boundary=json.loads((step1/'metadata/boundary_contract.json').read_text())
    geometry=json.loads((step2/'diagnostics/voxelization_summary.json').read_text())
    assert geometry['step2_status']=='PASS' and geometry['human_paraview_review']=='PASS'
    rho=float(cfg['physics']['density_kg_m3']);nu=float(cfg['physics']['kinematic_viscosity_m2_s'])
    mass=float(cfg['boundary_conditions']['target_mass_flow_kg_s']);q=float(cfg['boundary_conditions']['target_volume_flow_m3_s'])
    outlets=cfg['boundary_conditions']['outlet_gauge_pressures_pa']
    assert rho==constants['RHO0_KG_M3']==runtime['rho0_kg_m3']==promotion['rho0_kg_m3']==s3['fluid']['density_kg_m3']
    assert nu==constants['KINEMATIC_VISCOSITY_M2_S']==runtime['nu_m2_s']==promotion['kinematic_viscosity_m2_s']==s3['fluid']['kinematic_viscosity_m2_s']
    assert q==constants['TARGET_VOLUME_FLOW_M3_S']==runtime['target_volume_flow_m3_s']==promotion['target_volume_flow_m3_s']
    assert mass==constants['TARGET_MASS_FLOW_KG_S']==runtime['target_mass_flow_kg_s']==promotion['target_mass_flow_kg_s']
    assert abs(mass/rho-q)/q<1e-14
    assert outlets==constants['OUTLET_GAUGE_PRESSURES_PA']==runtime['outlet_gauge_pressure_pa']==promotion['outlet_gauge_pressures_pa']
    assert s3['pressure_reference']==s4['pressure_reference']=='GLOBAL_STRUCTURAL_LEAVES_ZERO_GAUGE'
    assert cfg3['boundary_transfer']['pressure_reference']['type']=='global_structural_leaves_zero_gauge'
    assert s3['wall']['type']==s4['wall']['type']=='NO_SLIP'
    assert abs(rho*nu-s3['fluid']['dynamic_viscosity_pa_s'])<1e-15
    evidence=[]
    def ev(item,value,unit,path,context,meaning,status):
        lines=path.read_text().splitlines();token=context.split('.')[-1]
        matching=[i for i,line in enumerate(lines,1) if token in line]
        evidence.append(dict(item=item,value=value,unit=unit,source_file=str(path),source_context=context,
            source_matching_lines=matching,meaning=meaning,status=status,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    ev('fluid.density',rho,'kg/m3',cfgpath,'physics.density_kg_m3','Current Newtonian fluid reference density','CONFIRMED_FROM_CURRENT_CONFIG')
    ev('fluid.kinematic_viscosity',nu,'m2/s',cfgpath,'physics.kinematic_viscosity_m2_s','Current physical viscosity; not old lattice viscosity','CONFIRMED_FROM_CURRENT_CONFIG')
    ev('fluid.dynamic_viscosity',rho*nu,'Pa s',s3path,'fluid.dynamic_viscosity_pa_s','mu=rho*nu; independently agrees with boundary package','DERIVED')
    ev('inlet.target_volume_flow',q,'m3/s',cfgpath,'boundary_conditions.target_volume_flow_m3_s','Positive flow into domain; current target supersedes earlier s3/s4 Q','CONFIRMED_FROM_CURRENT_CONFIG')
    ev('inlet.target_mass_flow',mass,'kg/s',runtimepath,'target_mass_flow_kg_s','Actual accepted downstream target; mass/rho agrees with Q','CONFIRMED_FROM_CURRENT_OUTPUT')
    output_outlets={}
    for b in boundary['outlets']:
        name=b['id'];s4port=next(p for p in s4['outlets'] if p['port_id']==b['source_id'])
        assert outlets[name]==s4port['P_solver_boundary_pa']
        output_outlets[name]=dict(id=name,role_status='ASSUMED',source_identity=b['source_id'],condition_type='pressure',
            pressure_pa=float(outlets[name]),pressure_type='gauge',pressure_reference='GLOBAL_STRUCTURAL_LEAVES_ZERO_GAUGE',
            outward_normal=b['normal'],center_m=b['center'],numerical_method=None)
        ev(f'{name}.pressure',float(outlets[name]),'Pa',cfgpath,f'boundary_conditions.outlet_gauge_pressures_pa.{name}',
            'Pressure prescribed at artificial distal cap; preserve signed gauge pressure, including negatives','CONFIRMED_FROM_CURRENT_CONFIG')
        ev(f'{name}.extension_pressure_provenance',float(outlets[name]),'Pa',s4path,'P_solver_boundary_pa',
            f"Original gauge {s4port['P_original_1D_pa']} minus artificial-extension correction {s4port['predicted_extension_pressure_drop_pa']}; retained by current downstream contract",'CONFIRMED_FROM_CURRENT_OUTPUT')
    ev('pressure.reference','GLOBAL_STRUCTURAL_LEAVES_ZERO_GAUGE','definition',s4path,'pressure_reference','Zero at assumed source-SWC structural leaves; physical absolute pressure not established','CONFIRMED_FROM_CURRENT_OUTPUT')
    ev('pressure.reference_meaning','GAUGE_PRESSURE_REFERENCE_NOT_PHYSIOLOGICAL_ZERO','definition',old/'utils/cfd_preprocess/pipeline.py','leaf_pressure_role','Current producer explicitly defines gauge origin, not measured physiological zero','CONFIRMED_FROM_CURRENT_CODE')
    ev('pressure.old_numerical_offset',promotion['pressure_reference_pa'],'Pa',promotedpath,'pressure_reference_role','OLD_SOLVER_NUMERICAL_METHOD; must not carry this LBM offset into HemoCell physical contract','CONFIRMED_FROM_CURRENT_OUTPUT')
    ev('flow.sign','outward Q=inlet negative; outlets positive; role inflow=-outward Q','convention',step1/'metadata/boundary_contract.json','normal_sign_convention','Step1 geometric outward normals; no absolute-value masking of backflow','CONFIRMED_FROM_CURRENT_OUTPUT')
    ev('wall.condition','no_slip','condition',s3path,'wall','Physical stationary rigid no-slip wall; wall_libb is old numerical realization','CONFIRMED_FROM_CURRENT_OUTPUT')
    contract=dict(status='PASS',scope='Current traceable physical BC candidate for frozen geometry; not experimental validation',
        fluid=dict(model='Newtonian',density_kg_m3=rho,kinematic_viscosity_m2_s=nu,dynamic_viscosity_pa_s=rho*nu),
        inlet=dict(id='inlet',role_status='ASSUMED',condition_type='target_volumetric_flow',target_volume_flow_m3_s=q,target_mass_flow_kg_s=mass,
            normal_convention='geometric outward; Qin_positive_into_domain=-integral(u dot n_out)dA',numerical_method=None,
            source_identity=boundary['inlets'][0]['source_id'],profile_status='NOT_YET_SELECTED; historical PARABOLIC is not a measured profile'),
        outlets=output_outlets,wall=dict(condition='no_slip',motion='stationary',rigidity='rigid',numerical_method=None),
        pressure_reference=dict(type='gauge',zero_definition='GLOBAL_STRUCTURAL_LEAVES_ZERO_GAUGE',zero_gauge_pa=0.0,
            absolute_physiological_reference_pa=None,absolute_physiological_reference_status='UNVERIFIED',
            old_lbm_offset_role='OLD_SOLVER_NUMERICAL_METHOD; not physical blood pressure'),
        evidence=evidence,old_solver_numerical_method=dict(inlet=cfg['boundary_conditions']['inlet_boundary'],
            outlets=cfg['boundary_conditions']['outlet_boundary'],wall=cfg['boundary_conditions']['wall_boundary'],
            dt_s=runtime['dt_s'],tau=runtime['tau'],pressure_offset_pa=runtime['pressure_reference_pa'],
            transfer_to_hemocell='NO; separate API and lattice-unit derivation required'),
        historical_conflicts=dict(s3_s4_inlet_q_m3_s=s3['inlet']['flow_rate_m3_s'],current_inlet_q_m3_s=q,
            resolution='Use currently validated config/constants and accepted+promoted output target, not historical s3/s4 inlet number',
            s3_s4_profile='PARABOLIC; historical numerical profile, not physiological evidence',
            current_config_execution_mode=cfg['execution']['mode'],inspected_promoted_output_mode='VALIDATED_BASE_PROMOTION_REPLAY; previously executed, not rerun'),
        gate_checks=dict(density=True,viscosity=True,target_Q=True,three_pressures=True,pressure_reference=True,units=True,
            config_constants_runtime_promotion_agree=True,port_identities_match_step1=True,mass_to_volume_flow_agree=True),
        unverified_items=['physiological inlet/outlet identity','absolute physiological pressure','experimental viscosity/flow calibration'])
    write(run/'contracts/physical_bc_contract.json',contract)
    used=['diagnostics/closed_flag_matrix.u8','diagnostics/opened_flag_matrix.u8','diagnostics/port_label_field.u8',
        'diagnostics/cap_triangle_ids.npz','diagnostics/cap_identification.json','diagnostics/transform_check.json',
        'diagnostics/palabos_geometry.json','diagnostics/voxelization_summary.json','diagnostics/human_paraview_review.json',
        'ports.tsv','run_settings.json']
    reuse=dict(status='PASS',method='Read frozen Step2 arrays/contracts directly; no voxelization or port-identification rerun',
        input_stl=geometry['input_stl'],input_stl_sha256=geometry['input_sha256'],lattice_shape=geometry['lattice_shape'],
        dx_status='CANDIDATE_BASELINE',dx_requested_m=geometry['requested_dx_m'],dx_effective_m=geometry['effective_dx_m'],
        ref_dir=geometry['ref_dir'],ref_dir_n=geometry['ref_dir_n'],physical_origin_m=geometry['physical_origin_m'],
        ports={name:{key:p[key] for key in ['label','center','normal','lattice_center','opened_port_voxels','triangle_ids']} for name,p in geometry['ports'].items()},
        files=[dict(path=str(step2/name),sha256=hashlib.sha256((step2/name).read_bytes()).hexdigest()) for name in used])
    write(run/'contracts/geometry_reuse_contract.json',reuse)
    paths=[cfgpath,cfg3path,old/'configs/cfd_surface_prepare.yaml',constpath,runtimepath,promotedpath,s3path,s4path,
        old/'s3_cfd_1D_data_preprocess.py',old/'s4_cfd_surface_prepare.py',old/'s5_cfd_flow_solve.py',
        old/'utils/cfd_preprocess/pipeline.py',old/'utils/cfd_preprocess/port_transfer.py',
        old/'utils/cfd_surface_prepare/vmtk_pipeline.py',old/'utils/cfd_flow/production.py',old/'utils/cfd_flow/config.py',
        promoted/'qc/run_summary.json',promoted/'solver_smoke/musubi.lua',base/'segments/segment_0000000_to_0239502/musubi.lua']
    records=[]
    for p in paths:
        target=run/'provenance/source_records'/p.relative_to(old);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
        records.append(dict(path=str(p),snapshot=str(target.relative_to(run)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    write(run/'provenance/physical_source_records.json',records)
    print('PHYSICAL_BC_CONTRACT=PASS; geometry reused; rho=',rho,'nu=',nu,'Q=',q,'gauge pressures=',outlets,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);a=ap.parse_args();prepare(a.run)
