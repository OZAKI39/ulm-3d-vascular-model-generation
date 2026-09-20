#!/usr/bin/env python3
"""Freeze the native case and independently integrate the written inlet prescription."""
import json
import shutil
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import pyvista as pv
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from sv_validation.mesh_diagnostics import unique_edges
from sv_validation.geometry import ROLE_IDS, triangles_of
from sv_validation.validation import reference_condition, triangle_flux, normalized_inlet
from sv_validation.provenance import write_json, sha256, now
from sv_validation.visuals import geometry_figure


def add(parent, tag, value=None, **attributes):
    node = ET.SubElement(parent, tag, attributes)
    if value is not None:
        node.text = str(value).lower() if isinstance(value, bool) else str(value)
    return node


def main():
    require_path = ROOT / 'reports/sv1/mesh_validity.json'
    assert json.loads(require_path.read_text())['status'] == 'PASS'
    assert not (ROOT / 'configs/time_policy.json').exists(), 'Time policy already frozen'
    reference = yaml.safe_load((ROOT / 'configs/sv_reference.yaml').read_text())
    physical = reference_condition(reference)
    Q = physical['inlet_volume_flow_m3_s']
    mesh_dir = ROOT / 'outputs/sv1/SV_MESH'
    case = ROOT / 'outputs/sv1/vascular_flow'
    case.mkdir(parents=True, exist_ok=True)
    with np.load(mesh_dir / 'mesh_arrays.npz') as data:
        points, tetra, boundary, tags = [data[key] for key in ('points_m', 'tetra', 'boundary_triangles', 'facet_tags')]
    inlet = pv.read(mesh_dir / 'mesh-surfaces/INLET.vtp')
    ip = np.asarray(inlet.points, np.float64)
    it = triangles_of(inlet)
    xyz = ip[it]
    area_vectors = .5*np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])
    triangle_areas = np.linalg.norm(area_vectors, axis=1)
    A = float(triangle_areas.sum())
    wall_ids = set(np.unique(boundary[tags == 1]).tolist())
    global_ids = np.asarray(inlet.point_data['GlobalNodeID'], int)-1
    shape = np.array([0. if int(i) in wall_ids else 1. for i in global_ids])
    # Official Flat profile with shared wall-rim nodes zeroed before flux normalization.
    effective_area = float(np.sum(triangle_areas * shape[it].mean(axis=1)))
    nodal_normals = np.zeros_like(ip)
    for j in range(3):
        np.add.at(nodal_normals, it[:, j], area_vectors / 3.)
    nodal_normals /= np.linalg.norm(nodal_normals, axis=1)[:, None]
    unit_profile = -shape[:, None] * nodal_normals / effective_area
    alignment = -triangle_flux(ip, it, unit_profile)
    # Correct only the roundoff-level loss due to normals on the discrete nearly planar cap.
    scalar_flow_parameter = -Q / alignment
    prescribed_velocity = -scalar_flow_parameter * unit_profile
    prescribed_flux = -triangle_flux(ip, it, prescribed_velocity)
    normalized_inlet(prescribed_flux, Q)
    inlet.point_data['PrescribedVelocity'] = prescribed_velocity
    inlet.save(case / 'prescribed_inlet.vtp')
    geometry = json.loads((ROOT / 'reports/sv1/geometry_qc.json').read_text())
    Dh = geometry['ports']['INLET']['actual']['hydraulic_diameter_m']
    Umean = Q / A
    nu = physical['dynamic_viscosity_pa_s'] / physical['density_kg_m3']
    tnu = Dh**2 / nu
    edges = unique_edges(tetra)
    h10 = float(np.quantile(np.linalg.norm(points[edges[:, 1]]-points[edges[:, 0]], axis=1), .1))
    dt = min(.5*h10/Umean, .05*tnu)
    first_steps = min(1000, int(np.ceil(20.*tnu/dt-1e-12)))
    interval = max(1, first_steps//40)
    time_policy = {'timestamp': now(), 'mode': 'transient-to-steady', 'equation': 'fluid',
                   'dt_s': dt, 'dt_rule': 'min(0.5*h10/Umean, 0.05*Dh^2/nu)',
                   'h10_m': h10, 'A_in_m2': A, 'Umean_m_s': Umean, 'Dh_m': Dh, 'nu_m2_s': nu,
                   'viscous_time_s': tnu, 'target_block_duration_s': 20.*tnu,
                   'first_block_steps': first_steps, 'save_interval_steps': interval,
                   'maximum_extensions': 1, 'extension_steps': first_steps,
                   'steady_last_intervals': 5, 'velocity_change_limit': 1e-5, 'flow_change_limit': 1e-6,
                   'Re': physical['density_kg_m3']*Umean*Dh/physical['dynamic_viscosity_pa_s']}
    write_json(ROOT / 'configs/time_policy.json', time_policy)
    xml = ET.Element('svMultiPhysicsFile', version='0.1')
    general = add(xml, 'GeneralSimulationParameters')
    for key, value in {
        'Continue_previous_simulation': False, 'Number_of_spatial_dimensions': 3,
        'Number_of_time_steps': first_steps, 'Time_step_size': format(dt, '.17g'),
        'Spectral_radius_of_infinite_time_step': .5, 'Searched_file_name_to_trigger_stop': 'STOP_SIM',
        'Save_results_to_VTK_format': True, 'Name_prefix_of_saved_VTK_files': 'result',
        'Increment_in_saving_VTK_files': interval, 'Start_saving_after_time_step': 1,
        'Increment_in_saving_restart_files': interval, 'Convert_BIN_to_VTK_format': False,
        'Verbose': True, 'Warning': True, 'Debug': False}.items():
        add(general, key, value)
    mesh = add(xml, 'Add_mesh', name='vascular')
    add(mesh, 'Mesh_file_path', '../SV_MESH/mesh-complete.mesh.vtu')
    for role in ROLE_IDS.values():
        face = add(mesh, 'Add_face', name=role)
        add(face, 'Face_file_path', '../SV_MESH/mesh-surfaces/' + role + '.vtp')
    equation = add(xml, 'Add_equation', type='fluid')
    for key, value in {'Coupled': True, 'Min_iterations': 2, 'Max_iterations': 12,
                       'Tolerance': 1e-10, 'Backflow_stabilization_coefficient': 0.,
                       'Density': physical['density_kg_m3']}.items():
        add(equation, key, value)
    add(add(equation, 'Viscosity', model='Constant'), 'Value', physical['dynamic_viscosity_pa_s'])
    for output_type in ('Spatial', 'Boundary_integral'):
        output = add(equation, 'Output', type=output_type)
        add(output, 'Velocity', True)
        add(output, 'Pressure', True)
    linear = add(equation, 'LS', type='NS')
    add(add(linear, 'Linear_algebra', type='fsils'), 'Preconditioner', 'fsils')
    for key, value in {'Max_iterations': 30, 'Tolerance': 1e-10, 'Absolute_tolerance': 1e-24,
                       'Krylov_space_dimension': 100, 'NS_GM_max_iterations': 20,
                       'NS_CG_max_iterations': 1000, 'NS_GM_tolerance': 1e-3,
                       'NS_CG_tolerance': 1e-3}.items():
        add(linear, key, value)
    bc = add(equation, 'Add_BC', name='INLET')
    for key, value in {'Type': 'Dir', 'Time_dependence': 'Steady', 'Profile': 'Flat',
                       'Impose_flux': True, 'Zero_out_perimeter': True,
                       'Value': format(scalar_flow_parameter, '.17g')}.items():
        add(bc, key, value)
    for role in ('OUTLET_01', 'OUTLET_02', 'OUTLET_03'):
        bc = add(equation, 'Add_BC', name=role)
        for key, value in {'Type': 'Neu', 'Time_dependence': 'Steady', 'Value': 0.}.items():
            add(bc, key, value)
    bc = add(equation, 'Add_BC', name='WALL')
    for key, value in {'Type': 'Dir', 'Time_dependence': 'Steady', 'Value': 0.}.items():
        add(bc, key, value)
    ET.indent(xml, space='  ')
    ET.ElementTree(xml).write(case / 'solver.xml', encoding='utf-8', xml_declaration=True)
    shutil.copy2(case / 'solver.xml', ROOT / 'configs/sv_flow.xml')
    # Read both the prescription artifact and the XML back before accepting its flux.
    written = pv.read(case / 'prescribed_inlet.vtp')
    Q_written = -triangle_flux(written.points, triangles_of(written), written.point_data['PrescribedVelocity'])
    error = normalized_inlet(Q_written, Q)
    parsed = ET.parse(case / 'solver.xml')
    assert float(parsed.find(".//Add_BC[@name='INLET']/Value").text) == scalar_flow_parameter
    native = json.loads((ROOT / 'reports/sv1/native_solver.json').read_text())
    commit = native['commit']
    evidence = {'commit': commit, 'mode': 'transient-to-steady', 'equation_type': 'fluid',
                'evidence': [
                    'https://github.com/SimVascular/svMultiPhysics/blob/' + commit + '/Code/Source/solver/' + name
                    for name in ('Parameters.cpp', 'main.cpp', 'baf_ini.cpp', 'set_bc.cpp', 'read_files.cpp')],
                'pressure_convention': 'Three Neu Value=0 natural reference loads; no nodal pressure Dirichlet or pressure pin; outlet pressures are reported as area averages',
                'inlet_convention': 'Official Flat + Impose_flux + Zero_out_perimeter; constant interior magnitude, zero shared wall rim, discrete P1 transition at the rim',
                'normal_roundoff_alignment_factor': alignment,
                'input_bcs_time_dependence_is_not_a_steady_pde_switch': True,
                'solver_source_modified': False, 'source_case_for_builtin_linear_solver': 'tests/cases/fluid/iliac_artery/solver.xml'}
    write_json(ROOT / 'reports/sv1/solver_capabilities.json', evidence)
    write_json(ROOT / 'reports/sv1/inlet_normalization.json', {
        'status': 'PASS', 'Q_target_m3_s': Q, 'Q_written_m3_s': Q_written, 'relative_error': error,
        'actual_area_m2': A, 'Umean_m_s': Umean, 'effective_profile_area_m2': effective_area,
        'interior_speed_m_s': -scalar_flow_parameter/effective_area,
        'rim_nodes': int(np.count_nonzero(shape == 0)), 'inlet_nodes': len(shape),
        'native_scalar_flow_parameter': scalar_flow_parameter, 'alignment_factor': alignment,
        'artifact': str((case/'prescribed_inlet.vtp').relative_to(ROOT)),
        'artifact_is_a_boundary_prescription_not_a_solution': True,
        'written_xml_sha256': sha256(case/'solver.xml'), 'written_field_sha256': sha256(case/'prescribed_inlet.vtp')})
    write_json(ROOT / 'reports/sv1/flow_config_provenance.json', {
        'xml_sha256': sha256(ROOT/'configs/sv_flow.xml'),
        'time_policy_sha256': sha256(ROOT/'configs/time_policy.json'),
        'reference_physics_sha256': sha256(ROOT/'configs/sv_reference.yaml')})
    geometry_figure(points, boundary, tags, 'real_geometry_and_bc.png', '流体从哪里进入、从哪里流出？', bc=True)
    print(json.dumps({'time_policy': time_policy, 'Q_written': Q_written, 'normalization_error': error,
                      'rim_nodes': int(np.count_nonzero(shape == 0)), 'alignment_factor': alignment}, indent=2))


if __name__ == '__main__':
    main()
