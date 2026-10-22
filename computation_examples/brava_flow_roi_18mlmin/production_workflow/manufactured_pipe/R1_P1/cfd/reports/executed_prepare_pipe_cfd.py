"""Prepare one gated pipe CFD case with production svMultiPhysics and CPU PETSc."""
from pathlib import Path
import argparse
import hashlib
import json
import shlex
import shutil
import xml.etree.ElementTree as ET

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(name, mesh_dir, manufactured=False):
    qa = json.loads((mesh_dir / 'NATIVE_MESH_QA.json').read_text())
    assert qa['status'] == 'PASS', qa['checks']
    proto = json.loads((ROOT / ('manufactured_pipe/PROTOCOL.json' if manufactured else 'straight_pipe/BENCHMARK_PROTOCOL.json')).read_text())
    inputs = json.loads((ROOT / ('straight_pipe/manufactured' if manufactured else 'straight_pipe') / name / 'INPUT.json').read_text())
    case = ROOT / 'manufactured_pipe' / name / 'cfd' if manufactured else ROOT / 'straight_pipe/calibration' / name
    if manufactured:
        fidelity = json.loads((mesh_dir / 'FIDELITY_QA.json').read_text())
        assert fidelity['status'] == 'PASS'
        assert qa['all_WALL_tetra_count'] == 0
        assert qa['tetra'] <= proto['maximum_tetra_per_pipe']
    case.mkdir(parents=True, exist_ok=False)
    run, reports, mesh = case / 'run', case / 'reports', case / 'SV_MESH'
    run.mkdir(); reports.mkdir(); mesh.mkdir()
    (mesh / 'mesh-surfaces').mkdir()
    shutil.copy2(mesh_dir / 'volume.vtu', mesh / 'mesh-complete.mesh.vtu')
    for role in ['WALL', 'INLET', 'OUTLET_01']:
        shutil.copy2(mesh_dir / (role + '.vtp'), mesh / 'mesh-surfaces' / (role + '.vtp'))
    with np.load(mesh_dir / 'BOUNDARY_TOPOLOGY.npz') as data:
        np.savez_compressed(mesh / 'mesh_arrays.npz', **{key: data[key] for key in ['points_m', 'tetra', 'boundary_triangles', 'facet_tags']})
    grid = pv.read(mesh / 'mesh-complete.mesh.vtu')
    assert 'GlobalNodeID' in grid.point_data and 'GlobalElementID' in grid.cell_data
    xml = ET.parse(CASE / 'cases/balanced_pressure_final/run/solver.xml')
    add_mesh = xml.find('.//Add_mesh')
    for face in list(add_mesh.findall('Add_face')):
        if face.get('name') not in ['WALL', 'INLET', 'OUTLET_01']:
            add_mesh.remove(face)
    equation = xml.find('.//Add_equation')
    for bc in list(equation.findall('Add_BC')):
        role = bc.get('name')
        if role not in ['WALL', 'INLET', 'OUTLET_01']:
            equation.remove(bc)
        elif role == 'INLET':
            bc.find('Profile').text = 'Parabolic'
            bc.find('Value').text = format(-inputs['analytic']['Q_m3_s'], '.17g')
        elif role == 'OUTLET_01':
            bc.find('Value').text = '0'
    if manufactured:
        from manufactured_boundary import prepare, configure_xml
        boundary = prepare(mesh, run, inputs, proto['mean_velocity_m_s'], proto['mu_Pa_s'])
        configure_xml(equation)
    ET.indent(xml, space='  ')
    xml.write(run / 'solver.xml', encoding='utf-8', xml_declaration=True)
    # Use the installed production binary; CPU matrix/vector types only.
    options = shlex.split((CASE / 'cases/balanced_pressure_final/run/PETSC_OPTIONS.txt').read_text())
    for flag, value in [('-mat_type', 'aij'), ('-vec_type', 'standard')]:
        options[options.index(flag) + 1] = value
    options = [item for item in options if item != '-log_view_gpu_time']
    (run / 'PETSC_OPTIONS.txt').write_text(' '.join(options) + '\n')
    policy = json.loads((CASE / 'cases/balanced_pressure_final/policy.json').read_text())
    inlet_area = qa['surface']['ports']['INLET']['cap_area_m2']
    policy.update(Q_target_m3_s=inputs['analytic']['Q_m3_s'], Umean_m_s=proto['mean_velocity_m_s'],
                  A_in_m2=inlet_area, MPI_ranks=8, maximum_wall_time_s=proto['maximum_CFD_s_per_case'],
                  resource_plan='8 existing CPU MPI ranks; PETSc aij/standard; production executable and tolerances unchanged',
                  mode='pipe_parabolic_flux_pressure_outlet')
    if manufactured:
        policy.update(Q_target_m3_s=boundary['Q_monitor_m3_s'], Q_analytical_m3_s=inputs['analytic']['Q_m3_s'],
                      mode='manufactured_exact_nodal_velocity_cauchy_traction', nodal_flux_scaling=False)
    (case / 'policy.json').write_text(json.dumps(policy, indent=2) + '\n')
    baseline = json.loads((CASE / 'cases/balanced_pressure_final/reports/execution.json').read_text())
    ratio = qa['points'] / 135852
    # Conservative startup overhead and a 3x empirical scaling envelope, not a guarantee.
    ram_gib = 1.5 + baseline['peak_tree_RSS_MiB_sampled'] / 1024 * ratio * 3
    seconds = 30 + baseline['wall_time_s'] * ratio * (100 / baseline['final_step']) * 3
    budget = proto['maximum_CFD_s_per_case']
    resource = dict(status='PASS' if ram_gib < 24 and seconds < budget else 'FAIL_ESTIMATED_BUDGET',
                    nodes=qa['points'], tetra=qa['tetra'], velocity_pressure_DOF_upper=4 * qa['points'],
                    estimated_tree_RAM_GiB=ram_gib, estimated_wall_time_s=seconds,
                    hard_wall_budget_s=budget, hard_tree_RSS_GiB=min(24, max(4, 2 * ram_gib)),
                    working_set_limit_GiB=37, baseline='cases/balanced_pressure_final/reports/execution.json',
                    method='Measured 135852-node/8-rank production run, node-linear extrapolation with 3x envelope and startup allowance; actual logs determine acceptance.')
    if manufactured:
        resource.pop('working_set_limit_GiB')
        resource.update(nonreclaimable_upper_limit_GiB=34,
                        guard_reason='Active file cache is reclaimable; retain 24 GiB job RSS cap and conservative anon/shmem/kernel/dirty upper estimate guard with OOM monitoring.')
    (reports / 'RESOURCE_PLAN.json').write_text(json.dumps(resource, indent=2) + '\n')
    assert resource['status'] == 'PASS', resource
    purpose = dict(physical_parameters_from='straight_pipe/BENCHMARK_PROTOCOL.json',
                   inlet='Steady parabolic profile normalized by native solver to exact analytical Q; nodal profile approximation retained as mesh error.',
                   outlet='Zero scalar pressure/traction, native natural outflow; analyze only interior z/L=0.4,0.5,0.6.',
                   wall='Rigid no-slip', end_effect_check='Compare all three interior planes; no end-plane solution is used for calibration. Natural outflow is not an exact full viscous-traction prescription.',
                   native_mesh_sha256=sha(mesh / 'mesh-complete.mesh.vtu'),
                   production_runner_sha256=sha(CASE / 'scripts/solve_remote.py'),
                   unchanged_formal_human_BC=True, initial_state='zero')
    if manufactured:
        purpose.update(inlet='Analytical u_z(x,y) at native physical coordinates; rim zero; General vector Dirichlet; no flux normalization',
                       outlet='Full analytical Cauchy traction VTP, not scalar natural zero traction',
                       end_effect_check='Exact continuum inlet/traction data makes the analytical finite circular-pipe solution admissible; residual polygon and trace errors are independently measured before CFD.',
                       boundary_fidelity='reports/BOUNDARY_FIDELITY.json',
                       conservation_reference='Discrete prescribed inlet trace Q; analytical Q still used in scientific calibration')
    (case / 'purpose.json').write_text(json.dumps(purpose, indent=2) + '\n')
    files = [run / 'solver.xml', run / 'PETSC_OPTIONS.txt', case / 'policy.json'] + [p for p in mesh.rglob('*') if p.is_file()]
    if manufactured:
        files += [run / 'inlet_exact_velocity.dat', run / 'outlet_exact_traction.vtp']
    (case / 'input_hashes.json').write_text(json.dumps({str(p.relative_to(case)): sha(p) for p in files}, indent=2) + '\n')
    print(json.dumps(dict(case=str(case), resource=resource, boundary_conditions=purpose), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case')
    parser.add_argument('mesh_dir', type=Path)
    parser.add_argument('--manufactured', action='store_true')
    args = parser.parse_args()
    main(args.case, args.mesh_dir, args.manufactured)
