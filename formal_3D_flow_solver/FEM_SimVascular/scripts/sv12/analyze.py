#!/usr/bin/env python3
"""Summarize the completed native run and integrate accepted pressure sections."""
import json
import os
import sys
from pathlib import Path
os.environ['LIBGL_ALWAYS_SOFTWARE'] = '1'
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
from sv_validation.sv12 import REPORT, OUTPUT, load
from sv_validation.provenance import write_json, sha256

flow = load('flow_execution')
runs = [load(name+'_execution') for name in flow['runs']]
seed = load('petsc_short_execution', 'sv1_1')
def compact(row, origin):
    return {**{k: v for k, v in row.items() if k != 'petsc_monitor'}, 'origin': origin}
linear = [compact(row, 'inherited SV1.1 steps 1–10') for row in seed['history']['linear_solves']]
nonlinear = []
for step in range(1, 11):
    rows = [r for r in seed['history']['linear_solves'] if r['step'] == step]
    first, last = rows[0], rows[-1]
    nonlinear.append({'step': step, 'time_s': last['time_s'], 'iterations': last['nonlinear_iteration'],
        'initial_residual': first['petsc_monitor'][0]['residual_norm'],
        'final_residual': last['petsc_monitor'][0]['residual_norm'],
        'Ri_over_R0': last['nonlinear_Ri_over_R0'], 'Ri_over_R1': last['nonlinear_Ri_over_R1'],
        'converged': last['nonlinear_iteration'] >= 2 and min(last['nonlinear_Ri_over_R0'], last['nonlinear_Ri_over_R1']) <= 1e-10,
        'origin': 'inherited SV1.1 steps 1–10'})
for name, run in zip(flow['runs'], runs):
    linear.extend(compact(row, 'SV1.2 '+name) for row in run['linear_solves'])
    nonlinear.extend({**row, 'origin': 'SV1.2 '+name} for row in run['nonlinear_history'])
def distribution(values):
    if not values:
        return {'count': 0, 'total': 0, 'mean': None, 'median': None, 'P95': None, 'max': None}
    return {'count': len(values), 'total': int(sum(values)), 'mean': float(np.mean(values)),
            'median': float(np.median(values)), 'P95': float(np.quantile(values, .95)), 'max': int(max(values))}
new_linear = [r for r in linear if r['step'] > 10]
new_nonlinear = [r for r in nonlinear if r['step'] > 10]
history = {
    'linear': linear, 'nonlinear': nonlinear,
    'linear_failures': sum(r['linear_failures'] for r in runs),
    'nonlinear_failures': sum(r['nonlinear_failures'] for r in runs),
    'ill_conditioned_warnings': sum(r['ill_conditioned_warnings'] for r in runs),
    'new_linear_iteration_distribution': distribution([r['linear_iterations'] for r in new_linear]),
    'complete_linear_iteration_distribution': distribution([r['linear_iterations'] for r in linear]),
    'new_nonlinear_iteration_distribution': distribution([r['iterations'] for r in new_nonlinear]),
    'complete_nonlinear_iteration_distribution': distribution([r['iterations'] for r in nonlinear]),
    'coverage_note': 'Steps 1–10 inherited, SHA-verified SV1.1 evidence; steps 11 onward newly computed in SV1.2. No benchmark repeated.'}
write_json(REPORT/'solver_history.json', history)
qc = load('saved_state_qc')
first_joint = next((qc['intervals'][i]['step'] for i in range(4, len(qc['intervals']))
                   if all(r['E_u'] <= 1e-5 and r['E_Q'] <= 1e-6
                          for r in qc['intervals'][i-4:i+1])), None)
first_mass = next((s['step'] for s in qc['states']
                   if s['epsilon_Q'] <= 1e-6 and s['epsilon_mass'] <= 1e-6), None)
write_json(REPORT/'steady_assessment.json', {
    'first_saved_state_passing_mass_step': first_mass,
    'first_five_joint_intervals_end_step': first_joint,
    'last_saved_step': qc['last_step'], 'final_joint_gate': qc['steady_reached'],
    'last_five_intervals': qc['intervals'][-5:],
    'consecutive_passing_intervals': qc['consecutive_passing_intervals'],
    'acceptance_note': 'First threshold crossing never stops the required full block. Only the final saved state is eligible for acceptance.'})
resource_rows = [{
    'name': name, 'solver_wall_s': run['GNU_wall_s'],
    'python_monotonic_s': run['elapsed_monotonic_s'],
    'peak_single_process_RSS_KiB': run['peak_single_process_RSS_KiB'],
    'native_reported_elapsed_initial_s': run['linear_solves'][0]['reported_elapsed_s'] if run['linear_solves'] else None,
    'native_reported_elapsed_final_s': run['linear_solves'][-1]['reported_elapsed_s'] if run['linear_solves'] else None,
    'raw_log': run['resource_log']} for name, run in zip(flow['runs'], runs)]
data_files = [p for p in (OUTPUT/'vascular_flow').rglob('*') if p.is_file()]
resources = {
    'measurements': resource_rows,
    'total_solver_wall_s': sum(r['solver_wall_s'] or 0 for r in resource_rows),
    'total_python_monotonic_s': sum(r['python_monotonic_s'] for r in resource_rows),
    'peak_single_process_RSS_KiB': max((r['peak_single_process_RSS_KiB'] or 0) for r in resource_rows),
    'memory_scope': 'GNU time maximum individual process RSS; not sum of MPI ranks',
    'clock_note': 'Native elapsed includes timer inherited from the native checkpoint; fresh GNU wall and Python monotonic clocks are both retained. No speed ratio is inferred.',
    'inherited_first_ten_python_wall_s': seed['elapsed_s'],
    'mpi_ranks': 4, 'OMP_NUM_THREADS': 1,
    'result_file_count': len(data_files), 'result_logical_size_bytes': sum(p.stat().st_size for p in data_files),
    'size_note': 'Logical case-file sizes including native checkpoints and copied seed; frozen mesh symlink excluded',
    'PETSc_iterations': history['new_linear_iteration_distribution'],
    'nonlinear_iterations': history['new_nonlinear_iteration_distribution']}
write_json(REPORT/'solver_resource_usage.json', resources)

if flow['accepted_solution_available']:
    import pyvista as pv
    final = load('accepted_solution')
    path = ROOT/final['path']
    assert sha256(path) == final['sha256']
    grid = pv.read(path)
    points = np.asarray(grid.points, float)
    center = points.mean(axis=0)
    _, _, vectors = np.linalg.svd(points-center, full_matrices=False)
    axis = vectors[0]
    if axis[np.argmax(abs(axis))] < 0:
        axis = -axis
    projected = (points-center)@axis
    limits = [float(projected.min()), float(projected.max())]
    sections = []
    for distance in np.linspace(*limits, 14)[1:-1]:
        cut = grid.slice(normal=axis, origin=center+distance*axis).triangulate()
        if cut.n_cells == 0:
            continue
        faces = cut.faces.reshape(-1, 4)
        assert np.all(faces[:, 0] == 3)
        tri = faces[:, 1:]
        xyz = np.asarray(cut.points, float)[tri]
        areas = .5*np.linalg.norm(np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0]), axis=1)
        values = np.asarray(cut['Pressure'], float)[tri].mean(axis=1)
        assert areas.sum() > 0 and np.isfinite(values).all()
        sections.append({'offset_m': float(distance), 'area_m2': float(areas.sum()),
                         'triangles': len(tri), 'mean_pressure_pa': float(np.dot(areas, values)/areas.sum())})
    assert len(sections) == 12
    write_json(REPORT/'pressure_sections.json', {
        'path': final['path'], 'sha256': final['sha256'],
        'port_area_averages_pa': final['area_average_pressure_pa'],
        'raw_pressure_range_pa': final['pressure_range_pa'], 'raw_pressure_offset_applied_pa': 0,
        'pressure_reference': 'Frozen zero-Neumann outlets; no postprocessing offset or pin',
        'axis': axis.tolist(), 'origin_m': center.tolist(), 'projection_limits_m': limits,
        'view_up': vectors[1].tolist(), 'sections': sections,
        'method': 'Area-weighted P1 pressure on geometric principal-axis planes. Planes may cross multiple branches; not a centerline pressure-drop curve.'})
    assert sha256(path) == final['sha256']
print('SV1.2 solver history, resources and accepted pressure analysis complete', flush=True)
