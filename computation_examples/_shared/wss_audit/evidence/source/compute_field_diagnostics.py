"""Read-only diagnostics of the frozen 2 mm/s FEM solution and actual KSP log.

WSS is derived from P1 element gradients on WALL triangles only. Raw facet
values and an explicitly identified area-weighted nodal display field coexist.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / 'flow_cases/mean-2p0-mmps'
OUT = CASE / 'field_diagnostics'
sys.path.insert(0, str(ROOT / 'scripts/sv13q'))
from flow_parser import parse_solver_log


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def write_csv(path, rows):
    with Path(path).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def p1_gradients(points, tetra, velocity):
    x = points[tetra]
    u = velocity[tetra]
    return np.swapaxes(np.linalg.solve(x[:, 1:] - x[:, :1], u[:, 1:] - u[:, :1]), 1, 2)


def face_keys(faces):
    a = np.ascontiguousarray(np.sort(faces, axis=1), dtype=np.int64)
    return a.view(np.dtype([('a', '<i8'), ('b', '<i8'), ('c', '<i8')])).ravel()


def boundary_owners(tetra, triangles):
    faces = tetra[:, [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]].reshape(-1, 3)
    keys, first, counts = np.unique(face_keys(faces), return_index=True, return_counts=True)
    query = face_keys(triangles)
    index = np.searchsorted(keys, query)
    assert np.all(index < len(keys)), 'Boundary face absent from tetrahedra'
    assert np.array_equal(keys[index], query)
    assert np.all(counts[index] == 1), 'Boundary triangles must have exactly one adjacent element'
    return first[index] // 4


def wall_geometry(points, tetra, triangles, owners):
    xyz = points[triangles]
    cross = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    norm = np.linalg.norm(cross, axis=1)
    assert np.all(norm > 0)
    normal = cross / norm[:, None]
    center = xyz.mean(axis=1)
    direction = center - points[tetra[owners]].mean(axis=1)
    normal[np.einsum('ij,ij->i', normal, direction) < 0] *= -1
    assert np.all(np.einsum('ij,ij->i', normal, direction) > 0)
    return center, norm / 2, normal


def tangential_traction(gradient, normal, viscosity):
    """Tangential part of viscous stress times outward fluid normal [Pa]."""
    traction = viscosity * np.einsum('nij,nj->ni', gradient + gradient.swapaxes(1, 2), normal)
    return traction - np.einsum('ij,ij->i', traction, normal)[:, None] * normal


def nodal_average(triangles, values, areas, npoints):
    weight = np.bincount(triangles.ravel(), weights=np.repeat(areas, 3), minlength=npoints)
    total = np.bincount(triangles.ravel(), weights=np.repeat(areas * values, 3), minlength=npoints)
    result = np.zeros(npoints)
    np.divide(total, weight, out=result, where=weight > 0)
    return result, weight


def surface(points, triangles):
    ids, inv = np.unique(triangles, return_inverse=True)
    faces = np.column_stack((np.full(len(triangles), 3), inv.reshape(-1, 3)))
    mesh = pv.PolyData(points[ids], faces)
    mesh.point_data['GlobalNodeID_zero_based'] = ids
    return mesh, ids


def make_residual_figures(history, summary, steps, rtol, atol):
    # One source-audited plotter serves both full diagnostics and plot-only edits.
    from residual_figures import render_residuals
    return render_residuals(history, summary, steps, rtol, atol)


def main():
    start = time.time()
    for name in ['data', 'figures', 'animations', 'inspection']:
        (OUT / name).mkdir(parents=True, exist_ok=True)
    sources = ['frozen_flow/flow_arrays_si.npz', 'frozen_flow/steady_flow_mean_2p0_mmps.vtu',
        'frozen_flow/manifest.json', 'run/solver.xml', 'run/solver.log', 'run/PETSC_OPTIONS.txt',
        'reports/execution.json', 'reports/physics_validation.json', 'streamlines/data/streamlines_si.vtp',
        'figures/Figure_06_streamlines_overview.png', 'animations/Animation_03_streamlines_full_vessel.mp4']
    lock = {name: sha(CASE / name) for name in sources}
    if (OUT / 'SOURCE_LOCK.json').exists():
        assert lock == json.loads((OUT / 'SOURCE_LOCK.json').read_text())
    else:
        dump(OUT / 'SOURCE_LOCK.json', lock)
    tree = ET.parse(CASE / 'run/solver.xml')
    dt = float(tree.find('.//Time_step_size').text)
    viscosity = float(tree.find('.//Viscosity/Value').text)
    data = np.load(CASE / 'frozen_flow/flow_arrays_si.npz')
    points, tetra, triangles, tags, velocity, pressure = [data[k] for k in
        ['points_m', 'tetra', 'boundary_triangles', 'facet_tags', 'velocity_m_s', 'pressure_pa']]
    wall_ids = np.flatnonzero(tags == 1)
    wall = triangles[wall_ids]
    owners = boundary_owners(tetra, wall)
    centers, area, normal = wall_geometry(points, tetra, wall, owners)
    gradients = p1_gradients(points, tetra, velocity)
    traction = tangential_traction(gradients[owners], normal, viscosity)
    magnitude = np.linalg.norm(traction, axis=1)
    display, weight = nodal_average(wall, magnitude, area, len(points))
    mesh, ids = surface(points, wall)
    mesh.point_data['WSS_display_Pa'] = display[ids]
    mesh.point_data['Pressure_Pa'] = pressure[ids]
    for key, value in {'WSS_raw_Pa': magnitude, 'Tangential_viscous_traction_Pa': traction,
        'Outward_normal': normal, 'Area_m2': area, 'Parent_tetra_zero_based': owners,
        'Global_boundary_facet_zero_based': wall_ids}.items():
        mesh.cell_data[key] = value
    mesh.save(OUT / 'data/wall_wss_si.vtp')
    pressure_surface, pressure_ids = surface(points, triangles)
    pressure_surface.point_data['Pressure_Pa'] = pressure[pressure_ids]
    pressure_surface.cell_data['Boundary_tag'] = tags
    pressure_surface.save(OUT / 'data/pressure_surface_si.vtp')
    rows = [dict(facet_id=int(fid), parent_tetra=int(owner), x_m=float(c[0]), y_m=float(c[1]), z_m=float(c[2]),
        area_m2=float(a), wss_Pa=float(m), traction_x_Pa=float(t[0]), traction_y_Pa=float(t[1]), traction_z_Pa=float(t[2]),
        normal_x=float(n[0]), normal_y=float(n[1]), normal_z=float(n[2]))
        for fid, owner, c, a, m, t, n in zip(wall_ids, owners, centers, area, magnitude, traction, normal)]
    write_csv(OUT / 'data/wall_wss_facets.csv', rows)
    # Independent affine 4x4 solves check the gradient layout and SI units.
    selected = np.unique(owners)[np.linspace(0, len(np.unique(owners))-1, 128, dtype=int)]
    affine = np.concatenate([np.ones((len(selected), 4, 1)), points[tetra[selected]]], axis=2)
    independent = np.swapaxes(np.linalg.solve(affine, velocity[tetra[selected]])[:, 1:], 1, 2)
    gradient_error = float(np.max(np.abs(independent - gradients[selected])))
    assert np.allclose(independent, gradients[selected], rtol=2e-10, atol=1e-7)
    tangential_error = float(np.max(np.abs(np.einsum('ij,ij->i', traction, normal))))
    assert tangential_error < max(1., magnitude.max()) * 1e-12
    assert np.isfinite(magnitude).all() and np.all(weight[ids] > 0)
    assert np.max(np.linalg.norm(velocity[ids], axis=1)) == 0
    history = parse_solver_log((CASE / 'run/solver.log').read_text(), dt)
    saved = json.loads((CASE / 'reports/execution.json').read_text())['history']
    solves = history['linear_solves']
    assert len(solves) == len(saved['linear_solves']) == 167
    assert not history['unparsed_rows'] and not history['failed_linear_solves'] and history['recovered_attempts'] == 0
    assert not history['unassigned_petsc_monitor'] and not history['nonlinear_failure_messages']
    options = (CASE / 'run/PETSC_OPTIONS.txt').read_text().split()
    rtol, atol = [float(options[options.index(k)+1]) for k in ['-ksp_rtol', '-ksp_atol']]
    summary, monitors = [], []
    for solve, original in zip(solves, saved['linear_solves']):
        assert solve['step'] == original['step'] and solve['petsc_monitor'] == original['petsc_monitor']
        ms = solve['petsc_monitor']; first, last = ms[0], ms[-1]
        assert all(m['finite'] for m in ms) and last['iteration'] == solve['linear_iterations']
        tolerance = max(atol, rtol * first['true_residual_norm'])
        record = dict(linear_solve_index=solve['linear_solve_index'], step=solve['step'], time_s=solve['time_s'],
            nonlinear_iteration=solve['nonlinear_iteration'], log_line=solve['log_line'],
            nonlinear_Ri_over_R0=solve['nonlinear_Ri_over_R0'], nonlinear_Ri_over_R1=solve['nonlinear_Ri_over_R1'],
            linear_iterations=solve['linear_iterations'], reason=solve['petsc_reason']['reason'],
            initial_true_residual=first['true_residual_norm'], final_true_residual=last['true_residual_norm'],
            final_true_relative_residual=last['true_relative_residual'], effective_absolute_tolerance=tolerance,
            final_true_over_effective_tolerance=last['true_residual_norm']/tolerance)
        summary.append(record)
        for m in ms:
            monitors.append(dict(linear_solve_index=solve['linear_solve_index'], step=solve['step'],
                nonlinear_iteration=solve['nonlinear_iteration'], **m))
    assert all(r['final_true_over_effective_tolerance'] <= 1.01 for r in summary)
    steps = []
    for step in sorted(set(r['step'] for r in summary)):
        rs = [r for r in summary if r['step'] == step]
        steps.append(dict(step=step, time_s=step*dt, nonlinear_iterations=len(rs),
            first_Ri_over_R0=rs[0]['nonlinear_Ri_over_R0'], last_Ri_over_R0=rs[-1]['nonlinear_Ri_over_R0']))
    write_csv(OUT / 'data/residual_linear_solves.csv', summary)
    write_csv(OUT / 'data/residual_true_monitor.csv', monitors)
    write_csv(OUT / 'data/residual_time_steps.csv', steps)
    make_residual_figures(history, summary, steps, rtol, atol)
    validation = dict(all_pass=True, source_hashes_unchanged=all(sha(CASE/k)==v for k,v in lock.items()),
        case='mean-2p0-mmps', dt_s=dt, last_step=steps[-1]['step'], frozen_time_s=steps[-1]['time_s'],
        viscosity_Pa_s=viscosity, nodes=len(points), tetrahedra=len(tetra), wall_facets=len(wall),
        pressure_Pa=dict(min=float(pressure.min()), max=float(pressure.max())),
        wss_Pa=dict(raw_min=float(magnitude.min()), raw_max=float(magnitude.max()),
            area_weighted_mean=float(np.sum(area*magnitude)/np.sum(area)),
            display_min=float(display[ids].min()), display_max=float(display[ids].max())),
        wss_method='Exact P1 tetra velocity gradients; tangential viscous traction on outward planar WALL facets',
        wss_display='Area-weighted mean of incident facet magnitudes at original wall nodes; geometry unchanged',
        wss_vector_convention='(I-n*n^T)*mu*(grad(u)+grad(u)^T)*n; svMultiPhysics bpost outputs its negative; magnitudes agree',
        wss_not_native_saved_solver_output=True, mesh_convergence_not_assessed=True,
        independent_gradient_max_absolute_error_s_inv=gradient_error, tangency_max_error_Pa=tangential_error,
        unique_parent_for_all_wall_facets=True, max_wall_speed_m_s=0.,
        residuals=dict(linear_solves=len(summary), monitor_rows=len(monitors), time_steps=len(steps),
            reasons=dict(Counter(r['reason'] for r in summary)), rtol=rtol, atol=atol,
            final_nonlinear_Ri_over_R0=steps[-1]['last_Ri_over_R0'],
            max_true_over_effective_tolerance=max(r['final_true_over_effective_tolerance'] for r in summary),
            log_matches_saved_history=True, recovered_attempts=0, failed_linear_solves=0,
            linear_residuals_are_of_solver_scaled_system=True),
        elapsed_seconds=time.time()-start, script_sha256=sha(__file__))
    assert validation['source_hashes_unchanged']
    generated = list((OUT/'data').iterdir()) + list((OUT/'figures').glob('residual_*'))
    validation['outputs_sha256'] = {str(p.relative_to(OUT)):sha(p) for p in generated if p.is_file()}
    dump(OUT/'COMPUTE_VALIDATION.json', validation)
    print(json.dumps({k:v for k,v in validation.items() if k!='outputs_sha256'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
