"""Compare a converged native field at the frozen exact-cylinder probe positions."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / 'vendor'))
from pipe_measurements import PipeMeasurements
from poiseuille import reference
from pipe_probe_geometry import frozen_probes, strict_cells


def error(actual, exact):
    return float(np.linalg.norm(actual - exact) / np.linalg.norm(exact))


def main(name, manufactured=False):
    case = ROOT / 'manufactured_pipe' / name / 'cfd' if manufactured else ROOT / 'straight_pipe/calibration' / name
    execution = json.loads((case / 'reports/execution.json').read_text())
    proto = json.loads((ROOT / ('manufactured_pipe/PROTOCOL.json' if manufactured else 'straight_pipe/BENCHMARK_PROTOCOL.json')).read_text())
    inputs = json.loads((ROOT / ('straight_pipe/manufactured' if manufactured else 'straight_pipe') / name / 'INPUT.json').read_text())
    out = case / 'analysis'
    out.mkdir(exist_ok=True)
    if execution['status'] != 'PASS':
        record = dict(status='FAIL_CFD_CONVERGENCE', case=name, solver_status=execution['status'],
                      reason='Do not certify errors from an unconverged or failed native solve.')
        (out / 'CALIBRATION.json').write_text(json.dumps(record, indent=2) + '\n')
        print(json.dumps(record, indent=2))
        return
    radius, length = inputs['radius_m'], inputs['length_m']
    mean, mu = proto['mean_velocity_m_s'], proto['mu_Pa_s']
    measure = PipeMeasurements(case / 'SV_MESH/mesh_arrays.npz', inputs['analytic']['Q_m3_s'], mean)
    final_path = case / execution['final_vtu']
    velocity, pressure = measure.read(final_path)
    grid = pv.read(case / 'SV_MESH/mesh-complete.mesh.vtu')
    points, tetra = measure.points, measure.tetra
    dx = points[tetra[:, 1:]] - points[tetra[:, :1]]
    du = velocity[tetra[:, 1:]] - velocity[tetra[:, :1]]
    gradient = np.linalg.solve(dx, du).transpose(0, 2, 1)
    levels, probes, groups, planes, azimuths = frozen_probes(proto, radius, length)
    angles = np.unique(azimuths)
    cells, locator_cells = strict_cells(grid, points, tetra, probes)
    valid = cells >= 0
    sampled_u = np.full((len(probes), 3), np.nan)
    sampled_grad = np.full(len(probes), np.nan)
    selected = cells[valid]
    bary = np.linalg.solve(dx[selected].transpose(0, 2, 1), (probes[valid] - points[tetra[selected, 0]])[..., None])[..., 0]
    weights = np.c_[1 - bary.sum(axis=1), bary]
    sampled_u[valid] = np.einsum('ni,nij->nj', weights, velocity[tetra[selected]])
    er = np.c_[np.cos(azimuths), np.sin(azimuths), np.zeros(len(probes))]
    sampled_grad[valid] = -np.einsum('ni,ni->n', gradient[selected, 2, :], er[valid])
    exact = reference(probes, radius, length, mean, mu)
    exact_u = exact['u_m_s']
    exact_grad = -exact['du_dr_s_inv']
    if manufactured:
        # Diagnostic interpolation error is kept separate from the actual CFD field.
        nodal_exact = np.zeros_like(velocity)
        nodal_exact[:,2] = 2 * mean * (1 - np.sum(points[:,:2]**2,axis=1)/radius**2)
        nodal_exact[measure.wall_nodes] = 0.  # Exact physical rim/wall value, up to coordinate storage roundoff.
        gradient_I = np.linalg.solve(dx, nodal_exact[tetra[:,1:]]-nodal_exact[tetra[:,:1]]).transpose(0,2,1)
        interpolated = np.full_like(sampled_u,np.nan)
        interpolated[valid] = np.einsum('ni,nij->nj',weights,nodal_exact[tetra[selected]])
        interpolated_grad = np.full(len(probes),np.nan)
        interpolated_grad[valid] = -np.einsum('ni,ni->n',gradient_I[selected,2,:],er[valid])
    # Raw traction from each adjacent tetra; no nodal/display smoothing.
    all_faces = np.concatenate([tetra[:, f] for f in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]])
    owners = np.tile(np.arange(len(tetra)), 4)
    lookup = {tuple(sorted(face)): int(owner) for face, owner in zip(all_faces, owners)}
    wall = measure.boundary[measure.tags == 1]
    wall_owner = np.array([lookup[tuple(sorted(face))] for face in wall])
    xyz = points[wall]
    normals = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    areas = np.linalg.norm(normals, axis=1) / 2
    normals /= (2 * areas[:, None])
    stress = mu * (gradient[wall_owner] + gradient[wall_owner].transpose(0, 2, 1))
    traction = np.einsum('nij,nj->ni', stress, normals)
    tangent = traction - np.einsum('ni,ni->n', traction, normals)[:, None] * normals
    wss = np.linalg.norm(tangent, axis=1)
    wall_poly = pv.PolyData(points, np.c_[np.full(len(wall), 3), wall].ravel())
    wall_targets = np.c_[radius * np.cos(azimuths), radius * np.sin(azimuths), probes[:, 2]]
    nearest_wall = wall_poly.find_closest_cell(wall_targets)
    probe_wss = wss[nearest_wall]
    exact_wss = np.full(len(probes), float(exact['WSS_Pa']))
    records = []
    for index, distance in enumerate(levels):
        mask = groups == index
        complete = bool(np.all(valid[mask]))
        record = dict(wall_distance_m=distance, samples=int(mask.sum()), inside=int(valid[mask].sum()),
                      complete_coverage=complete,
                      velocity_relative_L2=error(sampled_u[mask], exact_u[mask]) if complete else None,
                      gradient_relative_L2=error(sampled_grad[mask], exact_grad[mask]) if complete else None,
                      raw_wall_WSS_relative_L2=error(probe_wss[mask], exact_wss[mask]),
                      velocity_PASS=False, gradient_PASS=False)
        if complete:
            record['velocity_PASS'] = record['velocity_relative_L2'] <= proto['velocity_relative_L2_limit']
            record['gradient_PASS'] = record['gradient_relative_L2'] <= proto['gradient_WSS_relative_L2_limit']
        if manufactured:
            role = proto['probe_roles'][index]
            record.update(role=role, hard_gate=role!='WALL_STRESS_DIAGNOSTIC',
                          velocity_is_hard_gate=role in ('PARTICLE_RELEVANT','BULK'),
                          gradient_is_hard_gate=role=='PARTICLE_RELEVANT',
                          analytical_nodal_interpolation_velocity_error=error(interpolated[mask],exact_u[mask]) if complete else None,
                          CFD_minus_interpolant_velocity_relative=float(np.linalg.norm(sampled_u[mask]-interpolated[mask])/np.linalg.norm(exact_u[mask])) if complete else None,
                          analytical_nodal_interpolation_gradient_error=error(interpolated_grad[mask],exact_grad[mask]) if complete else None)
        records.append(record)
    # Integrate actual native cut triangles, never fit an assumed parabolic profile.
    grid.point_data['calibration_velocity'] = velocity
    grid.point_data['calibration_pressure'] = pressure
    sections = []
    for z in proto['analysis_z_over_L']:
        section = grid.slice(normal=[0, 0, 1], origin=[0, 0, z * length]).triangulate()
        tri = section.faces.reshape(-1, 4)[:, 1:]
        vertices = section.points[tri]
        area = np.linalg.norm(np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]), axis=1) / 2
        flow = float(np.dot(area, section['calibration_velocity'][tri, 2].mean(axis=1)))
        average = float(np.dot(area, section['calibration_pressure'][tri].mean(axis=1)) / area.sum())
        sections.append(dict(z_over_L=z, z_m=z * length, area_m2=float(area.sum()), Q_m3_s=flow,
                             Q_relative_error=abs(flow / inputs['analytic']['Q_m3_s'] - 1), pressure_mean_Pa=average))
    slope, intercept = np.polyfit([s['z_m'] for s in sections], [s['pressure_mean_Pa'] for s in sections], 1)
    pressure_gradient_error = abs(float(slope) / inputs['analytic']['pressure_gradient_Pa_m'] - 1)
    # Interior consistency measures at 0.25R and 0.5R; near-wall bands stay separate.
    central_variation = []
    for index in [len(levels) - 2, len(levels) - 1]:
        mask = groups == index
        values = sampled_u[mask].reshape(3, -1, 3)
        central_variation.append(dict(wall_distance_m=levels[index],
                                      maximum_plane_difference_relative=float(max(np.linalg.norm(values[i] - values[j]) / np.linalg.norm(exact_u[mask][:len(angles)]) for i in range(3) for j in range(i)))))
    interior_pass = all(s['maximum_plane_difference_relative'] <= .05 for s in central_variation)
    wss_error = records[0]['raw_wall_WSS_relative_L2']
    gates = dict(solver=execution['status'] == 'PASS', complete_probe_coverage=bool(np.all(valid)),
                 velocity=all(row['velocity_PASS'] for row in records), gradient=all(row['gradient_PASS'] for row in records),
                 WSS=wss_error <= .08, Q=max(s['Q_relative_error'] for s in sections) <= .05,
                 pressure_gradient=pressure_gradient_error <= .08, central_section_consistency=interior_pass)
    if manufactured:
        geometry = json.loads((case.parent/'mesh/FIDELITY_QA.json').read_text())
        boundary = json.loads((case/'reports/BOUNDARY_FIDELITY.json').read_text())
        particle_rows = [row for row in records if row['role']=='PARTICLE_RELEVANT']
        bulk_rows = [row for row in records if row['role']=='BULK']
        from audit_pipe_fidelity import boundary_profile
        inlet = pv.read(case/'SV_MESH/mesh-surfaces/INLET.vtp')
        _, expected_inlet = boundary_profile(inlet,radius,mean)
        index = {int(gid):i for i,gid in enumerate(grid.point_data['GlobalNodeID'])}
        inlet_nodes = np.array([index[int(gid)] for gid in inlet.point_data['GlobalNodeID']])
        expected_vector = np.c_[np.zeros((len(inlet_nodes),2)),expected_inlet]
        inlet_error = float(abs(velocity[inlet_nodes]-expected_vector).max())
        gates = dict(solver=execution['status']=='PASS',surface_fidelity=geometry['surface']['PASS'],
                     boundary_fidelity=boundary['status']=='PASS',applied_inlet_trace=inlet_error<=5e-8,
                     particle_probe_coverage=all(row['complete_coverage'] for row in particle_rows),
                     particle_velocity=all(row['velocity_PASS'] for row in particle_rows),
                     particle_gradient=all(row['gradient_PASS'] for row in particle_rows),
                     WSS=wss_error<=.08,internal_Q=max(s['Q_relative_error'] for s in sections)<=.05,
                     pressure_gradient=pressure_gradient_error<=.08,bulk_velocity=all(row['velocity_PASS'] for row in bulk_rows))
    summary = dict(status='PASS' if all(gates.values()) else 'FAIL_ENGINEERING_CALIBRATION', case=name, gates=gates,
                   samples=records, interior_sections=sections, central_section_variation=central_variation,
                   fitted_dp_dz_Pa_m=float(slope), exact_dp_dz_Pa_m=inputs['analytic']['pressure_gradient_Pa_m'],
                   pressure_gradient_relative_error=pressure_gradient_error, raw_WSS_relative_L2=wss_error,
                   tetra=len(tetra), nodes=len(points), CFD_wall_time_s=execution['wall_time_s'],
                   locator_tolerance_false_positives_rejected=int(np.sum((locator_cells >= 0) & (cells < 0))),
                   locator_interior_neighbor_corrections=int(np.sum((cells >= 0) & (cells != locator_cells))),
                   containment='VTK candidate locator followed by barycentric bounds with 64*double-epsilon*condition(DX) roundoff allowance; never extrapolate outside a tetra.',
                   CFD_peak_tree_RSS_MiB=execution.get('peak_tree_RSS_MiB_sampled'),
                   source_sha256={str(final_path.relative_to(ROOT)): hashlib.sha256(final_path.read_bytes()).hexdigest(),
                                  str((case / 'SV_MESH/mesh-complete.mesh.vtu').relative_to(ROOT)): hashlib.sha256((case / 'SV_MESH/mesh-complete.mesh.vtu').read_bytes()).hexdigest()},
                   missing_probe_policy='Outside probes retain NaN, full group error is null and coverage gate FAIL; no mesh-relative shifts or deletion.',
                   WSS_method='Raw adjacent-tetra symmetric viscous traction projected onto actual wall facets; nearest facet to each exact-cylinder wall target.',
                   acceptance='Engineering 5% velocity/Q, 8% gradient/WSS/dp-dz; 5% central bulk plane consistency. Not a universal physiological standard.')
    if manufactured:
        centers = points[tetra].mean(axis=1)
        central = (centers[:,2]>=.4*length)&(centers[:,2]<=.6*length)
        volumes = np.abs(np.linalg.det(dx))/6
        divergence = np.trace(gradient,axis1=1,axis2=2)
        rms = float(np.sqrt(np.dot(volumes[central],divergence[central]**2)/volumes[central].sum()))
        summary.update(benchmark='MANUFACTURED_PIPE_V1',probe_roles=proto['probe_roles'],wall_stress_hard_gate=False,
                       applied_inlet_profile_max_error_m_s=inlet_error,
                       div_diagnostic=dict(central_volume_RMS_s_inv=rms,normalized_by_U_over_R=rms*radius/mean,
                                           maximum_absolute_s_inv=float(abs(divergence[central]).max()),hard_gate=False),
                       analytical_Q_m3_s=inputs['analytic']['Q_m3_s'],prescribed_discrete_Q_m3_s=boundary['Q_monitor_m3_s'],
                       interpolation_diagnostic='u_h-u_exact=(u_h-Iu_exact)+(Iu_exact-u_exact) pointwise. Norms do not add. The interpolant is diagnostic only and never replaces CFD.',
                       acceptance='Particle-relevant u<=5%, G/WSS<=8%, internal Q<=5%, dp/dz<=8%, bulk u<=5%. 0.1um, div and cross-plane variation are diagnostics only.')
    (out / 'CALIBRATION.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    np.savez_compressed(out / 'FIXED_PROBES.npz', points_m=probes, cell=cells, locator_cell=locator_cells, groups=groups, planes=planes,
                        velocity=sampled_u, exact_velocity=exact_u, du_dn=sampled_grad, exact_du_dn=exact_grad,
                        raw_WSS_Pa=probe_wss, exact_WSS_Pa=exact_wss, gradient_all_tetra=gradient)
    if manufactured:
        np.savez_compressed(out/'INTERPOLATION_AND_DIV_DIAGNOSTICS.npz',interpolated_exact_velocity=interpolated,
                            interpolated_exact_du_dn=interpolated_grad,divergence_s_inv=divergence,
                            central_tetra=central,volume_m3=volumes)
    with (out / 'probe_error_groups.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    plt.rcParams.update({'font.size':11, 'figure.constrained_layout.use':True, 'pdf.fonttype':42})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.5))
    x = np.array(levels) * 1e6
    for key, label in [('velocity_relative_L2', 'Velocity'), ('gradient_relative_L2', 'Gradient')]:
        y = np.array([np.nan if row[key] is None else row[key] * 100 for row in records])
        ax[0].plot(x, y, 'o-', label=label)
    ax[0].axhline(5, ls='--', c='.5', label='u gate: 5%'); ax[0].axhline(8, ls=':', c='.5', label='G gate: 8%')
    if manufactured:
        ax[0].axvspan(min(x)/1.2, .2, color='.8', alpha=.4)
        ax[0].text(.02,.97,'0.1 um: diagnostic only',transform=ax[0].transAxes,va='top',fontsize=9)
    ax[0].set(xscale='log', xlim=(min(x) / 1.2, max(x) * 1.2), xlabel='Exact wall distance (um)', ylabel='Relative L2 error (%)', title='CFD error; fixed physical probes'); ax[0].legend(loc='best')
    ax[1].bar([format(v, '.3g') for v in x], [row['inside'] / row['samples'] * 100 for row in records], color='#367e9e')
    ax[1].set(xlabel='Wall distance (um)', ylabel='Inside original native volume (%)', ylim=(0, 105), title='Fixed analytical probe coverage'); ax[1].tick_params(axis='x', rotation=45)
    ax[2].bar(['WSS', 'Q (worst)', 'dp/dz'], [wss_error * 100, max(s['Q_relative_error'] for s in sections) * 100, pressure_gradient_error * 100], color='#8a6d9e')
    ax[2].set(ylabel='Relative error (%)', title='Interior physical quantities')
    fig.suptitle(('Manufactured ' if manufactured else '') + name + ' | ' + summary['status'])
    for ext in ['png', 'pdf']:
        fig.savefig(ROOT / ('figures/' + ('manufactured_' if manufactured else 'poiseuille_') + name + '.' + ext), dpi=300)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('case'); parser.add_argument('--manufactured',action='store_true')
    args = parser.parse_args(); main(args.case,args.manufactured)
