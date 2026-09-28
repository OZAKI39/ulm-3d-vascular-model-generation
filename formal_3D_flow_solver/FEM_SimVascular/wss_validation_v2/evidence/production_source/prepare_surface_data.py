"""Recover production P1 WSS from an explicitly identified SI flow and case XML."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
import pyvista as pv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'solver_support/src'))
from flow_solver_support import wss
from flow_solver_support.wss_case import material, recover, sha, coordinate_identity


def main():
    started = time.perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, help='Case owning the solver configuration and policy')
    parser.add_argument('--config', type=Path, help='Default: CASE/run/solver.xml')
    parser.add_argument('--flow', type=Path, help='Explicit solver VTU; default: case frozen VTU')
    parser.add_argument('--arrays', type=Path, help='Default: CASE/frozen_flow/flow_arrays_si.npz')
    parser.add_argument('--mesh', type=Path, help='Default: CASE/SV_MESH/mesh_arrays.npz')
    parser.add_argument('--output', type=Path, default=HERE/'input_data/field_diagnostics')
    args = parser.parse_args()
    context = json.loads((HERE/'BUILD_CONTEXT.json').read_text()) if args.case is None else None
    case = (args.case or Path(context['source_case'])).resolve()
    flow = args.flow
    if flow is None:
        candidates = sorted((case/'frozen_flow').glob('*.vtu'))
        if len(candidates) != 1:
            raise ValueError('Specify --flow when frozen_flow does not contain exactly one VTU')
        flow = candidates[0]
    flow = flow.resolve()
    arrays = (args.arrays or case/'frozen_flow/flow_arrays_si.npz').resolve()
    mesh_path = (args.mesh or case/'SV_MESH/mesh_arrays.npz').resolve()
    props = material(case, args.config)
    manifest_path = case/'frozen_flow/manifest.json'
    if manifest_path.exists() and flow.parent == manifest_path.parent.resolve():
        manifest = json.loads(manifest_path.read_text())
        if 'input_configuration_sha256' in manifest:
            assert props['config_sha256'] == manifest['input_configuration_sha256'], 'Configuration differs from frozen solve'
        if 'mesh_arrays_sha256' in manifest:
            assert sha(mesh_path) == manifest['mesh_arrays_sha256'], 'Mesh differs from frozen solve'
        if 'policy_sha256' in manifest:
            assert props['policy_sha256'] == manifest['policy_sha256'], 'Policy differs from frozen solve'
        for name, entry in manifest['files'].items():
            assert sha(manifest_path.parent/name) == entry['sha256'], name
    input_lock = case/'input_hashes.json'
    if input_lock.exists():
        locked = json.loads(input_lock.read_text())
        for key, path in [('run/solver.xml', Path(props['config'])),
                          ('policy.json', case/'policy.json'), ('SV_MESH/mesh_arrays.npz', mesh_path)]:
            assert key in locked and sha(path) == locked[key], 'Case input lock mismatch: '+key
    if context is not None:
        assert sha(flow) == context['source_field_sha256']
    frozen = np.load(arrays)
    points, tetra, triangles, tags, velocity, pressure = [frozen[k] for k in
        ['points_m', 'tetra', 'boundary_triangles', 'facet_tags', 'velocity_m_s', 'pressure_pa']]
    mesh_arrays = np.load(mesh_path)
    for key in ['points_m', 'tetra', 'boundary_triangles', 'facet_tags']:
        assert np.array_equal(frozen[key], mesh_arrays[key]), key
    grid = pv.read(flow)
    coordinate_record = coordinate_identity(grid.points, points)
    assert np.array_equal(np.sort(grid.cells_dict[pv.CellType.TETRA], axis=1), np.sort(tetra, axis=1))
    assert np.array_equal(grid['Velocity'], velocity) and np.array_equal(grid['Pressure'].reshape(-1), pressure)
    assert np.isfinite(velocity).all() and np.isfinite(pressure).all()
    xyz = points[triangles]
    oriented_area = .5*np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])
    facet_flux = np.einsum('ij,ij->i', oriented_area, velocity[triangles].mean(axis=1))
    flows = {int(tag): float(facet_flux[tags == tag].sum()) for tag in np.unique(tags)}
    target = float(json.loads((case/'policy.json').read_text())['Q_target_m3_s'])
    assert target > 0 and 4 in flows, 'SI target flow and inlet tag 4 are required'
    epsilon_Q = abs(-flows[4]-target)/target
    epsilon_mass = abs(sum(flows.values()))/target
    assert epsilon_Q <= 1e-6 and epsilon_mass <= 1e-6, (epsilon_Q, epsilon_mass)
    mesh, d = recover(points, tetra, triangles, tags, velocity, pressure, props['mu_Pa_s'])
    ids = d['ids']; magnitude = d['magnitude']; area = d['area']
    assert np.linalg.norm(velocity[ids], axis=1).max() <= 1e-12
    tangent_error = float(abs(np.einsum('ij,ij->i', d['traction'], d['normal'])).max())
    assert tangent_error < max(1., magnitude.max()) * 1e-12
    out = args.output.resolve(); (out/'data').mkdir(parents=True, exist_ok=True)
    mesh.save(out/'data/wall_wss_si.vtp')
    full, pressure_ids = wss.surface(points, triangles)
    full.point_data['Pressure_Pa'] = pressure[pressure_ids]; full.cell_data['Boundary_tag'] = tags
    full.save(out/'data/pressure_surface_si.vtp')
    record = dict(all_pass=True, case=case.name, source_case=str(case), source_field=str(flow),
        source_field_sha256=sha(flow), flow_arrays=str(arrays), flow_arrays_sha256=sha(arrays),
        mesh_arrays=str(mesh_path), mesh_arrays_sha256=sha(mesh_path), material=props,
        coordinate_identity=coordinate_record,
        viscosity_Pa_s=props['mu_Pa_s'], core_implementation=str(Path(wss.__file__).resolve()),
        core_sha256=sha(wss.__file__), entry_sha256=sha(__file__),
        case_adapter_sha256=sha(HERE.parent/'solver_support/src/flow_solver_support/wss_case.py'),
        pressure_Pa=dict(min=float(pressure.min()), max=float(pressure.max())),
        wss_Pa=dict(raw_min=float(magnitude.min()), raw_max=float(magnitude.max()),
            area_weighted_mean=float(np.average(magnitude, weights=area)),
            display_min=float(d['display'][ids].min()), display_max=float(d['display'][ids].max())),
        wss_method='Exact P1 tetra velocity gradients; tangential viscous traction on outward planar WALL facets',
        wss_display='Area-weighted mean of incident facet magnitudes at original wall nodes; geometry unchanged',
        wss_vector_convention='(I-n*n^T)*mu*(grad(u)+grad(u)^T)*n',
        wall_facets=len(magnitude), tangency_max_error_Pa=tangent_error,
        max_wall_speed_m_s=float(np.linalg.norm(velocity[ids], axis=1).max()),
        source_values_clipped=False, mesh_convergence_not_assessed=True,
        measurements=dict(outward_flows_by_tag_m3_s=flows, epsilon_Q=epsilon_Q, epsilon_mass=epsilon_mass),
        outputs_sha256={str(p.relative_to(out)):sha(p) for p in (out/'data').glob('*.vtp')},
        elapsed_seconds=time.perf_counter()-started)
    (out/'COMPUTE_VALIDATION.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps(record,indent=2),flush=True)


if __name__ == '__main__':
    main()
