"""Fresh-read frozen FEM artifacts and validate their common origin; no CFD.

All operational paths resolve within this checkout. This is an integrity reader,
not a particle locator, interpolator, or flow solver.
"""
import argparse
import hashlib
import json
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from sv_validation.sv13n import checkpoint_one_rank
from sv_validation.postprocess import SolutionMeasurements


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((ROOT / name).read_text())


def file_record(path):
    return dict(path=path.relative_to(ROOT).as_posix(), size_bytes=path.stat().st_size, sha256=sha(path))


def arrays(attributes):
    return [dict(name=k, dtype=str(v.dtype), components=1 if v.ndim == 1 else v.shape[1],
                 tuples=len(v)) for k, v in attributes.items()]


def canonical_cells(cells):
    a = np.sort(cells.reshape(-1, 5)[:, 1:], axis=1)
    return a[np.lexsort(a[:, ::-1].T)]


def validate():
    q = 'reports/sv1_3q/'
    accepted = read(q + 'REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance.json')
    winner = read(q + 'winner_steady_candidate.json')
    execution = read(q + 'remote/REAL_VASCULAR_GPU_ILU_REUSE_WINNER_execution.json')
    steady = read(q + 'steady_history.json')
    stop = read(q + 'stop_latency.json')
    policy = read('configs/sv1_3q/policy.json')
    status = read(q + 'stage_status.json')
    assert all(d['status'] == 'PASS' for d in [accepted, winner, steady, stop])
    assert status['completed'] and status['status'] == 'ILU_REBUILD_POLICY_WINNER_FOUND'
    assert execution['exit_code'] == 0
    base = ROOT / 'frozen_reference'
    flow_path = base / 'flow/steady_flow_stage_sv1_3q.vtu'
    mesh_path = base / 'SV_MESH/mesh-complete.mesh.vtu'
    xml_path = base / 'run/solver.xml'
    cp_path = base / 'fem_checkpoint/stFile_071.bin'
    flow, mesh = pv.read(flow_path), pv.read(mesh_path)
    assert sha(flow_path) == accepted['reload']['sha256'] == winner['VTU_sha256']
    assert any(f['sha256'] == sha(flow_path) and f['bytes'] == flow_path.stat().st_size for f in execution['results'])
    assert sha(xml_path) == accepted['config_sha256'] == execution['config_sha256'] == winner['config_sha256']
    assert np.all(mesh.celltypes == 10) and np.all(flow.celltypes == 10)
    assert mesh.n_points == flow.n_points and mesh.n_cells == flow.n_cells
    assert np.array_equal(mesh.points, flow.points)
    assert np.array_equal(canonical_cells(mesh.cells), canonical_cells(flow.cells))
    # The actual export swaps local vertex order; tetra cell rows are unchanged.
    mesh_cells = mesh.cells.reshape(-1, 5)[:, 1:]
    flow_cells = flow.cells.reshape(-1, 5)[:, 1:]
    rows_aligned = bool(np.array_equal(np.sort(mesh_cells, axis=1), np.sort(flow_cells, axis=1)))
    import itertools
    permutation = next((list(p) for p in itertools.permutations(range(4))
                        if np.array_equal(mesh_cells[:, p], flow_cells)), None)
    assert np.unique(mesh.cells.reshape(-1, 5)[:, 0]).tolist() == [4]
    assert flow.point_data['Velocity'].shape == (flow.n_points, 3)
    assert flow.point_data['Pressure'].shape == (flow.n_points,)
    assert set(flow.point_data.keys()) == {'Velocity', 'Pressure'}
    assert not flow.cell_data.keys()
    velocity, pressure = flow['Velocity'], flow['Pressure']
    assert np.isfinite(velocity).all() and np.isfinite(pressure).all()
    time = float(flow.field_data['TimeValue'][0])
    xml = ET.parse(xml_path).getroot()
    dt = float(xml.findtext('.//Time_step_size'))
    step = accepted['stop_step']
    assert step == winner['stop_step'] == stop['actual_stop_step'] == steady['stop_step']
    assert math.isclose(time, step * dt, rel_tol=1e-12)
    checkpoint = checkpoint_one_rank(cp_path, step, dt)
    assert checkpoint['sha256'] == accepted['checkpoint']['sha256']
    assert checkpoint['nodes'] == flow.n_points and math.isclose(checkpoint['time_s'], time, rel_tol=1e-12)
    checkpoint['path'] = cp_path.relative_to(ROOT).as_posix()
    checkpoint['source_path'] = accepted['checkpoint']['path']
    checkpoint['solver_sha256'] = accepted['solver_sha256']
    checkpoint['role'] = 'FEM_REPRODUCIBILITY_EVIDENCE_NOT_REQUIRED_FOR_PARTICLE_SAMPLING'
    assert steady['first_full_steady_step'] == winner['first_full_steady_step'] == stop['first_qualifying_saved_step']
    assert steady['states'][-1]['step'] == winner['first_full_steady_step']
    assert sha(ROOT / policy['production_policy_path']) == steady['policy_sha256']
    assert all(d['E_u'] <= policy['production_policy']['velocity_change_limit'] and
               d['E_Q'] <= policy['production_policy']['flow_change_limit'] for d in steady['intervals'][-5:])
    # Exact files named in the actual run XML must resolve without the WSL source tree.
    for node in xml.findall('.//Mesh_file_path') + xml.findall('.//Face_file_path'):
        assert (xml_path.parent / node.text).resolve().is_relative_to(ROOT)
        assert (xml_path.parent / node.text).is_file()
    integrity = read(q + 'remote/final_input_integrity.json')
    for f in integrity['files']:
        relative = f['path'].split('/outputs/SV_MESH/')[1]
        assert sha(base / 'SV_MESH' / relative) == f['sha256']
    import yaml
    units = yaml.safe_load((ROOT / 'configs/sv_reference.yaml').read_text())
    assert units['units']['length'] == read('reports/sv1/geometry_import.json')['units'] == 'm'
    physics = units['physics']
    assert float(xml.findtext('.//Density')) == physics['density_kg_m3']
    assert float(xml.findtext('.//Viscosity/Value')) == physics['dynamic_viscosity_pa_s']
    assert sha(base/'SV_MESH/mesh_arrays.npz') == winner['mesh_sha256']
    measurement = SolutionMeasurements(base/'SV_MESH/mesh_arrays.npz', physics['inlet_volume_flow_m3_s'],
                                       policy['production_policy']['Umean_m_s'])
    u, p = measurement.read(flow_path)
    measured = measurement.measure(u, p)
    for name in ['Q_in_m3_s', 'Q_out_total_m3_s', 'velocity_max_m_s', 'wall_velocity_max_m_s', 'epsilon_mass']:
        assert math.isclose(measured[name], accepted['measurement'][name], rel_tol=1e-12, abs_tol=1e-30), name
    assert measured['wall_noslip_pass'] and measured['epsilon_mass'] <= policy['production_policy']['mass_limit']
    fmap = read('configs/face_map.json')
    tetra = mesh.cells.reshape(-1, 5)[:, 1:]
    assert np.array_equal(mesh['GlobalNodeID'], np.arange(1, mesh.n_points+1))
    assert np.array_equal(mesh['GlobalElementID'], np.arange(1, mesh.n_cells+1))
    exterior = pv.read(base/'SV_MESH/mesh-complete.exterior.vtp')
    boundaries = {}
    boundary_keys = []
    for f in fmap['faces']:
        role, face_id = f['role'], f['sv_face_id']
        path = base/'SV_MESH/mesh-surfaces'/f'{role}.vtp'
        surface = pv.read(path)
        faces = surface.faces.reshape(-1, 4)
        assert np.all(faces[:, 0] == 3)
        tri = faces[:, 1:]
        xyz = surface.points[tri]
        global_nodes = surface['GlobalNodeID'] - 1
        assert np.array_equal(surface.points, mesh.points[global_nodes])
        elem = surface['GlobalElementID'] - 1
        global_tri = global_nodes[tri]
        assert np.all(np.any(global_tri[:, :, None] == tetra[elem, None, :], axis=2))
        area_vector = .5 * np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
        areas = np.linalg.norm(area_vector, axis=1)
        orientation = np.einsum('ij,ij->i', area_vector, xyz.mean(1) - mesh.points[tetra[elem]].mean(1))
        assert np.all(areas > 0) and np.all(orientation > 0), role
        boundary_keys.extend(map(tuple, np.sort(global_tri, axis=1).tolist()))
        assert np.count_nonzero(exterior['ModelFaceID'] == face_id) == surface.n_cells
        integrated = area_vector.sum(0)
        signed_flux = float(np.sum(np.einsum('ij,ij->i', velocity[global_tri].mean(1), area_vector)))
        if role != 'WALL':
            assert math.isclose(signed_flux, measured['signed_outward_boundary_flows_m3_s'][role], rel_tol=1e-12)
        boundaries[role] = dict(**file_record(path), sv_face_id=face_id, units='m', area_m2=float(areas.sum()),
            nodes=surface.n_points, facets=surface.n_cells, normal_array_stored=False,
            raw_triangle_winding='OUTWARD_FROM_FLUID_LUMEN', outward_facets=int(np.count_nonzero(orientation > 0)),
            normal_definition='n_out = cross(x1-x0,x2-x0)/norm(cross); verified against owning tetra centroid',
            inward_normal='-n_out', signed_flux_m3_s=signed_flux,
            integrated_area_vector_m2=integrated.tolist(),
            area_weighted_mean_unit_normal=(integrated/areas.sum()).tolist(),
            note='Use per-triangle outward normals; mean is not a replacement for local wall normals.')
    assert len(set(boundary_keys)) == len(boundary_keys) == exterior.n_cells
    ex_keys = set(map(tuple, np.sort((exterior['GlobalNodeID']-1)[exterior.faces.reshape(-1,4)[:,1:]], axis=1).tolist()))
    assert ex_keys == set(boundary_keys)
    boundary = dict(status='PASS', coordinates='Same unshifted right-handed Cartesian x/y/z as volume/flow; no anatomical axis labels assigned',
                    units='m', flux_sign='integral(u dot n_out dA): INLET negative; OUTLET positive; Qin = -inlet_signed_flux',
                    source_face_map='configs/face_map.json', boundaries=boundaries,
                    exterior=file_record(base/'SV_MESH/mesh-complete.exterior.vtp'), complete_disjoint_surface_partition=True)
    mesh_manifest = dict(status='PASS', **file_record(mesh_path), canonical_role='FROZEN_PARTICLE_VOLUME_MESH',
        nodes=mesh.n_points, tetra=mesh.n_cells, vtk_cell_type=10, nodes_per_cell=4,
        surface_facets=exterior.n_cells, units='m', bounding_box_m=list(mesh.bounds),
        boundary_manifest='frozen_reference/boundary_manifest.json',
        face_ids={r:v['sv_face_id'] for r,v in boundaries.items()},
        arrays=dict(point=arrays(mesh.point_data), cell=arrays(mesh.cell_data)),
        source='outputs/sv1/SV_MESH/mesh-complete.mesh.vtu',
        actual_stage_q_input_verified=True, solver_input_hash_source=q+'remote/final_input_integrity.json',
        coordinate_system=boundary['coordinates'],
        tetra_id_convention='zero-based cell index in canonical volume mesh; its GlobalElementID is one-based',
        flow_points_identical=True, flow_tetra_sets_identical=True,
        flow_cell_connectivity_order_identical=bool(np.array_equal(mesh.cells, flow.cells)),
        flow_cell_rows_identical_ignoring_local_vertex_order=rows_aligned,
        flow_local_vertex_permutation_from_volume=permutation,
        ordering_note='This frozen pair has identical tetra cell rows and point order. Only local vertex order differs; use canonical volume cell index for tetra_id and revalidate correspondence for any future file.')
    flow_manifest = dict(status='PASS', **file_record(flow_path), canonical_role='FROZEN_PARTICLE_BACKGROUND_FLOW',
        source=accepted['reload']['path'], run=accepted['name'], stage='SV1.3Q',
        points=flow.n_points, cells=flow.n_cells, TimeValue=time, time_units='s', step=step,
        point_arrays=arrays(flow.point_data), cell_arrays=arrays(flow.cell_data), field_arrays=arrays(flow.field_data),
        velocity=dict(name='Velocity', association='POINT', components=3, dtype=str(velocity.dtype), units='m/s', finite=True),
        pressure=dict(name='Pressure', association='POINT', components=1, dtype=str(pressure.dtype), units='Pa', finite=True),
        velocity_representation='Exported nodal values on four-node linear tetrahedra; no cell velocity array. This records the export, not an assumption about all internal solver spaces.',
        velocity_gradient='NOT STORED IN FROZEN VTU', vorticity='NOT STORED IN FROZEN VTU', strain_rate='NOT STORED IN FROZEN VTU',
        derived_field_action='Particle-0 must derive gradients from tetra nodal velocity and define interpolation and boundary behavior; no such implementation is included.',
        max_velocity_m_s=float(np.linalg.norm(velocity, axis=1).max()),
        pressure_range_pa=[float(pressure.min()), float(pressure.max())], coordinate_units='m',
        coordinate_system=boundary['coordinates'], mesh='frozen_reference/SV_MESH/mesh-complete.mesh.vtu',
        mesh_sha256=sha(mesh_path), solver_config=file_record(xml_path), solver_sha256=accepted['solver_sha256'],
        first_formal_steady_step=winner['first_full_steady_step'], safe_stop_step=step,
        checkpoint_time_s=checkpoint['time_s'], monitor_last_saved_step=steady['states'][-1]['step'])
    build = read(q+'remote/svmp_reuse_build.json')
    mpi = read('reports/sv1_3l/mpi_fortran_build.json')
    environment = read('frozen_reference/run/gpu_environment_manifest.json')
    gpu_text = '\n'.join(p.get('stdout','') for p in environment['probes'] if p['command'] == ['nvidia-smi'])
    import re
    gpu_match = re.search(r'NVIDIA GeForce RTX \d+', gpu_text)
    assert gpu_match is not None
    runtime = dict(stage='SV1.3Q', status=status['status'], solver='svMultiPhysics fluid + PETSc GMRES/right + ASM overlap 2 + ILU(2)',
        GPU=gpu_match[0], GPU_evidence='frozen_reference/run/gpu_environment_manifest.json',
        upstream=read('sync_metadata/upstream_source.json'), solver_sha256=accepted['solver_sha256'],
        PETSc_version=build['PETSc_build']['version'], PETSc_commit=build['PETSc_build']['commit'],
        PETSc_library_sha256=accepted['PETSc_library_sha256'], CUDA_version=build['PETSc_build']['CUDA_version'],
        MPI_version=mpi['version'], MPI_bindings=mpi['bindings_requested'],
        MPI_ranks=execution['MPI_ranks'], GPUs=execution['GPUs'], OMP_NUM_THREADS=execution['OMP_NUM_THREADS'],
        command=execution['command'], PETSC_OPTIONS=execution['PETSC_OPTIONS'],
        runtime_semantics=accepted['runtime_semantics'][0], matrix_type=accepted['Mat'], vector_type=accepted['Vec'],
        ILU_policy=winner['winner'], adaptive_threshold=policy['adaptive_threshold'], max_reuse_age=policy['max_reuse_age'],
        reference_rule=policy['I_ref_rule'], actual_rebuild_reasons=sorted({t['rebuild_reason'] for t in accepted['reuse']['trace'] if not t['reuse']}),
        raw_execution=q+'remote/REAL_VASCULAR_GPU_ILU_REUSE_WINNER_execution.json')
    summary = dict(status='PASS', stage='SV1.3Q', stage_status=status['status'],
        physics=physics, dt_s=dt, time_integrator='generalized-alpha', spectral_radius=float(xml.findtext('.//Spectral_radius_of_infinite_time_step')),
        BCs={bc.attrib['name']:{child.tag:child.text for child in bc} for bc in xml.findall('.//Add_BC')},
        production_policy=policy['production_policy'], measurement=measured, statistics=accepted['statistics'],
        ILU_rebuild_count=accepted['reuse']['ILU_rebuild_count'], ILU_reuse_count=accepted['reuse']['ILU_reuse_count'],
        linear_unrecovered_failures=accepted['linear_failures'], nonlinear_failures=accepted['nonlinear_failures'],
        recovery_count=accepted['reuse']['recovery_count'], wall_time_s=accepted['wall_time_s'],
        first_formal_steady_step=winner['first_full_steady_step'], safe_stop_step=step,
        all_artifacts_from_same_run=True, no_CFD_executed=True,
        deferred=dict(mesh_convergence='NOT PERFORMED BY USER DECISION / NOT PLANNED BEFORE PARTICLE DEVELOPMENT',
                      time_step_sensitivity='NOT PERFORMED BY USER DECISION / NOT PLANNED BEFORE PARTICLE DEVELOPMENT',
                      CPU_GPU_field_equivalence='DEFERRED BY USER DECISION'),
        interpretation='Stage Q is the frozen particle-development background-flow baseline; not new FEM production approval.')
    return {'frozen_reference/mesh_manifest.json':mesh_manifest,
            'frozen_reference/boundary_manifest.json':boundary,
            'frozen_reference/flow/flow_field_manifest.json':flow_manifest,
            'frozen_reference/fem_checkpoint/checkpoint_manifest.json':checkpoint,
            'frozen_reference/run/runtime_manifest.json':runtime,
            'frozen_reference/baseline_summary.json':summary}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Write regenerated manifests; default compares existing manifests')
    args = parser.parse_args()
    for name, value in validate().items():
        path = ROOT / name
        if args.write:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
        else:
            assert json.loads(path.read_text()) == value, name
    print('PASS: fresh process loaded flow, volume, all boundaries and checkpoint; same-run hashes/config/time/flux/normal checks passed. No CFD.')


if __name__ == '__main__':
    main()
