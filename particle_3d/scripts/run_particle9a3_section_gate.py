"""Search the frozen pre-junction root on CPU; stop before any births if invalid."""
from pathlib import Path
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import time
import csv
import platform

import numpy as np
import pyvista as pv

from particle_3d.interior_section import (root_topology, candidate_sequence, evaluate_candidate,
    locator_for, ALGORITHM_VERSION, cut_tetrahedra, root_component)
from particle_3d.audit import read_frozen

REPO = Path(__file__).resolve().parents[2]
FROZEN = REPO/'formal_3D_flow_solver/FEM_SimVascular/frozen_reference'
REPORT = REPO/'particle_3d/reports/particle9a3_interior_inlet'
OUTPUT = REPO/'particle_3d/outputs/particle9a3_inlet_audit'
ENV = None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def encode(obj):
    if isinstance(obj, np.ndarray): return obj.tolist()
    if isinstance(obj, np.generic): return obj.item()
    raise TypeError(type(obj).__name__)


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, default=encode, indent=2, allow_nan=False)+'\n')


def environment():
    global ENV
    if ENV is not None: return ENV
    _, grid, flow, boundaries = read_frozen(FROZEN.parent)
    grid.points = np.asarray(grid.points, dtype=np.float64)
    # Exactly the coordinate convention used by FrozenFEMField.from_grids.
    grid.point_data['Velocity'] = np.asarray(flow['Velocity'], float)
    grid.cell_data['canonical_tetra_id'] = np.arange(grid.n_cells)
    surfaces = {}
    for role in ['WALL', 'INLET', 'OUTLET_01', 'OUTLET_02', 'OUTLET_03']:
        s = boundaries[role].copy()
        s.points = grid.points[np.asarray(s.point_data['GlobalNodeID'], int)-1]
        surfaces[role] = s
    inlet = surfaces['INLET']
    ids = inlet.faces.reshape(-1, 4)[:, 1:]
    xyz = inlet.points[ids]
    vector_area = np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])/2
    area = np.linalg.norm(vector_area, axis=1)
    center = np.sum(area[:, None]*xyz.mean(axis=1), axis=0)/area.sum()
    swc = np.loadtxt(REPORT/'reference/roi_core.swc')
    topo = root_topology(swc, center)
    normals = vector_area/area[:, None]
    inward = topo['points_m'][1]-center
    normals[normals@inward < 0] *= -1
    vi = grid['Velocity'][np.asarray(inlet.point_data['GlobalNodeID'], int)-1][ids]
    q = float(np.sum(area*np.einsum('tvi,ti->tv', vi, normals).mean(axis=1)))
    q_saved = json.loads((FROZEN/'new_flow_flux.json').read_text())['integrated_inlet_Q_m3_s']
    if not np.isclose(q, q_saved, atol=0, rtol=1e-12):
        raise ValueError('Independent inlet integration does not match saved reference')
    policy = json.loads((FROZEN/'source_evidence/policy.json').read_text())
    ENV = (grid, topo, locator_for(surfaces['WALL']),
           {k: locator_for(s) for k, s in surfaces.items() if k != 'WALL'}, q, policy, surfaces)
    return ENV


def job(candidate):
    grid, topo, wall, caps, q, policy, _ = environment()
    data, _ = evaluate_candidate(grid, candidate, topo, wall, caps, q, policy['mass_limit'])
    return data


def independent_flux(grid, section, center, normal):
    """Independent numpy edge intersection/fan integration, bypassing VTK cutting."""
    tetra = grid.cells.reshape(-1, 5)[:, 1:]
    ids = np.unique(section.cell_data['canonical_tetra_id']).astype(int)
    basis = np.eye(3)[np.argmin(np.abs(normal))]
    e1 = np.cross(normal, basis); e1 /= np.linalg.norm(e1)
    e2 = np.cross(normal, e1)
    flux, area_sum = 0., 0.
    for tid in ids:
        x = grid.points[tetra[tid]]; v = grid['Velocity'][tetra[tid]]
        distances = (x-center)@normal
        polygon, velocity = [], []
        for i in range(4):
            for j in range(i+1, 4):
                if distances[i]*distances[j] < 0:
                    f = distances[i]/(distances[i]-distances[j])
                    polygon.append(x[i]+f*(x[j]-x[i])); velocity.append(v[i]+f*(v[j]-v[i]))
        if len(polygon) < 3: continue
        polygon, velocity = np.array(polygon), np.array(velocity)
        offset = polygon-polygon.mean(axis=0)
        order = np.argsort(np.arctan2(offset@e2, offset@e1))
        polygon, velocity = polygon[order], velocity[order]
        for j in range(1, len(polygon)-1):
            t = [0, j, j+1]
            x3 = polygon[t]
            area = np.linalg.norm(np.cross(x3[1]-x3[0], x3[2]-x3[0]))/2
            area_sum += area
            flux += area*float((velocity[t]@normal).mean())
    return dict(independent_signed_Q_m3_s=flux, independent_area_m2=area_sum, tetra_count=len(ids))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    start = time.monotonic()
    grid, topo, wall, caps, q, policy, surfaces = environment()
    print('Canonical mesh loaded; beginning deterministic root search.', flush=True)
    h = policy['h10_m'] if 'h10_m' in policy else policy['h10_mesh_m']
    candidates = list(candidate_sequence(topo, h, 1))
    # Both schedules are fixed before examining velocity/particle results.
    refined = list(candidate_sequence(topo, h, 2))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(job, candidates))
        print('Primary candidates:', len(rows), 'passes:', sum(r['accepted'] for r in rows), flush=True)
        fine_rows = list(pool.map(job, refined))
    # Selection belongs to the primary schedule. Fine schedule is only a resolution audit.
    accepted = next((r for r in rows if r['accepted']), None)
    all_rows = rows+fine_rows
    for row in all_rows:
        for field in ['center_m', 'normal']:
            for j, axis in enumerate('xyz'):
                row[f'{field}_{axis}'] = float(row[field][j])
    (REPORT/'data').mkdir(parents=True, exist_ok=True)
    flat = [{k: v for k, v in row.items() if k not in ('center_m', 'normal')} for row in all_rows]
    with (REPORT/'data/interior_section_search.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0]))
        writer.writeheader(); writer.writerows(flat)
    audit_rows = []
    witnesses = [rows[0], min(rows, key=lambda r: r['flux_relative_error']),
                 rows[len(rows)//2], rows[-1], fine_rows[0]]
    for index, row in enumerate(witnesses):
        _, section = evaluate_candidate(grid, row, topo, wall, caps, q, policy['mass_limit'])
        independent = independent_flux(grid, section, row['center_m'], row['normal'])
        independent.update(candidate_id=row['candidate_id'], refinement=row['refinement'],
                           vtk_signed_Q_m3_s=row['signed_Q_m3_s'],
                           relative_method_difference=abs(independent['independent_signed_Q_m3_s']-row['signed_Q_m3_s'])/q)
        audit_rows.append(independent)
        section.save(REPORT/f'data/rejected_section_witness_{index:02d}.vtp')
    save_json(REPORT/'data/independent_section_flux_validation.json', audit_rows)
    save_json(REPORT/'data/root_topology.json', topo)
    sources = {str(p.relative_to(REPO)): sha(p) for p in [
        FROZEN/'flow/steady_flow_mean_2p0_mmps.vtu', FROZEN/'SV_MESH/mesh-complete.mesh.vtu',
        FROZEN/'SV_MESH/mesh-surfaces/WALL.vtp', FROZEN/'source_evidence/policy.json',
        REPORT/'reference/roi_core.swc', REPORT/'reference/source_contract.json']}
    manifest = dict(status='INTERIOR_SECTION_FOUND' if accepted else 'NO_VALID_INTERIOR_INJECTION_SECTION',
        human_review='PENDING', algorithm_version=ALGORITHM_VERSION, selected_section=accepted,
        source_sha256=sources, centerline_source='reference/roi_core.swc linked by reference/source_contract.json lineage.canonical_roi_swc',
        extension_source='Actual cap area centroid to canonical ROI root; source port_classification.csv documents the artificial extension',
        axis_fallback=False, root_length_to_first_junction_m=topo['arclength_m'][-1],
        junction_node_id=topo['junction_node_id'], inlet_center_m=topo['points_m'][0],
        first_junction_m=topo['junction_m'], reference_inlet_Q_m3_s=q,
        flux_relative_tolerance=policy['mass_limit'], tolerance_source='frozen_reference/source_evidence/policy.json:mass_limit',
        candidate_spacing_m=h, primary_candidates=len(rows), refined_candidates=len(fine_rows),
        primary_pass_count=sum(r['accepted'] for r in rows), refined_pass_count=sum(r['accepted'] for r in fine_rows),
        primary_geometry_pass_count=sum(r['geometry_pass'] for r in rows),
        rejection_counts=dict(Counter(reason for r in rows for reason in r['rejected_reason'].split(';') if reason)),
        best_flux_error_primary=min(r['flux_relative_error'] for r in rows),
        best_flux_error_refined=min(r['flux_relative_error'] for r in fine_rows),
        maximum_independent_method_difference=max(r['relative_method_difference'] for r in audit_rows),
        planned_birth_audit_count=2000, executed_birth_audit_count=0, executed_point_tracers=0,
        planned_smoke_count=30, executed_smoke_count=0, executed_formal_count=0,
        downstream_status='NOT_EXECUTED_SECTION_GATE_FAILED' if not accepted else 'AWAITING_JOINT_IMPLEMENTATION',
        workers=args.workers, runtime_s=time.monotonic()-start, python=platform.python_version(),
        host=platform.node(), execution_root=str(REPO), utc=datetime.now(timezone.utc).isoformat())
    if accepted:
        _, s = evaluate_candidate(grid, accepted, topo, wall, caps, q, policy['mass_limit'])
        s.save(REPORT/'data/interior_injection_section.vtp')
    save_json(REPORT/'data/interior_injection_section.json', manifest)
    save_json(OUTPUT/'section_gate.json', manifest)
    print(json.dumps(manifest, default=encode, indent=2), flush=True)


if __name__ == '__main__': main()
