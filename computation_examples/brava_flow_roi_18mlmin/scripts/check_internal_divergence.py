"""Read-only divergence theorem check of the actual P1 field in port slabs."""
from pathlib import Path
import csv, hashlib, json, sys
import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'vendor'))
from flow_solver_support.wss import p1_gradients

active = json.loads((ROOT/'reports/ACTIVE_FLOW.json').read_text())
geometry = json.loads((ROOT/'inputs/geometry.json').read_text())
flow_path = ROOT/active['case']/'frozen_flow/steady_flow.vtu'
assert hashlib.sha256(flow_path.read_bytes()).hexdigest() == active['flow_sha256']
grid = pv.read(ROOT/'mesh/SV_MESH/mesh-complete.mesh.vtu')
grid.points = np.asarray(grid.points, dtype=np.float64)
flow = pv.read(flow_path)
assert np.array_equal(flow.points, grid.points)
tetra = grid.cells.reshape(-1,5)[:,1:]
grad = p1_gradients(grid.points, tetra, np.asarray(flow.point_data['Velocity']))
grid.cell_data['divergence_per_s'] = np.trace(grad,axis1=1,axis2=2)
out = ROOT/'reports/internal_sections'
section_rows = list(csv.DictReader((out/'section_flux.csv').open()))
rows = []
for port in geometry['ports'].values():
    role = port['role']
    normal = np.asarray(port['outward_normal'],dtype=float)
    normal /= np.linalg.norm(normal)
    cap = pv.read(ROOT/'mesh/SV_MESH/mesh-surfaces'/f'{role}.vtp')
    center = np.asarray(cap.center)
    depth = .002
    cut_center = center-depth*normal
    # Keep the side pointing toward the external port. Other branches may be
    # disconnected fragments and are excluded by closest connected component.
    clipped = grid.clip(normal=normal, origin=cut_center, invert=False).connectivity()
    candidates = []
    for label in np.unique(clipped.cell_data['RegionId']):
        sub = clipped.extract_cells(clipped.cell_data['RegionId']==label)
        volume = sub.compute_cell_sizes(length=False,area=False,volume=True).cell_data['Volume']
        centroid = np.average(sub.cell_centers().points,axis=0,weights=volume)
        candidates.append((float(np.linalg.norm(centroid-(center-.5*depth*normal))), sub, volume, centroid))
    offset, slab, volume, centroid = min(candidates,key=lambda x:x[0])
    assert offset < port['radius_m']
    signed_distance = (slab.points-cut_center)@normal
    assert signed_distance.min()>-1e-8 and signed_distance.max()<depth+1e-8
    divergence_integral = float(np.sum(volume*slab.cell_data['divergence_per_s']))
    section = next(r for r in section_rows if r['boundary']==role and float(r['inward_offset_mm'])==2.)
    qcap = active['measurements']['signed_outward_boundary_flows_m3_s'][role]
    qsection = float(section['signed_outward_Q_m3_s'])
    flux_difference = qcap-qsection
    closure_error = abs(divergence_integral-flux_difference)
    rows.append(dict(boundary=role,depth_mm=2.,volume_m3=float(volume.sum()),
        signed_cap_flow_m3_s=qcap,signed_section_flow_m3_s=qsection,
        surface_flux_difference_m3_s=flux_difference,volume_divergence_integral_m3_s=divergence_integral,
        closure_error_m3_s=closure_error,closure_error_fraction_inlet=closure_error/active['measurements']['Q_in_m3_s'],
        divergence_defect_fraction_inlet=abs(divergence_integral)/active['measurements']['Q_in_m3_s'],
        clipped_cells=slab.n_cells))
with (out/'divergence_theorem.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
record = dict(source_flow_sha256=active['flow_sha256'], rows=rows,
    method='Exact constant P1 cell divergence integrated over clipped cap-side tube slabs; independent of triangle surface flux summation',
    maximum_closure_error_fraction_inlet=max(r['closure_error_fraction_inlet'] for r in rows),
    interpretation='A nonzero slab divergence integral confirms a local resolved-velocity conservation defect; this does not alter the solved field.')
(out/'DIVERGENCE_THEOREM_CHECK.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
