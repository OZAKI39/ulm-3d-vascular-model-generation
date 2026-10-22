"""Read-only plane integrals in the verified straight artificial port tubes."""
from pathlib import Path
import csv,hashlib,json
import numpy as np,pyvista as pv
ROOT=Path(__file__).resolve().parents[1]
active=json.loads((ROOT/'reports/ACTIVE_FLOW.json').read_text());assert active['status']=='PASS'
case=ROOT/active['case'];g=json.loads((ROOT/'inputs/geometry.json').read_text())
path=case/'frozen_flow/steady_flow.vtu';assert hashlib.sha256(path.read_bytes()).hexdigest()==active['flow_sha256']
grid=pv.read(ROOT/'mesh/SV_MESH/mesh-complete.mesh.vtu');grid.points=np.asarray(grid.points,dtype=float)
flow=pv.read(path);assert np.array_equal(flow.points,grid.points)
grid.point_data['Velocity']=flow.point_data['Velocity']
rows=[];out=ROOT/'reports/internal_sections';out.mkdir(exist_ok=True)
for name,spec in g['ports'].items():
    role=spec['role'];surface=pv.read(ROOT/'mesh/SV_MESH/mesh-surfaces'/f'{role}.vtp')
    center=np.array(surface.center);normal=np.array(spec['outward_normal']);normal/=np.linalg.norm(normal)
    reference=active['measurements']['signed_outward_boundary_flows_m3_s'][role]
    for depth_mm in [2.,5.,10.]:
        expected=center-depth_mm*.001*normal
        cut=grid.slice(normal=normal,origin=expected,generate_triangles=True).connectivity()
        candidates=[]
        for region in np.unique(cut.cell_data['RegionId']):
            sub=cut.extract_cells(cut.cell_data['RegionId']==region).extract_surface(algorithm='dataset_surface').triangulate()
            ids=sub.faces.reshape(-1,4)[:,1:];xyz=np.asarray(sub.points)[ids]
            area=np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)*.5
            centroid=np.average(xyz.mean(1),axis=0,weights=area)
            candidates.append((float(np.linalg.norm(centroid-expected)),sub,ids,area,centroid))
        distance,sub,ids,area,centroid=min(candidates,key=lambda item:item[0])
        # Reject a cut that follows a curve or accidentally selects another branch.
        geometry_valid=distance<.3*spec['radius_m']
        q=float(np.sum(area*(sub.point_data['Velocity'][ids].mean(1)@normal)))
        filename=f'{role}_inward_{depth_mm:g}mm.vtp';sub.save(out/filename)
        rows.append(dict(boundary=role,inward_offset_mm=depth_mm,geometry_valid=geometry_valid,
            center_offset_from_straight_axis_m=distance,x_m=float(centroid[0]),y_m=float(centroid[1]),z_m=float(centroid[2]),
            area_m2=float(area.sum()),signed_outward_Q_m3_s=q,flow_magnitude_uL_min=abs(q)*6e10,
            signed_boundary_reference_m3_s=reference,relative_deviation_from_port=abs(q-reference)/abs(reference),
            deviation_as_fraction_total_inlet=abs(q-reference)/active['measurements']['Q_in_m3_s'],file=filename))
with (out/'section_flux.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
valid=[row for row in rows if row['geometry_valid']]
record=dict(source_flow_sha256=active['flow_sha256'],method='Exact linear-triangle velocity integral on tetrahedral plane cuts; connected component closest to the straight artificial-port axis',
    coordinate_unit='m',velocity_unit='m/s',rows=rows,valid_sections=len(valid),
    maximum_deviation_fraction_total_inlet=max(r['deviation_as_fraction_total_inlet'] for r in valid),
    scope='Artificial straight port tube sections only; not a proof of local solenoidality, mesh independence, or physical trapping')
(out/'INTERNAL_SECTION_CHECK.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
