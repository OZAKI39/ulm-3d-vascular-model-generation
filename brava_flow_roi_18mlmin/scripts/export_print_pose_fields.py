"""Export the same accepted physical solution in each rigid print-pose SI frame."""
from pathlib import Path
import hashlib,json
import numpy as np,pyvista as pv
ROOT=Path(__file__).resolve().parents[1]
active=json.loads((ROOT/'reports/ACTIVE_FLOW.json').read_text());assert active['status']=='PASS'
case=ROOT/active['case'];g=json.loads((ROOT/'inputs/geometry.json').read_text())
original=case/'frozen_flow/steady_flow.vtu';assert hashlib.sha256(original.read_bytes()).hexdigest()==active['flow_sha256']
sources={'mesh':pv.read(ROOT/'mesh/SV_MESH/mesh-complete.mesh.vtu'),'flow':pv.read(original),'raw_wall_wss':pv.read(case/'postprocess/wall_wss_si.vtp'),'pressure_surface':pv.read(case/'postprocess/pressure_surface_si.vtp')}
for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03','WALL']:
    sources['mesh-surfaces/'+role]=pv.read(ROOT/'mesh/SV_MESH/mesh-surfaces'/f'{role}.vtp')
for pose in ['0','15']:
    T=np.array(g['pose_transforms'][pose]['source_mm_to_candidate_mm']);R=T[:3,:3]
    assert np.allclose(R.T@R,np.eye(3),atol=1e-14) and np.isclose(np.linalg.det(R),1)
    out=case/'frozen_flow'/('candidate_'+pose);out.mkdir(exist_ok=True)
    records={}
    for name,source in sources.items():
        transformed=source.copy(deep=True)
        transformed.points=((np.asarray(source.points,dtype=float)*1000+g['origin_source_mm'])@R.T+T[:3,3])*.001
        for source_values,target_values in [(source.point_data,transformed.point_data),(source.cell_data,transformed.cell_data)]:
            for field in ['Velocity','Tangential_viscous_traction_Pa','Outward_normal']:
                if field in source_values:
                    target_values[field]=source_values[field]@R.T
                    assert np.allclose(np.linalg.norm(target_values[field],axis=1),np.linalg.norm(source_values[field],axis=1),rtol=1e-13,atol=1e-14)
            for field in ['Pressure','Pressure_Pa','WSS_raw_Pa','WSS_display_Pa','Area_m2']:
                if field in source_values:assert np.array_equal(source_values[field],target_values[field])
        # Roundtrip is limited by binary64 transformation roundoff, not remeshing.
        restored=((transformed.points*1000-T[:3,3])@R-g['origin_source_mm'])*.001
        err=float(np.max(np.abs(restored-source.points)));assert err<1e-15
        path=out/(name+('.vtu' if name in ['flow','mesh'] else '.vtp'));path.parent.mkdir(exist_ok=True);transformed.save(path)
        records[str(path.relative_to(out))]=dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),coordinate_roundtrip_max_m=err)
    (out/'POSE_IDENTITY.json').write_text(json.dumps(dict(candidate_id=int(pose),source_flow_sha256=active['flow_sha256'],
        coordinate_unit='m',velocity_unit='m/s',pressure_and_WSS_unit='Pa',source_mm_to_pose_mm=T.tolist(),
        physical_interpretation='Same no-gravity solution under a rigid coordinate change; not a second independent CFD or gravity simulation',
        unchanged_scalar_fields=True,properly_rotated_vectors=True,files=records),indent=2)+'\n')
    print('POSE_FIELDS',pose,records,flush=True)
