import json
import numpy as np
import pytest
from particle_3d.particle_shapes import Sphere,Ellipsoid,Capsule
from particle_3d.rbc_orientation import rotation_matrix
from particle_3d.wall_geometry import WallGeometry
from particle_3d.wall_gap import wall_gap


@pytest.mark.parametrize('mode',['SPHERE_MB','FREE_OBLATE','CAPILLARY_DEFORMED'])
def test_real_wall_static_finite_shape_gap(p3_repo,real_wall,mode):
    rows=json.loads((p3_repo/'particle_3d/reports/particle3/data/08_real_wall_static.json').read_text())
    chosen=[r for r in rows if r['shape_mode']==mode];assert len(chosen)==5
    for row in chosen:
        tri=real_wall.triangles[row['triangle_id']];np.testing.assert_array_equal(tri,row['triangle_m'])
        patch=WallGeometry(tri[None]);g=row['geometry'];center=row['center_m']
        shape=Sphere(center,g['radius_m']) if mode=='SPHERE_MB' else Ellipsoid(center,g['axes_m'],rotation_matrix(g['q'])) if mode=='FREE_OBLATE' else Capsule(center,g['axis'],g['R_cap_m'],g['L_cap_m'])
        gap=wall_gap(shape,patch)
        assert abs(gap.gap_m-row['analytic_patch_gap_m'])<=gap.roundoff_m
        assert gap.gap_m>0
        assert abs(wall_gap(shape,real_wall).gap_m-row['full_wall_gap_m'])<=gap.roundoff_m


@pytest.mark.parametrize('particle',['MB','RBC'])
def test_real_wall_controlled_contact_and_tangential_sliding(p3_repo,particle):
    rows=json.loads((p3_repo/'particle_3d/reports/particle3/data/09_controlled_summary.json').read_text())
    rows=[r for r in rows if r['particle_type']==particle];assert len(rows)==20
    assert {r['direction'] for r in rows}=={'APPROACH','OBLIQUE','TANGENT','SEPARATING'}
    for row in rows:
        assert row['min_gap_m']>=-row['roundoff_m']
        speed=np.linalg.norm(row['free_velocity_m_s'])
        assert row['max_normal_error_m_s']<=2048*np.finfo(float).eps*speed
        assert row['max_tangential_error_m_s']<=2048*np.finfo(float).eps*speed
        assert row['time_coverage_error_s']<16*np.finfo(float).eps
