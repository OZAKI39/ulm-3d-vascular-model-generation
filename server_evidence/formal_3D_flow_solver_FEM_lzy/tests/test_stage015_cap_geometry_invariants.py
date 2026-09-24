import numpy as np
import pytest
from fem3d.cap_remesh import check_rim,port_geometry
from stage015_helpers import ROOT,OUT,REPORT,CANDIDATES,read,square


@pytest.mark.parametrize('candidate',CANDIDATES)
def test_geometry_decision_uses_unmodified_numeric_thresholds(candidate):
    result=read(OUT/candidate/'qc/surface_invariants.json');policy=read(REPORT/'acceptance_policy.json')
    mesh=np.load(OUT/candidate/'surface/tagged_surface_si.npz')
    contract=read(ROOT/'inputs/stage01_5/baseline_contract.json')
    for port in contract['ports']:
        saved=result['ports'][port['name']]
        tri=mesh['triangles'][mesh['facet_tags']==port['surface_entity_id']]
        fresh=port_geometry(mesh['points_m'],tri,port,np.asarray(saved['interior_vertex_ids']),policy)
        assert fresh['checks']==saved['geometry']['checks']
        assert abs(fresh['area_m2']-saved['geometry']['area_m2'])/fresh['area_m2']<1e-14
    all_pass=True
    for p in result['ports'].values():
        g=p['geometry']; expected=g['relative_area_error']<=1e-12
        assert g['checks']['area']==expected
        assert g['checks']['normal']==(g['normal_dot_source']>=1-1e-12)
        assert g['checks']['centroid']==(g['centroid_displacement_m']<=1e-12)
        assert p['rim']['maximum_rim_displacement_m']==0
        all_pass&=all(g['checks'].values())
    assert (result['status']=='PASS')==all_pass
    assert result['volume_meshing_permitted']==all_pass


def test_reversing_cap_normal_is_rejected():
    p,t=square();port={'plane_origin_m':[.5,.5,0],'outward_normal':[0,0,1],'area_m2':1.,'max_plane_deviation_m':0}
    policy=read(REPORT/'acceptance_policy.json')
    result=port_geometry(p,t[:,[0,2,1]],port,np.array([4]),policy)
    assert result['status']=='FAIL' and not result['checks']['normal']


def test_swapping_port_tags_cannot_match_the_original_rim():
    p,t=square();other=p+np.array([3.,0,0]);points=np.vstack([p,other])
    edges=np.array([[0,1],[1,2],[2,3],[3,0]])
    # Tag swap selects the other patch; a valid patch with a wrong identity must fail.
    with pytest.raises(ValueError,match='Rim edge'): check_rim(points,edges,points,t+len(p))
