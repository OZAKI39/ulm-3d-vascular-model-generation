from stage03_helpers import *
from fem3d.vascular_diagnostics import distribution

def test_physical_centroid_locations_and_face_neighbors_are_finite():
    m=qc()['residual_cells'];source=read('inputs/stage03/residual_source_evidence.json')['records']
    assert m['count']==len(source)==3 and m['all_finite'] and m['automatic_ratio_threshold'] is None
    assert sorted(tuple(r['source_centroid_m']) for r in m['cells'])==sorted(tuple(r['centroid_m']) for r in source)
    for r in m['cells']:
        assert r['centroid_match_displacement_m']<=1e-15 and not r['source_array_index_used_for_location']
        assert r['neighbor_count']>=1 and min(r['reference_barycentric_coordinates'])>0
        for name,value in r['values'].items():
            assert np.isfinite(value).all()
            stats=r['neighbors'][name]
            assert np.all(np.asarray(stats['min'])<=stats['median']) and np.all(np.asarray(stats['median'])<=stats['max'])

def test_zero_neighbor_median_is_undefined_not_infinite_or_automatic_failure():
    r=distribution([0,0,0],2)
    assert r['residual_over_neighbor_median'] is None and not r['ratio_defined']


def test_residual_cells_relocated_even_though_flow_monitor_cannot_run():
    d=read('outputs/stage03/reference/qc/singularity_diagnosis.json')
    assert d['all_three_cells_relocated'] and len(d['residual_cell_relocation'])==3
    for cell in d['residual_cell_relocation']:
        assert cell['displacement_m']<=1e-15 and not cell['source_array_index_used_for_location']
        assert cell['flow_samples']=='UNAVAILABLE_SOLVER_FAILED'


def test_pressure_support_diagnosis_is_recomputed_from_actual_connectivity():
    from fem3d.vascular_rank_diagnostics import constrained_pressure_support
    m=np.load(ROOT/config()['mesh']['source']/'mesh/volume_mesh.npz')
    current=constrained_pressure_support(m)
    saved=read('outputs/stage03/reference/qc/singularity_diagnosis.json')['topological_evidence']
    assert current==saved
    assert all(r['free_velocity_support_cells']==0 for r in current['uncoupled_pressure_vertices'])
