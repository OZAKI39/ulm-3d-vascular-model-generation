from particle_3d.particle6_validation import neighbor_case
from particle_3d.particle6_stepper import Particle6Stepper
import pytest

def test_six_shape_pairs_and_exact_filter(mu):
    d=neighbor_case(mu)
    assert len(d['groups'])==6 and len(d['standalone_exact_pairs'])==6
    assert d['standalone_candidates']==d['lammps_candidates']
    assert d['standalone_exact_pairs']==d['lammps_exact_pairs']
    assert d['false_negatives']==d['filtered_mismatch']==0 and d['extra_candidates']>0

def test_query_cannot_drop_required_physics():
    with pytest.raises(ValueError,match='DOES_NOT_COVER'):Particle6Stepper._require_coverage([(1,2)],[])
