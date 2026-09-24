import pytest
from particle_3d.particle82_batch import sorted_unique

def test_merge_does_not_use_completion_order():
    records=[dict(particle_id=i,value=str(i)) for i in [13,2,100,1]]
    assert sorted_unique(records)==sorted_unique(list(reversed(records)))
    assert [r['particle_id'] for r in sorted_unique(records)]==[1,2,13,100]

def test_duplicate_ids_rejected():
    with pytest.raises(ValueError,match='Duplicate'):sorted_unique([dict(particle_id=1),dict(particle_id=1)])
