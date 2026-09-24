import json
import pytest
from particle_3d.particle82_diagnostics import basin_rows

@pytest.mark.parametrize('nested',[False,True])
def test_real_shard_layouts(tmp_path,nested):
    folder=tmp_path/'shards' if nested else tmp_path
    folder.mkdir(exist_ok=True)
    (folder/'tracer_0000000_0000000.json').write_text(json.dumps({'rows':[{'tracer_id':0,'outlet':'OUTLET_01'}]}))
    assert basin_rows(tmp_path)==[{'tracer_id':0,'outlet':'OUTLET_01'}]

def test_missing_and_duplicate_shards_fail(tmp_path):
    with pytest.raises(ValueError,match='No point'):basin_rows(tmp_path)
    (tmp_path/'tracer_bad.json').write_text(json.dumps({'rows':[{'tracer_id':0},{'tracer_id':0}]}))
    with pytest.raises(ValueError,match='complete/unique'):basin_rows(tmp_path)
