import json
from sv_validation.validation import face_mapping

def test_original_semantic_roles_preserved(root):
    mapping=json.loads((root/"configs/face_map.json").read_text())["faces"]
    assert face_mapping(mapping, {"1":"WALL","2":"OUTLET_03","3":"OUTLET_01","4":"INLET","5":"OUTLET_02"})
