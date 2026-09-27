import ast
from pathlib import Path
import particle_3d.injection_method_c as m

def test_production_has_no_future_route_input_or_quota():
    p=Path(m.__file__);text=p.read_text();tree=ast.parse(text)
    for forbidden in ['OUTLET_01','OUTLET_02','OUTLET_03','point_basin','point_tracer','outlet_classifier','outlet_quota']:
        assert forbidden not in text
    imports=[n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
    assert not any(any(w in (x or '') for w in ['routing','tracer','classifier','particle82a_admission']) for x in imports)
    assert set(ast.get_source_segment(text,n) for n in ast.walk(tree) if isinstance(n,ast.Attribute) and n.attr=='classifier')==set()
