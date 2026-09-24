import json,hashlib
from pathlib import Path
from types import SimpleNamespace
import particle_3d.injection_method_c as new
from particle_3d.particle82a_admission import common_event,method_b
from particle_3d.particle82a_geometry import segment_distance,perimeter_edges
from method_c_test_helpers import source

def test_legacy_hash_and_replay_unchanged(distribution):
    root=Path(new.__file__).resolve().parents[3]
    before=json.loads((root/'particle_3d/reports/particle9a2_inlet_sampling/data/protected_before.json').read_text())
    old=root/'particle_3d/src/particle_3d/particle82a_admission.py'
    assert hashlib.sha256(old.read_bytes()).hexdigest()==next(v for k,v in before['files'].items() if k.endswith('/particle_3d/src/particle_3d/particle82a_admission.py'))
    s=source(distribution);geometry=SimpleNamespace(distances=lambda p:(s.wall.nearest_center_triangle(p)[1],s.wall.nearest_center_triangle(p)[1],0.))
    ctx=SimpleNamespace(env=SimpleNamespace(sampler=s.sampler,wall=s.wall,field=SimpleNamespace(locate=lambda p:(0,None))),geometry=geometry,distribution=distribution.original)
    old_event=common_event(41,ctx=ctx);a=method_b(old_event,ctx)
    s.event(1)
    b=method_b(old_event,ctx)
    assert a==b and a['method']=='B'
    assert new.LEGACY_METHOD=='LEGACY_METHOD_B_POSITION_ANCHORED_SIZE_CONDITIONING'
