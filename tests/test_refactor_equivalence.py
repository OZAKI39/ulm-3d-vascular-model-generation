"""Behavior regressions against explicitly captured old-source baselines, not reimplementations."""
import json
import numpy as np
import pytest
from mouse_vessel_flow.paths import project_path,new_output
from mouse_vessel_flow.particles.particle65_cases import run_case
from mouse_vessel_flow.particles.particle3_cases import json_safe

def test_full_small_scene_outputs_match_original_bytes():
    expected=project_path('validation/original_small_scene_baseline.json').read_text()
    result={name:run_case(name,.004,.032) for name in ('wall','pair')}
    assert json.dumps(result,default=json_safe,sort_keys=True)+'\n'==expected

def test_native_handoff_kernel_matches_original_frozen_wall_bytes():
    from mouse_vessel_flow.particles.wall_geometry import WallGeometry
    from mouse_vessel_flow.particles.native_adapter import native_closest
    from mouse_vessel_flow.particles.convex_triangle import triangle_closest_many
    wall=WallGeometry.from_frozen(project_path('inputs/reference/fem'))
    triangles=wall.triangles[np.linspace(0,len(wall.triangles)-1,128,dtype=int)]
    rng=np.random.default_rng(20261022)
    for i in range(16):
        point=triangles[i].mean(axis=0)+rng.normal(size=3)*1e-6
        a=triangle_closest_many(point,triangles);b=native_closest(point,triangles)
        assert a[0].tobytes()==b[0].tobytes() and a[1].tobytes()==b[1].tobytes()

def test_retained_geometry_is_the_flow_roi():
    from mouse_vessel_flow.geometry import retained_roi
    roi=retained_roi()
    assert roi.anchor_id==3274
    assert 'fMOST_0_5_6_0_0_6_0001_02_01' in roi.roi_id

def test_frozen_output_guard_rejects_write():
    with pytest.raises(ValueError):new_output(project_path('results/roi_only/new_run'))
    with pytest.raises(ValueError):new_output(project_path('inputs/new_run'))

def test_embedded_tetgen_keeps_original_body_and_legacy_syntax():
    import ast,zipfile
    path=project_path('mouse_vessel_flow/mesh_tools/tetgen_native.py')
    ast.parse(path.read_text(),feature_version=(3,5))
    with zipfile.ZipFile(project_path('inputs/reference/source_identity/original_source_identity.zip')) as z:
        old=ast.parse(z.read('selected/mesh_tools/tetgen_native.py'))
    new=ast.parse(path.read_text());fn=next(n for n in new.body if isinstance(n,ast.FunctionDef) and n.name=='generate')
    start=next(i for i,n in enumerate(old.body) if isinstance(n,ast.With))
    expected=[n for n in old.body[start:] if not (isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id=='factor')]
    begin=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.With))
    assert [ast.dump(n,include_attributes=False) for n in expected]==[ast.dump(n,include_attributes=False) for n in fn.body[begin:-1]]
    from mouse_vessel_flow.mesh_tools.tetgen_native import validate_root
    with pytest.raises(ValueError):validate_root(project_path('results/roi_only'))

@pytest.mark.parametrize('case',['roi_only','best_feasible','network_h0'])
def test_real_source_admission_and_job_binding(case):
    from mouse_vessel_flow.flow_field import load_case
    from mouse_vessel_flow.particles.injection import saved_events
    from mouse_vessel_flow.particles.dynamics import bind_job
    env=load_case(case);event=saved_events(case)[0]
    particle,status,detail=env.checker.check(event,np.asarray(event['birth_center_m']),{})
    assert particle is not None and status=='ACCEPTED'
    job=bind_job(env,{'role':'BIND_ONLY_NO_TRAJECTORY_EXECUTED'},env.case_config['particles']['dt_s'],
                 accelerated=env.case_config['particles']['accelerated'])
    assert callable(job)
    if case!='network_h0':
        import gzip
        from mouse_vessel_flow.particles.injection import source_for_case
        source=source_for_case(case,environment=env)
        report=project_path(env.case_config['particles']['report'])
        with gzip.open(report/'data/proposal_ledger.jsonl.gz','rt') as stream:
            for _ in range(2):
                old=json.loads(next(stream));new=source.proposal(old['source_event_id'])
                # Parent reducer assigns particle ID/time; code identity is deliberately new.
                for key in new:
                    if key not in ('identity_sha256','particle_id'):assert new[key]==old[key],key
