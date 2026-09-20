import ast
from fem_sync_support import ROOT, read


def test_no_particle_solver_or_sampler_was_implemented():
    forbidden={'FrozenFEMField','FlowSampler','FlowSample','RBC','MB','tetra_point_locator','interpolate_velocity','resistance_solver'}
    for top in ['src','scripts/fem_freeze_sync']:
        for p in (ROOT/top).rglob('*.py'):
            tree=ast.parse(p.read_text())
            assert not any(isinstance(n,(ast.ClassDef,ast.FunctionDef,ast.AsyncFunctionDef)) and n.name in forbidden for n in ast.walk(tree))
    copied={r['relative_path']:r for r in read('sync_metadata/copied_from_wsl.json')['files']}
    assert all(p.relative_to(ROOT).as_posix() in copied for p in (ROOT/'src').rglob('*.py'))


def test_no_new_compute_or_two_dimensional_source():
    b=read('frozen_reference/baseline_summary.json')
    assert b['no_CFD_executed']
    assert not (ROOT/'particle').exists() and not (ROOT/'particles').exists()
    assert not list((ROOT/'src').rglob('particle_*.py'))
