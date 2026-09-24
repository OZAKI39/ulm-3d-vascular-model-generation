import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
read=lambda p:json.loads(p.read_text())


def test_user_requested_cleanup_removed_sources_tests_and_isolated_environment():
    plan=read(ROOT/'reports/stage03/stage018_cleanup_plan.json')
    for name in plan['remove_files']+plan['remove_directories']:
        assert not (ROOT/name).exists() and not (ROOT/name).is_symlink(),name
    assert not list((ROOT/'tests').glob('test_stage018_*.py'))
    for site in ('wsl','remote'):
        result=read(ROOT/f'outputs/stage03/cleanup/{site}.json')
        audit=read(ROOT/f'outputs/stage03/cleanup/{site}_audit.json')
        assert result['all_planned_paths_absent'] and result['production_environment_tool_dependencies_absent']
        assert result['fem_environment_unchanged'] and audit['status']=='PASS'
        assert not audit['active_source_scan_matches']


def test_selected_mesh_core_and_final_historical_science_are_unchanged():
    frozen=read(ROOT/'reports/stage03/history_baseline.json')
    for name,h in frozen['existing_source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    for name,item in frozen['files'].items():
        if name.startswith(('outputs/stage01_7/selected/','reports/stage01_8/')):
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==item['sha256'],name
    assert (ROOT/'reports/stage01_8/REPORT.md').is_file()
