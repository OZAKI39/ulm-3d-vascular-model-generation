"""Explicit version checks for old integration tests after approved UI revisions."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASELINE=json.loads((ROOT/'research_retrospective/production_protected_before.json').read_text())
MAINTENANCE=json.loads((ROOT/'research_retrospective/retained_test_maintenance/changes.json').read_text())
def assert_current_version(path):
    path=Path(path).resolve();name=str(path.relative_to(ROOT));expected=MAINTENANCE.get(name,{}).get('after_sha256',BASELINE.get(str(path)))
    assert expected, 'No frozen version: '+str(path)
    assert hashlib.sha256(path.read_bytes()).hexdigest()==expected, name

def assert_authorized_later_ui_changes(paths):
    allowed={'s1-2_swc_roi_generate_human.py','s1-3_swc_roi_generate_MeVO.py','utils/rodent_vasculature/interactive.py'}
    for path in paths:
        p=Path(path).resolve();assert str(p.relative_to(ROOT)) in allowed
        assert_current_version(p)

def assert_maintained_source_changes(paths,snapshot):
    for path in paths:
        p=Path(path).resolve();name=str(p.relative_to(ROOT));row=MAINTENANCE[name]
        assert row['before_sha256']==snapshot[path]
        assert hashlib.sha256((ROOT/row['original_copy']).read_bytes()).hexdigest()==row['before_sha256']
        assert_current_version(p)
