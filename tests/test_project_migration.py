import json
from pathlib import Path
import pytest
from vascular_processing.project_paths import frozen_hash,digest,verify_migrated_snapshot,ROOT


def test_saved_digest_accepts_old_symlink_name(tmp_path):
    new=tmp_path/'new';new.write_text('frozen geometry')
    old=tmp_path/'old';old.symlink_to(new)
    assert frozen_hash({str(old):digest(new)},new)==digest(new)


def test_conflicting_alias_digests_rejected(tmp_path):
    new=tmp_path/'new';new.write_text('data');a=tmp_path/'a';b=tmp_path/'b';a.symlink_to(new);b.symlink_to(new)
    with pytest.raises(KeyError):frozen_hash({str(a):'first',str(b):'second'},new)


def test_migration_ledger_requires_original_and_current_bytes(tmp_path):
    code=tmp_path/'loader.py';code.write_text('before');expected=digest(code)
    original=tmp_path/'original.txt';original.write_bytes(code.read_bytes());code.write_text('path-only edit')
    ledger=tmp_path/'migration/20260926/path_compatibility_edits.json';ledger.parent.mkdir(parents=True)
    ledger.write_text(json.dumps({'loader.py':dict(before_sha256=expected,after_sha256=digest(code),original_copy='original.txt')}))
    result=verify_migrated_snapshot({str(code):expected},tmp_path)
    assert result['passed'] and result['migration_only_changes']==[str(code)] and not result['all_content_unchanged']
    original.write_text('changed evidence')
    assert not verify_migrated_snapshot({str(code):expected},tmp_path)['passed']


def test_unlisted_binary_change_is_never_accepted(tmp_path):
    model=tmp_path/'model.stl';model.write_bytes(b'accepted model');expected=digest(model);model.write_bytes(b'changed')
    result=verify_migrated_snapshot({str(model):expected},tmp_path)
    assert not result['passed'] and result['unexpected_changes']==[str(model)]


def test_real_candidate_models_and_slices_still_match_migration_inventory():
    folder=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/final_abs_casting_mold/candidate_0_and_15'
    for line in (folder/'SHA256SUMS').read_text().splitlines():
        expected,name=line.split('  ',1);assert digest(folder/name)==expected


def test_topbrain_alias_is_internal_to_new_project():
    assert (ROOT/'data/TopBrain').resolve().is_relative_to(ROOT)
    assert (ROOT/'data/TopBrain/imagesTr_topbrain_mr').is_dir()
