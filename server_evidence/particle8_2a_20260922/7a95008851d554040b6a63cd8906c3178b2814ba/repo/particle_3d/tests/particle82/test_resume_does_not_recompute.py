import json
from pathlib import Path
import pytest
from particle_3d.particle82_batch import task_run
from particle_3d.particle82_provenance import atomic_json,sha256
from particle_3d.particle8_replay import canonical_hash

def test_resume_reuses_committed_id(tmp_path,remote,monkeypatch):
    import particle_3d.particle82_batch as b
    config=dict(dt_s=.00025,guard_factor=1);event=dict(particle_id=1)
    def solver(event,**kwargs):
        assert not kwargs.get('force',False),'Completed ID must never be forcibly recomputed'
        return dict(particle_id=1,sample_count=2,samples_sha256='saved')
    monkeypatch.setattr(b,'integrate_one',solver)
    mp=tmp_path/'trajectories/mb_000001.json';atomic_json(mp,{'saved':True})
    receipt=dict(config_sha256=canonical_hash(config),git_commit='test-commit',hostname='test-remote',
                 metadata_sha256=sha256(mp),frozen_input_sha256=remote['frozen_input_sha256'])
    atomic_json(tmp_path/'receipts/mb_000001.json',receipt)
    before=mp.read_bytes();assert task_run(([event],str(tmp_path),config,remote,True))[0]['particle_id']==1
    assert mp.read_bytes()==before

def test_resume_corrupt_metadata_rejected(tmp_path,remote):
    config=dict(dt_s=.00025,guard_factor=1);mp=tmp_path/'trajectories/mb_000001.json';atomic_json(mp,{})
    atomic_json(tmp_path/'receipts/mb_000001.json',dict(config_sha256=canonical_hash(config),
        git_commit='test-commit',hostname='test-remote',metadata_sha256='wrong'))
    with pytest.raises(ValueError,match='metadata mismatch'):task_run(([dict(particle_id=1)],str(tmp_path),config,remote,True))
