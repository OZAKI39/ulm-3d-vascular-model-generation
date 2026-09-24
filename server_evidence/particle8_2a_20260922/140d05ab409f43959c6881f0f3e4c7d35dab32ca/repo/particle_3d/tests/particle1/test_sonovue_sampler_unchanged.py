from particle_3d.sonovue_adapter import read_sonovue,sample_single_validation_size


def test_entire_frozen_sonovue_inventory_unchanged(sonovue_root):
    _,_,before=read_sonovue(sonovue_root)
    sample_single_validation_size(sonovue_root,seed=20260920)
    _,_,after=read_sonovue(sonovue_root)
    assert before==after and len(after)==24
