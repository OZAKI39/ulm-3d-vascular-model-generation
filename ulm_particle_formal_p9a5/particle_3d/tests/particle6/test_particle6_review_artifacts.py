import json

def test_review_scope_and_evidence(repo):
    report=repo/'particle_3d/reports/particle6';text=(report/'PARTICLE6_REVIEW.md').read_text()
    assert sum(s.startswith('## ') for s in text.splitlines())==13
    for k in range(11):
        note=(report/f'{k:02d}_step_notes.md').read_text()
        for label in ['应该看什么','实际看到什么','有没有异常']:assert label in note and label in text
    v=json.loads((report/'PARTICLE6_VALIDATION.json').read_text())
    assert v['manual_visual_review']=='PENDING_USER_REVIEW' and not v['particle7_started']
    assert v['no_cfd_executed'] and v['real_rbc_lammps_dynamics']=='DEFERRED_DUE_TO_UPSTREAM_PHYSICS'
    for key in ['production_neighbor_cutoff_frozen','production_neighbor_skin_frozen','production_lubrication_cutoff_frozen','production_particle_timestep_frozen','parallel_lammps_validated']:assert v[key] is False
    for name in ['00_lammps_environment','01_state_roundtrip','02_neighbor_equivalence','03_one_step','04_mixed_multistep','05_resistance','06_force_audit','07_restart_mixed','08_metadata_restart','09_real_two_mb','10_neighbor_rebuild']:
        for suffix in ['.csv','.json']:assert (report/'data'/(name+suffix)).stat().st_size>0
