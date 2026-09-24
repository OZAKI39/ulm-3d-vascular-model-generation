from particle_3d.particle1_audit import check_dependency, P0_COMMIT


def test_accepted_p0_all_recorded_code_and_tests_are_byte_identical(repo):
    record = check_dependency(repo)
    assert record['git_commit'] == P0_COMMIT
    assert record['manual_visual_review'] == 'PASS'
    assert record['manual_review_evidence']['reviewer'] == 'USER_CHAT_REVIEW'
    assert record['manual_review_evidence']['no_enhanced_particle0_test_requested'] is True
    # Exact original file set: no extra Particle-0 tests have been added.
    assert {str(p.relative_to(repo)) for p in (repo/'particle_3d/tests/particle0').glob('test_*.py')} == set(record['test_files'])
