import csv

from vascular_processing.bambu_casting_mold import write_slice_comparison


def test_geometry_only_run_has_valid_empty_slice_table(tmp_path):
    output = tmp_path / 'comparison.csv'
    write_slice_comparison([], output)
    with output.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        assert {'candidate_id', 'status', 'support_mode'} <= set(reader.fieldnames)
        assert list(reader) == []


def test_pending_slice_record_is_not_reported_as_printed(tmp_path):
    output = tmp_path / 'comparison.csv'
    write_slice_comparison([dict(orientation_id=15, support_mode='ON',
                                status='BAMBU_VALIDATION_PENDING')], output)
    with output.open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream))
    assert rows[0]['candidate_id'] == '15'
    assert rows[0]['status'] == 'BAMBU_VALIDATION_PENDING'
    assert rows[0]['print_time_seconds'] == ''
