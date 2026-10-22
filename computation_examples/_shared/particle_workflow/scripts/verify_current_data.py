"""Check current retained scientific files and all formal trajectory receipts."""
from pathlib import Path
import argparse, json, sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import REL, DT, FLOW_SHA, digest, load_births, completion_matches

def verify():
    report = ROOT/REL
    manifest = report/'data/retained_science.sha256'
    checked = 0
    for line in manifest.read_text().splitlines():
        expected, relative = line.split('  ',1)
        path = ROOT/relative
        if not path.is_file() or digest(path) != expected:
            raise ValueError('Retained science missing or changed: '+relative)
        checked += 1
    cohort = json.loads((report/'data/FINAL_FORMAL_COHORT.json').read_text())
    rows = json.loads((report/'data/formal_metrics.json').read_text())
    identity = json.loads((report/'data/worker_scaling_results.json').read_text())['identity']
    assert cohort['events'] == load_births(ROOT)['events'][:len(rows)]
    assert cohort['count'] == len(rows) >= 500 and DT == .001
    assert len({row['particle_id'] for row in rows}) == len(rows)
    for event, row in zip(cohort['events'],rows):
        assert event['particle_id'] == row['particle_id']
        folder = report/row['track_relative_path']
        assert completion_matches(folder,identity,event), folder
        assert json.loads((folder/'trajectory.json').read_text())['integration_config']['dt_s'] == DT
    point = json.loads((report/'data/point_complete.json').read_text())
    assert point['count'] == len(rows) and point['cohort_unchanged']
    return dict(all_sha256_match=True, protected_files=checked, formal_trajectories=len(rows),
                paired_points=point['count'], nominal_dt_s=DT, flow_sha256=FLOW_SHA,
                scope='RETAINED_CURRENT_SCIENCE_AFTER_AUTHORIZED_WSL_CLEANUP',
                old_files_changed=0, historical_snapshot_enforced=False)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = verify()
    text = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text,end='')

if __name__ == '__main__':
    main()
