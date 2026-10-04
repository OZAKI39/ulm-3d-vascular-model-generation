"""Load existing C5/C7 states and reports only; never evaluate the model."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

AUDIT = Path(__file__).resolve().parent
PROJECT = AUDIT.parents[1]


def main(stage):
    paths = set()
    for name in ('pfpd_at1_c5_activation_stationarity_audit_v1', 'pfpd_at1_c7_early_branch_v1'):
        root = PROJECT / 'verification' / name
        bracket = json.loads((root / 'activation/BRACKET.json').read_text())
        for edge in ('lower', 'upper'):
            paths.add(root / 'intact' / bracket[edge]['label'] / 'STATE.npz')
        paths.update(root.glob('intact/*activation*/STATE.npz'))
        paths.update(root.glob('intact/*weak*/STATE.npz'))
    paths.update((PROJECT/'verification/pfpd_at1_c5_stability_v1/states').glob('*/STATE.npz'))
    paths.update((PROJECT/'verification/pfpd_at1_c7_early_branch_v1/branch').glob('*/STATE.npz'))
    states = []
    for path in sorted(paths):
        with np.load(path, allow_pickle=False) as data:
            arrays = {key: data[key] for key in data.files}
            assert 'x' in arrays and 's' in arrays, path
            assert np.isfinite(arrays['x']).all() and np.isfinite(arrays['s']).all(), path
            assert arrays['x'].shape in ((125, 3), (343, 3)), path
            assert arrays['s'].shape in ((5392,), (46059,)), path
            states.append({'path': str(path.relative_to(PROJECT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                           'arrays': {k: {'shape': list(v.shape), 'dtype': str(v.dtype)} for k,v in arrays.items()},
                           'maximum_s': float(arrays['s'].max())})
    reports = []
    for task in ['pfpd_at1_c5_stability_v1', 'pfpd_at1_c7_early_branch_v1', 'pfpd_at1_activation_mesh_dependence_audit_v1']:
        root = PROJECT / 'verification' / task
        for path in [root/'SUMMARY.json', *root.glob('REVIEW*.md')]:
            text = path.read_text()
            if path.suffix == '.json':
                assert isinstance(json.loads(text), dict)
            reports.append({'path': str(path.relative_to(PROJECT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    result = dict(stage=stage, loadable=True, simulations_run=0, states=states, reports=reports)
    if stage == 'AFTER':
        before = json.loads((AUDIT/'SCIENTIFIC_SMOKE_BEFORE.json').read_text())
        result['states_and_reports_identical'] = before['states'] == states and before['reports'] == reports
        assert result['states_and_reports_identical']
    (AUDIT/f'SCIENTIFIC_SMOKE_{stage}.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'stage': stage, 'state_count': len(states), 'report_count': len(reports), 'loadable': True}))


if __name__ == '__main__':
    main(sys.argv[1])
