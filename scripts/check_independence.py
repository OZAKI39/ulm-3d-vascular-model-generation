"""Exercise imports and frozen data with old workspaces denied by a Python audit hook.

Run with .venv/bin/python -I scripts/check_independence.py OUTPUT_JSON.
No CFD or trajectory integration is launched.
"""
from pathlib import Path
import json
import os
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
BLOCKED = tuple(str(ROOT.parent / p) for p in (
    'computation_examples', 'temp_storage', 'vascular_printing', 'new_human_vessel_flow',
    'new_vessel_printing', 'new_clot_simulation', 'new_microbubble_experiment',
    'new_microbubble_aggregation', 'github_sync_audits'))
denied = []


def audit(event, args):
    if event not in ('open', 'os.listdir', 'os.scandir') or not args:
        return
    p = args[0]
    if not isinstance(p, (str, bytes, os.PathLike)):
        return
    path = os.path.realpath(os.fsdecode(p))
    if any(path == base or path.startswith(base + '/') for base in BLOCKED):
        denied.append({'event': event, 'path': path})
        raise PermissionError('Historical/external project access denied: ' + path)


def main():
    output = Path(sys.argv[1]).resolve()
    if output.exists():
        raise FileExistsError(output)
    sys.addaudithook(audit)
    checks = {}
    from mouse_vessel_flow.validation import independent_imports, field_check, wss_check
    from mouse_vessel_flow.provenance import verify_source_integrity, protected_originals
    from mouse_vessel_flow.network import validate_design
    from mouse_vessel_flow.geometry import retained_roi
    from mouse_vessel_flow.paths import CASE_NAMES, case_config
    from mouse_vessel_flow.particles.analysis import verify_tracks

    def check(key, function):
        try:
            checks[key] = function()
        except Exception as exc:
            checks[key] = {'status': 'FAIL', 'error': repr(exc), 'traceback': traceback.format_exc()}
        print(key, checks[key]['status'], flush=True)

    check('imports', independent_imports)
    check('current_sources', lambda: {'status': 'PASS', 'count': len(verify_source_integrity())})
    check('protected_sources', lambda: {'status': 'PASS', 'count': len(protected_originals())})
    check('geometry', lambda: {'status': 'PASS', 'anchor_id': retained_roi().anchor_id})
    check('network', validate_design)
    check('wss', wss_check)
    for name in CASE_NAMES:
        check(name + '_field', lambda name=name: field_check(name))
        if case_config(name).get('particles'):
            check(name + '_particles', lambda name=name: verify_tracks(name))
    summary = {'status': 'PASS' if not denied and all(r['status'] == 'PASS' for r in checks.values()) else 'FAIL',
               'denied_access_attempts': denied, 'blocked_roots': BLOCKED, 'checks': checks,
               'new_CFD': 'NOT_RUN', 'new_formal_trajectories': 'NOT_RUN',
               'external_native_installation': '/home/lzy/projects/new_vessel_flow_solvers',
               'limit': 'Import, retained geometry/network, frozen field/trajectory/WSS reads only; no new full solve claimed.'}
    output.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    return 0 if summary['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
