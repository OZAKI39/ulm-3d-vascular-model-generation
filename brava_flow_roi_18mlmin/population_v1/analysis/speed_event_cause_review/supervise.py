"""Reuse the project supervisor, with all state/output confined to this folder."""
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]

if __name__ == '__main__':
    stage = sys.argv[1]
    elapsed = sum(json.loads(p.read_text())['elapsed_s'] for p in
                  (HERE / 'monitoring').glob('*/EXECUTION.json'))
    remaining = 3600 - elapsed
    if remaining <= 0:
        raise SystemExit('CUMULATIVE_COMPUTE_BUDGET_EXHAUSTED')
    spec = importlib.util.spec_from_file_location('existing_supervisor',
        BASE / 'production_workflow/scripts/managed_job.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.ROOT = HERE
    attempt = 1
    name = stage
    while (HERE / 'monitoring' / name).exists():
        attempt += 1
        name = f'{stage}_attempt{attempt}'
    task = dict(name=name, command=[sys.executable, str(HERE / (stage + '.py'))],
                cwd=str(HERE), max_seconds=min(remaining, 2400 if stage == 'diagnose' else 600),
                max_tree_RSS_GiB=12, max_nonreclaimable_GiB=36,
                environment={'PYTHONDONTWRITEBYTECODE': '1'})
    (HERE / 'TASK.json').write_text(json.dumps(task, indent=2) + '\n')
    module.main()
