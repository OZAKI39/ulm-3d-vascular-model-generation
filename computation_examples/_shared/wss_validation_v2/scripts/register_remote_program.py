"""Register one fresh validation CFD run with the existing private supervisor.

Invoke on vast4090; never changes a system supervisor or a production case.
"""
from pathlib import Path
import argparse, json, subprocess, time

p = argparse.ArgumentParser()
p.add_argument('--name', required=True)
p.add_argument('--case', required=True)
p.add_argument('--wait-for')
a = p.parse_args()
root = Path(__file__).resolve().parents[1]
assert root == Path('/workspace/wss_validation_v2_20260927T1230Z')
assert all(c.isalnum() or c == '_' for c in a.name)
case = (root / a.case).resolve()
assert case.is_relative_to(root) and (case / 'input_hashes.json').exists()
assert not ((case/'reports/user_cancellation.json').exists() or (case.parent/'user_cancellation.json').exists()), 'SKIPPED_BY_USER: cannot register a cancelled case'
assert not (case / 'run/solver.log').exists(), 'Fresh case required'
config = root / 'supervisord.conf'
old = config.read_text()
assert '[program:' + a.name + ']' not in old
python = '/root/particle8_2_runs/env/bin/python'
if a.wait_for:
    wait = (root / a.wait_for / 'reports/execution.json').resolve()
    assert wait.is_relative_to(root)
    cmd = f'{python} -B {root}/scripts/run_after.py --wait-for {wait} --case {case}'
else:
    cmd = f'{python} -B {root}/scripts/run_solver.py --case {case}'
block = f'''
[program:{a.name}]
command={cmd}
directory={root}
environment=PYTHONDONTWRITEBYTECODE="1",OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1"
autostart=true
autorestart=false
startsecs=1
stopasgroup=true
killasgroup=true
stdout_logfile={root}/logs/{a.name}.log
redirect_stderr=true
stdout_logfile_maxbytes=0
'''
backup = root / 'logs' / ('supervisord_before_' + a.name + '.conf')
backup.write_text(old)
config.write_text(old + block)
subprocess.run(['supervisorctl', '-c', str(config), 'reread'], check=True)
subprocess.run(['supervisorctl', '-c', str(config), 'update', a.name], check=True)
(root / 'logs' / ('registered_' + a.name + '.json')).write_text(json.dumps(
    dict(case=str(case), command=cmd, registered_unix=time.time(),
         supervisor_configuration=str(config), backup=str(backup)), indent=2) + '\n')
