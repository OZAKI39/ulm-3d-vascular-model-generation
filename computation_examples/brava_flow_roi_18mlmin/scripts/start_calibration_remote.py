from pathlib import Path
import subprocess
R=Path('/workspace/brava_flow_roi_18mlmin_20260928')
conf=R/'supervisord.conf'
assert not conf.exists(), 'Existing supervisor must be inspected, never overwritten'
conf.write_text(f'''[unix_http_server]
file={R}/supervisor.sock
chmod=0600
[supervisord]
logfile={R}/logs/supervisord.log
pidfile={R}/supervisord.pid
childlogdir={R}/logs
[rpcinterface:supervisor]
supervisor.rpcinterface_factory=supervisor.rpcinterface:make_main_rpcinterface
[supervisorctl]
serverurl=unix://{R}/supervisor.sock
[program:brava_calibration_gpu]
command=/root/particle8_2_runs/env/bin/python -B {R}/scripts/solve_remote.py --case {R}/cases/calibration_equal_split --reference-root {R}
directory={R}
environment=PYTHONDONTWRITEBYTECODE="1",OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1"
autostart=false
autorestart=false
startsecs=1
stopasgroup=true
killasgroup=true
stdout_logfile={R}/logs/calibration_runner.log
redirect_stderr=true
stdout_logfile_maxbytes=0
''')
subprocess.run(['supervisord','-c',str(conf)],check=True)
subprocess.run(['supervisorctl','-c',str(conf),'start','brava_calibration_gpu'],check=True)
