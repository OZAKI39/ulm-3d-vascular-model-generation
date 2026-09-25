from pathlib import Path
import subprocess,os
R=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617');W=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617');PY='/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019/venv/bin/python'
assert not (R/'PHASE_INITIAL_STATE.json').exists()
config=f'''[unix_http_server]
file={W}/s.sock
[supervisord]
logfile={W}/supervisor.log
pidfile={W}/supervisor.pid
childlogdir={W}
[rpcinterface:supervisor]
supervisor.rpcinterface_factory=supervisor.rpcinterface:make_main_rpcinterface
[supervisorctl]
serverurl=unix://{W}/s.sock
[program:initial]
command={PY} -B {R}/scripts/run_phase.py initial
directory={R}
autostart=true
autorestart=false
startsecs=0
stdout_logfile={R}/logs/PHASE_INITIAL.log
redirect_stderr=true
stdout_logfile_maxbytes=0
'''
(W/'supervisord.conf').write_text(config)
p=subprocess.run(['supervisord','-c',str(W/'supervisord.conf')],capture_output=True,text=True);print(p.returncode,p.stdout,p.stderr);p.check_returncode()
