from pathlib import Path
import subprocess
R=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617');W=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617');PY='/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019/venv/bin/python'
assert not (R/'PHASE_PASSIVE_PAIR_STATE.json').exists(), 'Phase already exists; inspect rather than restart'
assert not (R/'cases/J2_passive_separated_pair/RUN_RECEIPT.json').exists()
config=f'''[unix_http_server]
file={W}/j2.sock
[supervisord]
logfile={W}/j2_supervisor.log
pidfile={W}/j2_supervisor.pid
childlogdir={W}
[rpcinterface:supervisor]
supervisor.rpcinterface_factory=supervisor.rpcinterface:make_main_rpcinterface
[supervisorctl]
serverurl=unix://{W}/j2.sock
[program:passive_pair]
command={PY} -B {R}/scripts/run_phase.py passive_pair
directory={R}
autostart=true
autorestart=false
startsecs=0
stdout_logfile={R}/logs/PHASE_PASSIVE_PAIR.log
redirect_stderr=true
stdout_logfile_maxbytes=0
'''
(W/'j2_supervisord.conf').write_text(config)
p=subprocess.run(['supervisord','-c',str(W/'j2_supervisord.conf')],capture_output=True,text=True);print(p.returncode,p.stdout,p.stderr);p.check_returncode()
