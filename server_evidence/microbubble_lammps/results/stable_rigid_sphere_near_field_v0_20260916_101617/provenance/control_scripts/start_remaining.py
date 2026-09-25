from pathlib import Path
import subprocess
R=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617');W=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617');PY='/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019/venv/bin/python'
phase='remaining'
import json
assert json.loads((R/'contracts/SELECTED_C_GAP.json').read_text())['selected_C_gap']==.4
assert not (R/('PHASE_'+phase.upper()+'_STATE.json')).exists()
p=W/'supervisord.conf';s=p.read_text();assert '[program:'+phase+']' not in s
s+=f'''\n[program:{phase}]
command={PY} -B {R}/scripts/run_phase.py {phase}
directory={R}
autostart=true
autorestart=false
startsecs=0
stdout_logfile={R}/logs/PHASE_{phase.upper()}.log
redirect_stderr=true
stdout_logfile_maxbytes=0
''';p.write_text(s)
for command in ['reread','update']:
 r=subprocess.run(['supervisorctl','-c',str(p),command],capture_output=True,text=True);print(r.stdout,r.stderr);r.check_returncode()
