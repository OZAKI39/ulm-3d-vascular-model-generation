from common import *
import subprocess
name='clot_mode_i_'+C['run']+'_precheck';script=O/'managed_precheck.sh'
assert not script.exists()
script.write_text('#!/bin/bash\nset -e\nsource /workspace/y_clot/env.sh\ncd /workspace/y_clot/project\nexport OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1\nexec timeout --signal=TERM --kill-after=15 340 python -u analysis/mode_i_v1/precheck.py\n');script.chmod(0o755)
p=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not p.exists();p.write_text(f'[program:{name}]\ncommand=/bin/bash {script}\ndirectory={R}\nautostart=false\nautorestart=false\nstartsecs=0\nstartretries=0\nstdout_logfile={O}/supervisor_precheck.log\nredirect_stderr=true\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n')
code="""import xmlrpc.client,supervisor.xmlrpc
r=xmlrpc.client.ServerProxy('http://127.0.0.1',transport=supervisor.xmlrpc.SupervisorTransport(None,None,'unix:///var/run/supervisor.sock'))
print(r.supervisor.reloadConfig())
r.supervisor.addProcessGroup(NAME)
print(r.supervisor.startProcess(NAME,False))
print(r.supervisor.getProcessInfo(NAME))
""".replace('NAME',repr(name))
q=subprocess.run(['/usr/bin/python3','-c',code],capture_output=True,text=True);(O/'PRECHECK_MANAGED_START.txt').write_text(q.stdout+q.stderr);print(q.stdout,q.stderr);assert q.returncode==0
