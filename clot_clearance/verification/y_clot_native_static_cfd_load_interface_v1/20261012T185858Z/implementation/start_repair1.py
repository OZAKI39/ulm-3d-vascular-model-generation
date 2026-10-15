from common import *
import subprocess
assert (O/'STATIC_LOAD_DESIGN_FREEZE.json').exists() and not (O/'CAMPAIGN_REPAIR_LAUNCHED.json').exists()
name='clot_static_cfd_'+C['run']+'_repair1';script=O/'managed_repair1_run.sh';script.write_text('#!/bin/bash\nset -e\nsource /workspace/y_clot/env.sh\ncd /workspace/y_clot/project\nexec python -u analysis/static_cfd_v1/launch_repair1.py\n');script.chmod(0o755)
p=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not p.exists();p.write_text(f'[program:{name}]\ncommand=/bin/bash {script}\ndirectory={R}\nautostart=false\nautorestart=false\nstartsecs=0\nstartretries=0\nstdout_logfile={O}/supervisor_repair1.log\nredirect_stderr=true\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n')
code="""import xmlrpc.client,supervisor.xmlrpc
r=xmlrpc.client.ServerProxy('http://127.0.0.1',transport=supervisor.xmlrpc.SupervisorTransport(None,None,'unix:///var/run/supervisor.sock'))
print(r.supervisor.reloadConfig())
r.supervisor.addProcessGroup(NAME)
print(r.supervisor.startProcess(NAME,False))
print(r.supervisor.getProcessInfo(NAME))
""".replace('NAME',repr(name))
q=subprocess.run(['/usr/bin/python3','-c',code],capture_output=True,text=True);(O/'MANAGED_REPAIR_START.txt').write_text(q.stdout+q.stderr);print(q.stdout,q.stderr);assert q.returncode==0
