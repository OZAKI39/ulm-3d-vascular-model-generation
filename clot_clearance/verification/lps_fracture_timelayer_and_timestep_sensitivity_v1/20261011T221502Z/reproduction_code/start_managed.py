from common import *
import subprocess
assert (O/'TIMESTEP_DESIGN_FREEZE.json').exists() and not (O/'CAMPAIGN_LAUNCHED.json').exists()
name='clot_timelayer_'+C['run']+'_formal';script=O/'managed_formal_run.sh';script.write_text('#!/bin/bash\nset -e\nsource /workspace/y_clot/env.sh\ncd /workspace/y_clot/project\nexec python -u analysis/timelayer_v1/launch.py\n');script.chmod(0o755)
p=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not p.exists();p.write_text(f'[program:{name}]\ncommand=/bin/bash {script}\ndirectory={R}\nautostart=false\nautorestart=false\nstartsecs=0\nstartretries=0\nstdout_logfile={O}/supervisor.log\nredirect_stderr=true\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n')
code="""import xmlrpc.client,supervisor.xmlrpc
r=xmlrpc.client.ServerProxy('http://127.0.0.1',transport=supervisor.xmlrpc.SupervisorTransport(None,None,'unix:///var/run/supervisor.sock'))
print(r.supervisor.reloadConfig())
r.supervisor.addProcessGroup(NAME)
print(r.supervisor.startProcess(NAME,False))
print(r.supervisor.getProcessInfo(NAME))
""".replace('NAME',repr(name))
q=subprocess.run(['/usr/bin/python3','-c',code],capture_output=True,text=True);(O/'MANAGED_START.txt').write_text(q.stdout+q.stderr);print(q.stdout,q.stderr);assert q.returncode==0
