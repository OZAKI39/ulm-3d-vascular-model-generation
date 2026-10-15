from common import *
import subprocess
assert not (O/'FORMAL_PILOT_LAUNCHED.json').exists()
assert (O/'V3_PILOT_FREEZE.sha256').exists()
name='clot_v3_softening_'+RUN
script=O/'managed_run.sh';script.write_text('#!/bin/bash\nset -e\nsource /workspace/y_clot/env.sh\ncd '+str(A)+'\nexec python -u launch.py pilot\n');script.chmod(0o755)
config=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not config.exists()
config.write_text(f'[program:{name}]\ncommand=/bin/bash {script}\ndirectory={A}\nautostart=false\nautorestart=false\nstartsecs=0\nstartretries=0\nstdout_logfile={O}/supervisor.log\nredirect_stderr=true\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n')
code="""import xmlrpc.client,supervisor.xmlrpc
r=xmlrpc.client.ServerProxy('http://127.0.0.1',transport=supervisor.xmlrpc.SupervisorTransport(None,None,'unix:///var/run/supervisor.sock'))
print(r.supervisor.reloadConfig());r.supervisor.addProcessGroup(NAME);print(r.supervisor.startProcess(NAME,False));print(r.supervisor.getProcessInfo(NAME))
""".replace('NAME',repr(name))
p=subprocess.run(['/usr/bin/python3','-c',code],capture_output=True,text=True);(O/'MANAGED_START.txt').write_text(p.stdout+p.stderr);print(p.stdout,p.stderr);assert p.returncode==0
dump(O/'MANAGED_JOB.json',dict(name=name,config=str(config),script=str(script),autorestart=False,one_formal_launch=True))
