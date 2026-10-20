from common import *
import sys,subprocess,shlex
mode=sys.argv[1];assert mode in ['stokes','inverse']
if mode=='inverse':
 q=json.loads((O/'cases/stokes/SUMMARY.json').read_text());assert q['status']=='COMPLETED' and q['max_medium_components']==0
name='clot_shear_'+mode+'_20261017';wrapper=A/('run_'+mode+'.sh')
args=['python','-u',str(A/'run_shear_cases.py'),'--name',mode]
wrapper.write_text('#!/bin/bash\nset -euo pipefail\nsource /workspace/y_clot/env.sh\nexport OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\nexport LD_LIBRARY_PATH='+str(LIB.parent)+':${LD_LIBRARY_PATH:-}\ncd '+str(R)+'\nexec '+shlex.join(args)+'\n');wrapper.chmod(0o755)
config=f'[program:{name}]\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nredirect_stderr=true\nstdout_logfile={O}/{mode}.supervisor.log\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n'
path=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not path.exists();path.write_text(config);(O/(name+'.supervisor.conf')).write_text(config)
subprocess.run(['supervisorctl','reread'],check=True);subprocess.run(['supervisorctl','update',name],check=True);subprocess.run(['supervisorctl','start',name],check=True)
