from common import *
import subprocess
name='clot_multivortex_postprocess_'+RUN;wrapper=A/'run_postprocess.sh'
wrapper.write_text('#!/bin/bash\nset -euo pipefail\nsource /workspace/y_clot/env.sh\nexport OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\ncd '+str(R)+'\nexec python -u '+str(A/'postprocess.py')+'\n');wrapper.chmod(0o755)
config=f'[program:{name}]\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nredirect_stderr=true\nstdout_logfile={O}/postprocess.supervisor.log\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n'
path=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not path.exists();path.write_text(config)
(O/(name+'.supervisor.conf')).write_text(config)
subprocess.run(['supervisorctl','reread'],check=True);subprocess.run(['supervisorctl','update',name],check=True);subprocess.run(['supervisorctl','start',name],check=True)
