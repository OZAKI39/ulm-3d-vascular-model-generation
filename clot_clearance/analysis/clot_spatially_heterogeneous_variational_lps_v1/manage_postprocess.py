from common import *
import subprocess,shlex,sys
assert (O/'SOURCE_FREEZE.json').exists()
suffix=('_resume'+(sys.argv[2] if len(sys.argv)>2 else '1')) if len(sys.argv)>1 else '';extra=['--resume-from',sys.argv[1]] if suffix else []
name='clot_heterogeneous_postprocess_'+RUN+suffix
assert not (Path('/etc/supervisor/conf.d')/(name+'.conf')).exists()
wrapper=A/('run_postprocess'+suffix+'.sh');wrapper.write_text('#!/bin/bash\nset -euo pipefail\nsource /workspace/y_clot/env.sh\nexport OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\ncd '+str(R)+'\nexec '+shlex.join(['python','-u',str(A/'postprocess.py')]+extra)+'\n');wrapper.chmod(0o755)
config=f'[program:{name}]\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nredirect_stderr=true\nstdout_logfile={O}/postprocess{suffix}.supervisor.log\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n';p=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not p.exists();p.write_text(config);(O/(name+'.supervisor.conf')).write_text(config)
subprocess.run(['supervisorctl','reread'],check=True);subprocess.run(['supervisorctl','update',name],check=True);subprocess.run(['supervisorctl','start',name],check=True)
