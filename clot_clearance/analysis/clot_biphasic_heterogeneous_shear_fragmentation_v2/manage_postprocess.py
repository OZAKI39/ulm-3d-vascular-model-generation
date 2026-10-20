from common import *
import subprocess,sys
stage=sys.argv[1];assert stage in ['prepare','finish']
name='clot_biphasic_post_'+stage+'_'+RUN;wrapper=A/('run_post_'+stage+'.sh')
wrapper.write_text(f'''#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYVISTA_OFF_SCREEN=true
cd {A}
exec python -u postprocess.py {stage}
''');wrapper.chmod(0o755)
config=f'[program:{name}]\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nredirect_stderr=true\nstdout_logfile={O}/post_{stage}.supervisor.log\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n'
path=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not path.exists();path.write_text(config);(O/(name+'.supervisor.conf')).write_text(config)
subprocess.run(['supervisorctl','reread'],check=True);subprocess.run(['supervisorctl','update',name],check=True);subprocess.run(['supervisorctl','start',name],check=True)
