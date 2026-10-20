from common import *
import subprocess
name='clot_biphasic_validation_'+RUN;wrapper=A/'run_validation.sh'
wrapper.write_text(f'''#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd {A}
HET_LIB={O}/baseline_capacity/liblammps.so LD_LIBRARY_PATH={O}/baseline_capacity:${{LD_LIBRARY_PATH:-}} python -u uniform_probe.py old_uniform
export LD_LIBRARY_PATH={O}/build:${{LD_LIBRARY_PATH:-}}
python -u uniform_probe.py new_uniform
python -u checks.py
python -u path_audit.py
python -u freeze_sources.py
''');wrapper.chmod(0o755)
config=f'[program:{name}]\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nredirect_stderr=true\nstdout_logfile={O}/validation.supervisor.log\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n'
path=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not path.exists();path.write_text(config);(O/(name+'.supervisor.conf')).write_text(config)
subprocess.run(['supervisorctl','reread'],check=True);subprocess.run(['supervisorctl','update',name],check=True);subprocess.run(['supervisorctl','start',name],check=True)
