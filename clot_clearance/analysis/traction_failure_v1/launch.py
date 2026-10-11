"""Start only this round's explicitly named finite offline jobs, without restart."""
from analysis.traction_failure_v1.common import *
import subprocess,sys,shlex
label,module=sys.argv[1:];assert label in ['B_mapping','B_support_repair','A_offline','figures','figures_fixed','figures_final']
name=f'y_clot_traction_{label}_{RUN}';log=COMMON/'logs';log.mkdir(exist_ok=True);wrapper=COMMON/(label+'.sh');conf=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not conf.exists()
wrapper.write_text('#!/bin/bash\nset -e\nsource /workspace/y_clot/env.sh\ncd /workspace/y_clot/project\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1\nexec python -m '+shlex.quote(module)+'\n');wrapper.chmod(0o755)
conf.write_text(f'[program:{name}]\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nstopasgroup=true\nkillasgroup=true\nstdout_logfile={log/label}.log\nredirect_stderr=true\nstdout_logfile_maxbytes=0\n')
subprocess.run(['supervisorctl','reread'],check=True,capture_output=True);subprocess.run(['supervisorctl','update',name],check=True);subprocess.run(['supervisorctl','start',name],check=True)
save(COMMON/(label+'_LAUNCH.json'),dict(name=name,command=['python','-m',module],wrapper_sha256=sha(wrapper),supervisor_config=str(conf),autostart=False,autorestart=False))
