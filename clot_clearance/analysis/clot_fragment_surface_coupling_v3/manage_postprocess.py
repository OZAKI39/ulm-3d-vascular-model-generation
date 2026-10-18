"""One waiting postprocess job per phase; no extra physical run."""
from common import *
import argparse,subprocess,datetime

def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['prepare','finish']);args=p.parse_args();phase=args.phase
    request=O/f'POST_REQUEST_{phase}.json';assert not request.exists(),'Inspect prior postprocess attempt before recovery'
    name=f'clot_surface_v3_post_{phase}_{RUN}';wrapper=A/f'post_{phase}.sh'
    wrapper.write_text(f'#!/bin/bash\nset -euo pipefail\nsource /workspace/y_clot/env.sh\nexport OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\nexport LD_LIBRARY_PATH={V2_O}/build:${{LD_LIBRARY_PATH:-}}\ncd {R}\nexec python -u {A}/wait_and_postprocess.py {phase}\n');wrapper.chmod(0o755)
    cfg=f'[program:{name}]\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nredirect_stderr=true\nstdout_logfile={O}/post_{phase}.supervisor.log\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n'
    path=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not path.exists();path.write_text(cfg);(O/(name+'.supervisor.conf')).write_text(cfg)
    dump(request,dict(phase=phase,name=name,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),physical_runs_started=0))
    subprocess.run(['supervisorctl','reread'],check=True);subprocess.run(['supervisorctl','update',name],check=True);subprocess.run(['supervisorctl','start',name],check=True)

if __name__=='__main__':main()
