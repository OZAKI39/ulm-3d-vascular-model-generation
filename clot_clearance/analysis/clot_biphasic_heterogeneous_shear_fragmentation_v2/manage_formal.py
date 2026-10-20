from common import *
import subprocess
assert (O/'SOURCE_FREEZE.json').exists()
assert not (O/'FORMAL_SEQUENCE_REQUEST.json').exists(),'Exactly one sequential formal pair, no outcome retry'
assert not any((O/'cases'/m/'history.csv').exists() for m in CASES)
name='clot_biphasic_formal_'+RUN;wrapper=A/'run_formal_sequence.sh'
wrapper.write_text(f'''#!/bin/bash
set -uo pipefail
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH={O}/build:${{LD_LIBRARY_PATH:-}}
cd {R}
python -u {A}/run_cases.py --name uniform
uniform_status=$?
echo "UNIFORM_PROCESS_EXIT $uniform_status"
python -u {A}/run_cases.py --name biphasic
biphasic_status=$?
echo "BIPHASIC_PROCESS_EXIT $biphasic_status"
exit $((uniform_status | biphasic_status))
''');wrapper.chmod(0o755)
config=f'[program:{name}]\ncommand={wrapper}\nautostart=false\nautorestart=false\nstartsecs=0\nredirect_stderr=true\nstdout_logfile={O}/formal.supervisor.log\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n'
path=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not path.exists();path.write_text(config);(O/(name+'.supervisor.conf')).write_text(config)
dump(O/'FORMAL_SEQUENCE_REQUEST.json',dict(cases=CASES,maximum_new_formal_runs=2,each_once=True,source_freeze_sha256=sha(O/'SOURCE_FREEZE.json'),automatic_restart=False,wall_time_budget=None,memory_budget=None))
subprocess.run(['supervisorctl','reread'],check=True);subprocess.run(['supervisorctl','update',name],check=True);subprocess.run(['supervisorctl','start',name],check=True)
