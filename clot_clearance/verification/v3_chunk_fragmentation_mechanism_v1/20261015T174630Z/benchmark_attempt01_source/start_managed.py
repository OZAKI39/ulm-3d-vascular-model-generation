"""Run one explicitly named audit/benchmark under the existing supervisor."""
from common import *
import sys,subprocess
stage=sys.argv[1];assert stage in {'topology_audit','recover_states','spatial_audit','mixed_mode','benchmark','figures'}
name='clot_chunk_'+stage+'_'+RUN;conf=Path('/etc/supervisor/conf.d')/(name+'.conf')
assert not conf.exists(),'Already launched; inspect its status rather than rerun blindly'
wrapper=O/(stage+'_managed.sh')
wrapper.write_text('#!/bin/bash\nsource /workspace/y_clot/env.sh\nexport OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1\nexport LD_LIBRARY_PATH='+str(LIB.parent)+':"$LD_LIBRARY_PATH"\ncd '+str(A)+'\npython -u -B '+stage+'.py\nrc=$?\nprintf \'%s\\n\' "$rc" > '+str(O/(stage+'_exit_code.txt'))+'\nexit "$rc"\n')
conf.write_text('[program:'+name+']\ncommand=/bin/bash '+str(wrapper)+'\ndirectory='+str(A)+'\nautostart=false\nautorestart=false\nstartsecs=0\nstopasgroup=true\nkillasgroup=true\nredirect_stderr=true\nstdout_logfile='+str(O/(stage+'_process.log'))+'\nstdout_logfile_maxbytes=0\n')
subprocess.run(['supervisorctl','reread'],check=True);subprocess.run(['supervisorctl','update',name],check=True);subprocess.run(['supervisorctl','start',name],check=True)
dump(O/(stage+'_MANAGED.json'),dict(program=name,config_path=str(conf),wrapper_sha256=sha(wrapper),script_sha256=sha(A/(stage+'.py')),new_physical_case=stage=='benchmark'))
