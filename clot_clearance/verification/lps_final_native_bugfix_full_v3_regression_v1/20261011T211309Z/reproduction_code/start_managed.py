from common import *
import subprocess,sys
assert json.loads((O/'STATIC_PREFLIGHT.json').read_text())['status']=='PASS'
assert not (O/'FORMAL_V3_LAUNCHED.json').exists()
assert sha(B/'liblammps.so')==EXPECTED
name='clot_final_v3_20261011T211309Z'
script=O/'managed_run.sh';script.write_text('#!/bin/bash\nset -e\nsource /workspace/y_clot/env.sh\ncd /workspace/y_clot/project\nexec python -u analysis/final_v3_v1/launch.py\n');script.chmod(0o755)
config=Path('/etc/supervisor/conf.d')/(name+'.conf');assert not config.exists();config.write_text(f'[program:{name}]\ncommand=/bin/bash {script}\ndirectory={R}\nautostart=false\nautorestart=false\nstartsecs=0\nstartretries=0\nstdout_logfile={O}/supervisor.log\nredirect_stderr=true\nstdout_logfile_maxbytes=0\nstopasgroup=true\nkillasgroup=true\n')
dump(O/'IMPLEMENTATION_CHECKS.json',{'before_formal_run':'PASS','physics_module_SHA_frozen':True,'original_deck_byte_identical':True,'new_solver_source_changes':False,'readonly_monitor':'END_OF_STEP only; per-atom reads copied; pair single reads actual neighbor list into private buffer; no force/velocity/topology writes; MPI gather/barrier only; original external callbacks unchanged','sort_configuration':'unchanged default1000, actualnextsort and permutation inspected','guard_audit':'Only owned num_bond arrays inspected; every existing endpoint type and ID checked against independent full topology on changed steps','crash_safe':'CSV/JSONL line buffered, key states synchronous; rank exceptions abort with diagnostics','static_checks_actual_attempts':2,'static_harness_correction':'First static compute succeeded, Python energy.eval(False) scalar was incorrectly indexed during report collection. Corrected helper and repeated same tiny static check; no V3 initiated, no new material case.','formal_full_V3_limit':1,'observer_sha256':sha(O/'observer.so'),'dynamic_library_hard_gate':'actual maps+SHA per MPI rank before executing deck','production_sha256':sha(PROD)})
# Control only our newly created program. No broad update/restart.
code="""import xmlrpc.client,supervisor.xmlrpc,json
r=xmlrpc.client.ServerProxy('http://127.0.0.1',transport=supervisor.xmlrpc.SupervisorTransport(None,None,'unix:///var/run/supervisor.sock'))
print(r.supervisor.reloadConfig())
r.supervisor.addProcessGroup(NAME)
print(r.supervisor.startProcess(NAME,False))
print(r.supervisor.getProcessInfo(NAME))
""".replace('NAME',repr(name))
p=subprocess.run(['/usr/bin/python3','-c',code],capture_output=True,text=True);(O/'MANAGED_START.txt').write_text(p.stdout+p.stderr);print(p.stdout,p.stderr);assert p.returncode==0
