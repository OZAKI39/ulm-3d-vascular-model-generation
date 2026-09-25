from pathlib import Path
import os,sys,subprocess,threading,json,shutil
R=Path(__file__).resolve().parents[1];OLD=Path('/workspace/hemocell_gpu_poc/formal_step3c_gpu_20260914_181003')
sys.dont_write_bytecode=True;sys.path.insert(0,str(R/'scripts'))
from evaluator_service import serve_evaluations
source=R/'provenance/handshake_probe.cpp'
source.write_text('#include "../source/formal_control.hpp"\n#include <iostream>\nint main(int argc,char**argv){try{formalExternalEvaluate(argv[1],5000);std::cout<<"HANDSHAKE_PASS\\n";return 0;}catch(std::exception const&e){std::cout<<e.what()<<"\\n";return 2;}}\n')
exe=R/'provenance/handshake_probe'
subprocess.run(['/usr/bin/g++','-O2','-std=c++17',str(source),'-o',str(exe)],check=True)
results=[]
for kind in ['real_saved_5000','forced_evaluator_failure']:
    F=R/'provenance'/('cpu_handshake_'+kind);F.mkdir()
    for d in ['diagnostics','scripts','logs','provenance','evaluations','tmp']:(F/d).mkdir()
    (F/'contracts').symlink_to(OLD/'contracts',target_is_directory=True)
    (F/'frozen_contracts').symlink_to(F/'contracts',target_is_directory=True)
    (F/'reference').symlink_to(OLD/'reference',target_is_directory=True)
    (F/'provenance/INPUT_RECEIPT.json').symlink_to(OLD/'provenance/INPUT_RECEIPT.json')
    (F/'diagnostics/field_samples').symlink_to(OLD/'diagnostics/field_samples',target_is_directory=True)
    (F/'diagnostics/flow_history.csv').symlink_to(OLD/'diagnostics/flow_history.csv')
    for p in (R/'scripts').glob('*.py'):
        if p.name!='monitor_online.py':(F/'scripts'/p.name).symlink_to(p)
    if kind=='real_saved_5000':shutil.copyfile(R/'scripts/monitor_online.py',F/'scripts/monitor_online.py')
    else:(F/'scripts/monitor_online.py').write_text('raise RuntimeError("Intentional CPU-only failure-path test")\n')
    done=threading.Event();thread=threading.Thread(target=serve_evaluations,args=(F,done),daemon=True);thread.start()
    try:proc=subprocess.run([str(exe),str(F)],capture_output=True,text=True,timeout=30)
    finally:done.set();thread.join(timeout=5)
    expected=0 if kind=='real_saved_5000' else 2
    assert proc.returncode==expected,(kind,proc.stdout,proc.stderr)
    if expected==0:
        decision=(F/'diagnostics/decision.txt').read_text().strip()
        assert decision=='5000 0 0'
        assert json.loads((F/'diagnostics/SHORT_CASE_IDENTITY.json').read_text())['status']=='PASS'
        assert json.loads((F/'evaluations/evaluation_5000.json').read_text())['row']['eligible']==False
    results.append(dict(test=kind,returncode=proc.returncode,output=proc.stdout,status='PASS'))
result=dict(status='PASS',new_GPU_solver_launches=0,new_lattice_timesteps=0,tests=results,mechanism='Exact production C++ handshake header with normal CPU evaluator process; real saved5000fields and forced failure. GPU compute implementation unchanged.')
(R/'provenance/EXTERNAL_EVALUATOR_CPU_TESTS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
