#!/usr/bin/python3
import csv,sys
from pathlib import Path
sys.dont_write_bytecode=True
from remote_common import *
from convergence import *
from gpu_integrity import compare_short,inspect_snapshot
def main(root,run,step):
    c=contract(root);verify_run_inputs(root,run)
    data=history(run);check_observed_safety(data,c)
    # Independent host reconstruction from saved live GPU fields at original evaluation cadence.
    try:
        if step==5000:
            short=compare_short(root,run);write(run/'diagnostics/SHORT_CASE_IDENTITY.json',short)
        check=inspect_snapshot(root,run,step)
        write(run/f'diagnostics/gpu_integrity_{step}.json',check)
    except Exception as e:
        write(run/'diagnostics/GPU_NUMERICAL_INTEGRITY_FAILURE.json',{'status':'FAIL_GPU_NUMERICAL_INTEGRITY','step':step,'error':str(e)})
        raise
    path=run/'diagnostics/evaluator_state.json'
    state=read(path) if path.exists() else initial_state()
    expected=5000 if state['last_evaluation'] is None else state['last_evaluation']+5000
    if step!=expected:raise ValueError('Skipped/repeated evaluator invocation')
    fv=field_residual(run,step,c) if step>=c['execution']['first_possible_evaluation_step'] else None
    row,state,detail=evaluate_point(data,step,c,state,fv)
    with (run/'diagnostics/convergence_history.csv').open('a',newline='') as f:
        w=csv.DictWriter(f,fieldnames=COLUMNS)
        if f.tell()==0:w.writeheader()
        w.writerow(row);f.flush()
    write(run/f'evaluations/evaluation_{step}.json',{'row':row,'state':state,'detail':detail})
    write(path,state)
    # Decision released only after integrity and frozen scientific gates finish.
    (run/'diagnostics/decision.txt').write_text(f'{step} {int(detail["stop"])} {int(detail["new_candidate"])}\n')
    print('FORMAL_EVALUATION',step,'eligible',row['eligible'],'passed',row['all_gates_pass'],'R_inlet',row['R_inlet'],'failed',row['failed_gates'],'confirmations',state['consecutive_pass_count'],flush=True)
if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3]))
