from pathlib import Path
import csv,json
import numpy as np
S=Path(__file__).resolve().parents[1]
results=[]
for name in ['LOW','MEDIUM','CAPACITY_STRESS_LIMIT','FLUX_WEIGHTED']:
    a=S/'runs'/f'{name}_release_mpi1';b=S/'runs'/f'{name}_release_mpi4';checks={}
    for file in ['injection_events.csv','flux_timeseries.csv','INTEGRATION_STAGES.csv','PORT_EVENTS.csv','NEAR_WALL_EVENTS.csv','PENDING_FINAL.csv','PAIR_HYDRODYNAMIC_EVENTS.csv','RETRY_HISTORY.csv']:
        checks[file]=a.joinpath(file).read_bytes()==b.joinpath(file).read_bytes()
    with (a/'TRAJECTORIES.csv').open() as f:r1=list(csv.DictReader(f))
    with (b/'TRAJECTORIES.csv').open() as f:r4=list(csv.DictReader(f))
    checks['trajectory_rows_same']=len(r1)==len(r4)
    maxerr=0.;ownership_differences=0
    for x,y in zip(r1,r4):
        ownership_differences+=x.pop('owner_rank')!=y.pop('owner_rank')
        for k in ['x_m','y_m','z_m']:maxerr=max(maxerr,abs(float(x[k])-float(y[k])))
        assert x==y,(name,x,y)
    checks['trajectory_all_other_fields_bitwise_equal']=True
    state1=json.loads((a/'RUN_STATE.json').read_text());state4=json.loads((b/'RUN_STATE.json').read_text());state1.pop('mpi_ranks');state4.pop('mpi_ranks');checks['terminal_state_equal_except_rank_count']=state1==state4
    assert all(checks.values()),checks
    results.append(dict(case=name,status='PASS',checks=checks,max_position_difference_m=maxerr,owner_rank_label_differences=ownership_differences))
(S/'validation/MPI_CONSISTENCY.json').write_text(json.dumps(dict(status='PASS',comparison='MPI1 versus MPI4 same frozen inputs, source and position seeds; same rank0 controller',cases=results),indent=2)+'\n')
print(json.dumps(results,indent=2))
