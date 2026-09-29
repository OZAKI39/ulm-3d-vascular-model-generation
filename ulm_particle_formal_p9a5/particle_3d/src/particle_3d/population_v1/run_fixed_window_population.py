"""Formal fixed-window entry. Times come exclusively from the required contract."""
from pathlib import Path
import argparse,json,platform,sys,time,traceback
import numpy as np
from .environment import load_environment
from .source import PerfusionSource
from .engine import PopulationEngine
from .fixed_window import validate_contract,run_fixed_window,FixedWindowOutput,SCHEMA
from .provenance import digest,content_sha,write_json
from .pairs import contract as pair_contract
from .safety import POLICY


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['environment','output','contract']:ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--backend',choices=['cpu','cuda'],default='cuda');ap.add_argument('--git-commit',required=True)
    ap.add_argument('--master-seed',type=int);ap.add_argument('--realization-id',type=int)
    ap.add_argument('--concentration',type=float);ap.add_argument('--dt',type=float)
    ap.add_argument('--chunk-rows',type=int,default=200000)
    args=ap.parse_args()
    if args.output.exists():raise ValueError('A new run directory is required; historical output cannot be overwritten')
    if args.chunk_rows<1:raise ValueError('Positive chunk size required')
    c=validate_contract(json.loads(args.contract.read_text()));c=dict(c)
    for attr,key in [('master_seed','master_seed'),('realization_id','realization_id'),('concentration','concentration_m3'),('dt','nominal_dt_s')]:
        v=getattr(args,attr)
        if v is not None:c[key]=v
    c['lambda_source_s_inv']=c['concentration_m3']*c['Q_in_m3_s'];validate_contract(c)
    env=load_environment(json.loads(args.environment.read_text()))
    source=PerfusionSource(env,concentration_m3=c['concentration_m3'],seed=c['master_seed'],realization_id=c['realization_id']);validate_contract(c,source)
    identity=dict(git_commit=args.git_commit,current_protocol=SCHEMA,previous_protocol='OUTCOME_DEPENDENT_COVERAGE_STOP',
        experiment_contract=c,experiment_contract_sha256=content_sha(c),input_contract_file_sha256=digest(args.contract),
        flow_sha256=env.flow_sha256,flow_provenance=env.provenance,source_identity=source.identity,
        source_identity_sha256=source.identity_sha256,master_seed=c['master_seed'],realization_id=c['realization_id'],
        backend=args.backend,python=sys.version,platform=platform.platform(),numpy_version=np.__version__,
        unchanged_pair_contract=pair_contract(),unchanged_safety_policy=POLICY,
        code_files_sha256={str(p.relative_to(Path(__file__).parent)):digest(p) for p in sorted(Path(__file__).parent.rglob('*')) if p.is_file() and p.suffix in ['.py','.cu','.cpp']})
    if args.backend=='cuda':
        import torch
        from .gpu_engine import GPUPopulationEngine
        identity.update(torch_version=torch.__version__,gpu=torch.cuda.get_device_name(),cuda=torch.version.cuda)
        engine=GPUPopulationEngine(env,dt=c['nominal_dt_s'],pairs=True,realization_id=c['realization_id'],cache_path=args.environment.parent/'tetra_acceleration.npz')
    else:engine=PopulationEngine(env,dt=c['nominal_dt_s'],pairs=True,realization_id=c['realization_id'])
    writer=FixedWindowOutput(args.output,identity,chunk_rows=args.chunk_rows);engine.recorder=writer
    write_json(args.output/'FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json',c)
    start=time.perf_counter();last=start
    def progress(row):
        nonlocal last
        if time.perf_counter()-last>=5 or row['phase'] in ['INFUSION_OFF','EXPERIMENT_END']:
            write_json(args.output/'PROGRESS.json',dict(row,elapsed_wall_s=time.perf_counter()-start));last=time.perf_counter()
            print(json.dumps(row),flush=True)
    try:
        result=run_fixed_window(engine,source,c,progress)
        proposals=result.pop('source_proposals');timeline=result.pop('population_timeline')
        result.update(elapsed_wall_s=time.perf_counter()-start,backend_counts=dict(getattr(engine,'backend_counts',{})))
        write_json(args.output/'population_timeline.json',timeline)
        import csv
        keys=sorted(set().union(*(r.keys() for r in timeline)))
        with (args.output/'population_timeline.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,keys);w.writeheader();w.writerows(timeline)
        writer.close(result);write_json(args.output/'SUMMARY.json',result)
    except BaseException as error:
        if hasattr(engine,'flush_output'):engine.flush_output()
        (args.output/'error.log').write_text(traceback.format_exc())
        # Partial execution is never advertised as a completed fixed-window experiment.
        result=dict(status='FAILED_EXECUTION_PARTIAL_OUTPUT',error=repr(error),physical_time_s=engine.time,
            counts=dict(engine.counts),active_count=engine.pool.n,current_protocol=SCHEMA)
        writer.close(result);write_json(args.output/'FAILED_EXECUTION.json',result);raise
if __name__=='__main__':main()
