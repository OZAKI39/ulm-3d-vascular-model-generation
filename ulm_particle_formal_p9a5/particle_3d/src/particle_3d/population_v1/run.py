"""Authoritative-flow CLI: exact source events, observed outlet coverage, drain."""
from pathlib import Path
import argparse,json,platform,sys,time,traceback
import numpy as np
from .environment import load_environment
from .engine import PopulationEngine
from .source import PerfusionSource
from .output import PopulationOutput
from .safety import POLICY
from .pairs import contract as pair_contract
from .provenance import write_json,content_sha,digest


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--environment',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--backend',choices=['cpu','cuda'],default='cuda')
    ap.add_argument('--concentration',type=float,default=4e9);ap.add_argument('--seed',type=int,default=2026092909)
    ap.add_argument('--realization',type=int,default=0);ap.add_argument('--minimum-exits',type=int,default=20)
    ap.add_argument('--post-infusion',type=float,default=12.);ap.add_argument('--dt',type=float,default=.0005)
    ap.add_argument('--git-commit',required=True);ap.add_argument('--fixed-infusion',type=float)
    args=ap.parse_args()
    if args.output.exists():raise ValueError('Output must be a new directory; never overwrite an earlier run')
    if args.minimum_exits<1 or args.post_infusion<0:raise ValueError('Invalid stopping/observation parameters')
    env=load_environment(json.loads(args.environment.read_text()))
    source=PerfusionSource(env,concentration_m3=args.concentration,seed=args.seed,realization_id=args.realization)
    if source.rate==0 and args.fixed_infusion is None:raise ValueError('Zero source cannot reach an outlet coverage target')
    contract=dict(name='GLOBAL_TIME_BRAVA_POPULATION_V1',source=source.identity,pair=pair_contract(),
        safety=POLICY,nominal_dt_s=args.dt,minimum_actual_outlet_crossings=args.minimum_exits,
        stop_infusing='FIRST_ACTUAL_PHYSICAL_TIME_ALL_OFFICIAL_OUTLETS_REACH_TARGET',
        fixed_infusion_test_override_s=args.fixed_infusion,post_infusion_observation_s=args.post_infusion,
        completion='FIRST_CENTER_CROSSING_OF_OFFICIAL_OPEN_OUTLET',wall='AUTHORITATIVE_FINITE_TRIANGLES',
        precision='FLOAT64',long_residence='OBSERVATION_CENSORED_NOT_SUPPORTED_STATIONARY',
        pair_cluster_static_scope='SINGLE_PARTICLE_STATIC_GUARDS_REQUIRE_32_ISOLATED_SOLVES; COUPLED_COMPONENTS_NOT_RELABELED_AS_SINGLE_BODY_SUPPORT')
    files={str(p.relative_to(Path(__file__).parent)):digest(p) for p in sorted(Path(__file__).parent.rglob('*'))
           if p.is_file() and '__pycache__' not in p.parts and p.suffix in ['.py','.cpp','.cu']}
    identity=dict(git_commit=args.git_commit,flow_sha256=env.flow_sha256,
        population_contract_sha256=content_sha(contract),population_contract=contract,
        source_distribution_identity=env.distribution.contract,source_identity_sha256=source.identity_sha256,
        master_seed=args.seed,realization_id=args.realization,backend=args.backend,
        code_files_sha256=files,python=sys.version,platform=platform.platform(),numpy_version=np.__version__,
        flow_provenance=env.provenance)
    if args.backend=='cuda':
        import torch
        from .gpu_engine import GPUPopulationEngine
        identity.update(torch_version=torch.__version__,gpu=torch.cuda.get_device_name(),cuda=torch.version.cuda)
        engine=GPUPopulationEngine(env,dt=args.dt,pairs=True,realization_id=args.realization,
            cache_path=args.environment.parent/'tetra_acceleration.npz')
    else:engine=PopulationEngine(env,dt=args.dt,pairs=True,realization_id=args.realization)
    writer=PopulationOutput(args.output,identity);engine.recorder=writer
    write_json(args.output/'POPULATION_CONTRACT.json',contract)
    outlets=sorted(k for k in env.boundaries if k.startswith('OUTLET_'))
    started=time.perf_counter();last_update=started;infusion_stop=None;observed_terminals=0
    coverage={k:0 for k in outlets};terminal_cursor=0
    def progress(status):
        row=dict(status=status,**engine.accounting(),outlets=dict(coverage),infusion_stop_time_s=infusion_stop,
            elapsed_wall_s=time.perf_counter()-started,backend_counts=dict(getattr(engine,'backend_counts',{})))
        write_json(args.output/'PROGRESS.json',row)
        for s in writer.streams.values():
            if not s.closed:s.flush()
        print(json.dumps(row,ensure_ascii=False),flush=True)
    try:
        while True:
            if infusion_stop is not None and (not engine.pool.n or engine.time>=infusion_stop+args.post_infusion):break
            row=source.peek();ratio=engine.time/engine.dt
            tick=(int(np.floor(ratio+8*abs(np.spacing(ratio))))+1)*engine.dt
            target=min(tick,row['birth_time_s'] if row is not None else float('inf'))
            if infusion_stop is not None:target=min(target,infusion_stop+args.post_infusion)
            if args.fixed_infusion is not None and infusion_stop is None:target=min(target,args.fixed_infusion)
            engine.advance(target)
            # Terminal records are chronologically appended at the real crossing
            # time, which can occur before the end of this integration interval.
            for terminal in engine.terminals[terminal_cursor:]:
                if terminal['status']=='EXITED':coverage[terminal['outlet']]+=1
                if infusion_stop is None and all(v>=args.minimum_exits for v in coverage.values()):
                    infusion_stop=terminal['physical_time_s'];source.stop(infusion_stop)
            terminal_cursor=len(engine.terminals)
            if args.fixed_infusion is not None and infusion_stop is None and engine.time>=args.fixed_infusion:
                infusion_stop=args.fixed_infusion;source.stop(infusion_stop)
            if infusion_stop is None and row is not None and engine.time==row['birth_time_s']:source.admit(engine)
            if time.perf_counter()-last_update>=5:
                progress('INFUSING' if infusion_stop is None else 'OBSERVING_AFTER_INFUSION');last_update=time.perf_counter()
        p=engine.pool
        if p.n:
            ids=p.particle_id[:p.n].cpu().tolist() if args.backend=='cuda' else p.particle_id[:p.n].tolist()
            for pid in ids:
                birth=engine.births[pid];terminal=dict(particle_id=pid,source_event_id=birth['source_event_id'],
                    realization_id=args.realization,birth_time_s=birth['birth_time_s'],physical_time_s=engine.time,
                    residence_time_s=engine.time-birth['birth_time_s'],status='OBSERVATION_CENSORED',
                    terminal_reason='POST_INFUSION_FINITE_OBSERVATION_LIMIT',outlet=None,exit_time_s=None)
                engine.terminals.append(terminal);writer.terminal(terminal);engine.counts['censored_count']+=1
            if args.backend=='cuda':
                mask=__import__('torch').ones(p.n,device=p.device,dtype=__import__('torch').bool)
                engine.safety.remove(mask);p.remove_mask(mask)
            else:p.remove(ids)
        if hasattr(engine,'flush_output'):engine.flush_output()
        summary=dict(status='COMPLETED_WITH_EXPLICIT_TERMINALS',**engine.accounting(),outlets=coverage,
            target_coverage_pass=all(v>=args.minimum_exits for v in coverage.values()),infusion_stop_time_s=infusion_stop,
            source_stop_rule_is_data_dependent=True,elapsed_wall_s=time.perf_counter()-started,
            backend_counts=dict(getattr(engine,'backend_counts',{})),
            scientific_integrity_pass=engine.counts['solver_failed_count']==0)
        manifest=writer.close(summary);write_json(args.output/'SUMMARY.json',summary);progress(summary['status'])
    except BaseException as error:
        if hasattr(engine,'flush_output'):engine.flush_output()
        (args.output/'error.log').write_text(traceback.format_exc())
        try:failed_accounting=engine.accounting()
        except AssertionError:
            failed_accounting=dict(engine.counts,physical_time_s=engine.time,active_count=engine.pool.n,conservation_pass=False,
                reason='INTERRUPTED_DURING_NONATOMIC_TERMINAL_AUDIT; PARTIAL_RUN_NOT_VALID_OUTPUT')
        summary=dict(status='FAILED_EXECUTION_PARTIAL_STATES_PRESERVED',error=repr(error),**failed_accounting,
            backend_counts=dict(getattr(engine,'backend_counts',{})),outlets=coverage,infusion_stop_time_s=infusion_stop)
        writer.close(summary);write_json(args.output/'FAILED_EXECUTION.json',summary)
        raise


if __name__=='__main__':main()
