"""Fixed-state numerical equivalence and timing; no trajectory integration."""
import json,time,hashlib,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import h5py
from scripts.stage4l.dense_probe import load_case,sha,write
from mevo_mb.coupling.integrators import GlobalSynchronousIntegrator
from mevo_mb.bubbles.constrained_motion import ConstrainedLubrication,horizon_constraints
from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes
from mevo_mb.bubbles.auto_backend import AutoPairBackend
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'results/stage5/motion_solver'

def main():
    out=BASE/'benchmarks'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True,exist_ok=False)
    selected=json.load(open(BASE/'audit_v2/BASELINE_SELECTION.json'));rows=[];checks=[]
    cpuquota=Path('/sys/fs/cgroup/cpu.max').read_text().split();quota=int(cpuquota[0])/int(cpuquota[1]);workers=min(4,max(1,int(quota)))
    for case in ['B2','B3']:
        source=Path(next(r['run'] for r in selected if r['key']==case+'_5ms'))
        config,field,flow,_,inputs=load_case(Path('/workspace/mevo_acoustic_microbubbles'),case)
        config['particles']['translation_backend']='constrained_lubrication_v1'
        tracker=GlobalSynchronousIntegrator(config,field,flow)
        with h5py.File(source/'trajectories.h5') as h:
            for frame in [0,len(h['time_s'])//2]:
                mask=h['active'][frame];x=h['positions_m'][frame][mask];r=h['radii_m'][:][mask];ids=h['particle_id'][:][mask];t=float(h['time_s'][frame])
                f=tracker.evaluate(x,r,t);wall=tracker.cylinder.wall_gap_and_normal(x,r)
                answers={}
                for mode,components,warm,threads in [('global_cold',False,False,1),('global_warm',False,True,1),('components_cold',True,False,1),('components_warm',True,True,1),('components_parallel',True,True,workers)]:
                    model=ConstrainedLubrication(config['medium']['viscosity_Pa_s'],components=components,warm_start=warm,component_workers=threads)
                    result=model.solve(x,r,f['primary_force_N']+f['secondary_force_N'],f['flow_velocity_m_s'],wall)
                    constraints=horizon_constraints(x,r,wall,result.velocity,20e-6,config['particles']['contact_tolerance_m'],ids)
                    for repeat in range(4):
                        start=time.perf_counter();u,cache,d=model.constrained(result,constraints,ids);elapsed=time.perf_counter()-start
                        if warm:model.cache=cache
                        rows.append(dict(kind='qp',case=case,frame=frame,N=len(r),mode=mode,repeat=repeat,seconds=elapsed,diagnostics=d))
                    answers[mode]=u
                    if model.executor:model.executor.shutdown()
                for key,u in answers.items():
                    diff=float(np.max(abs(u-answers['global_cold'])));scale=max(float(np.max(abs(u))),1e-30)
                    checks.append(dict(case=case,frame=frame,comparison=key,relative_max=diff/scale,pass_gate=diff/scale<=1e-8))
                pressure=field.interpolate_complex_pressure(x)
                q=tracker.response.harmonic_response(r,field.frequency_Hz,pressure);vol=q['harmonics_m3'];vol=vol.get() if hasattr(vol,'get') else vol
                outputs={}
                for backend in ['numpy','cupy','auto']:
                    pair=AutoPairBackend(benchmark_path=config['particles']['pair_benchmark_path'],enforce_far_field=False) if backend=='auto' else ClassicalFarFieldBjerknes(backend,enforce_far_field=False)
                    for repeat in range(12):
                        start=time.perf_counter();f=pair.forces(x,r,vol,field.frequency_Hz,config['medium']['density_kg_m3']);f=f.get() if hasattr(f,'get') else f;elapsed=time.perf_counter()-start
                        rows.append(dict(kind='pair',case=case,frame=frame,N=len(r),mode=backend,repeat=repeat,seconds=elapsed,decisions=getattr(pair,'decisions',None)))
                    outputs[backend]=f
                for key,f in outputs.items():checks.append(dict(case=case,frame=frame,comparison='pair_'+key,relative_L2=float(np.linalg.norm(f-outputs['numpy'])/np.linalg.norm(outputs['numpy'])),pass_gate=bool(np.allclose(f,outputs['numpy'],rtol=2e-12,atol=1e-26))))
                write(out/'RAW_TIMING_ROWS.json',rows);write(out/'EQUIVALENCE.json',checks)
                print(case,frame,'benchmarked',flush=True)
    timing={}
    for kind in ['qp','pair']:
        timing[kind]={mode:float(np.median([r['seconds'] for r in rows if r['kind']==kind and r['mode']==mode and r['repeat']>0])) for mode in sorted(set(r['mode'] for r in rows if r['kind']==kind))}
    pair_choice=min(['numpy','cupy'],key=lambda mode:timing['pair'][mode])
    summary=dict(status='PASS' if all(c['pass_gate'] for c in checks) else 'FAIL',timings_median_s=timing,selected_pair_backend=pair_choice,selected_component_workers=workers if timing['qp']['components_parallel']<timing['qp']['components_warm'] else 1,qp_device='CPU preferred for contact QP; GPU QP not implemented or timed',source_sha256=sha(ROOT/'src/mevo_mb/bubbles/constrained_motion.py'),raw_rows=str(out/'RAW_TIMING_ROWS.json'),equivalence=str(out/'EQUIVALENCE.json'),warm_speedup=timing['qp']['components_cold']/timing['qp']['components_warm'],component_speedup=timing['qp']['global_warm']/timing['qp']['components_warm'],resource_cpu_quota=quota)
    write(out/'SUMMARY.json',summary);write(BASE/'benchmarks/BENCHMARK_SUMMARY.json',summary)
    print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
