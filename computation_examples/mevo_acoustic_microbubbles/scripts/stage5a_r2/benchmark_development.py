"""Development layout and synchronized pair-force benchmark, not speedup evidence."""
from time import perf_counter
import numpy as np
from scripts.stage5a_r2.run import load
from scripts.stage5a_r2.seal import BASE
from scripts.stage4l.dense_probe import write
from mevo_mb.bubbles.motion_r2 import MotionR2,anchored_constraints
from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes

def main():
    rows=[]
    for case in ['B2','B3']:
        tracker,sweep,cfg,state,inputs=load(case,'S1')
        mask=state['active'];x=state['positions_m'][mask];r=state['radii_m'][mask];ids=state['particle_id'][mask]
        t=float(state['time_s']);f=tracker.evaluate(x,r,t);wall=tracker.cylinder.wall_gap_and_normal(x,r)
        answers={}
        for mode,components,workers in [('global',False,1),('components',True,1),('parallel',True,6)]:
            m=MotionR2(cfg['medium']['viscosity_Pa_s'],components=components,component_workers=workers)
            a=m.solve(x,r,f['primary_force_N']+f['secondary_force_N'],f['flow_velocity_m_s'],wall)
            c=anchored_constraints(x,x,r,wall,wall,1e-7,ids,1e-11,a.velocity)
            for rep in range(5):
                start=perf_counter();u,cache,d=m.constrained_stage(a,c,ids,'benchmark');m.cache=cache
                rows.append(dict(case=case,kind='qp',mode=mode,repeat=rep,seconds=perf_counter()-start,
                                 kkt=d['kkt'],iterations=d['iterations'],components=d['component_sizes']))
            answers[mode]=u
            if m.executor:m.executor.shutdown()
        for mode in answers:
            np.testing.assert_allclose(answers['global'],answers[mode],atol=1e-9,rtol=1e-8)
        pressure=tracker.field.interpolate_complex_pressure(x)
        response=tracker.response.harmonic_response(r,tracker.field.frequency_Hz,pressure)
        vol=response['harmonics_m3'];vol=vol.get() if hasattr(vol,'get') else vol
        force={}
        for backend in ['numpy','cupy']:
            engine=ClassicalFarFieldBjerknes(backend,enforce_far_field=False)
            for rep in range(20):
                start=perf_counter();value=engine.forces(x,r,vol,tracker.field.frequency_Hz,cfg['medium']['density_kg_m3'])
                value=value.get() if hasattr(value,'get') else value
                rows.append(dict(case=case,kind='pair',mode=backend,repeat=rep,seconds=perf_counter()-start))
            force[backend]=value
        np.testing.assert_allclose(force['numpy'],force['cupy'],atol=1e-26,rtol=2e-12)
        write(BASE/'benchmarks/DEVELOPMENT_LAYOUT.json',dict(rows=rows,scope='development diagnostic while prefix process uses one CPU; not formal isolated performance'))
    for kind in ['qp','pair']:
        for mode in sorted({r['mode'] for r in rows if r['kind']==kind}):
            print(kind,mode,np.median([r['seconds'] for r in rows if r['kind']==kind and r['mode']==mode and r['repeat']>0]),flush=True)

if __name__=='__main__':main()
