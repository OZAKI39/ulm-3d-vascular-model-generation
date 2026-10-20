from common import *
from native import native
from work_controller import WorkController
from current_surface import CurrentSurface
import time,resource
start=time.perf_counter();cpu=time.process_time();results=[]
assert json.loads((O/'GEOMETRY_CHECKS.json').read_text())['status']=='PASS'
for mode in CASES:
    ns,C=native(O/'checks'/mode,mode);s=ns.dynamics();version=ns.l.version();ref=ns.ref
    ns.l.command(f'reset_timestep 0 time {HOLD:.17g}')
    vel=np.zeros_like(ref['X']);vel[:ns.n]=[.02,.003,0]
    ids=ns.l.numpy.extract_atom('id')[:ns.N]-1;ns.l.numpy.extract_atom('v')[:ns.N]=vel[ids]
    ns.l.command('run 0');p=ns.state();h=MIN_H
    controller=WorkController(ns.tool.field);ctrl=controller.prepare(p,0,h)
    pred=p['x'].copy();pred[:ns.n]+=h*(p['velocity'][:ns.n]+h*p['force'][:ns.n]/(2*RHO*V))
    before=ns.tool.ever_loaded.copy()
    # Trial geometry must not commit cumulative history or native damage.
    ns.tool.evaluate(pred,p['velocity'],p['d'],p['time']+h)
    assert np.array_equal(before,ns.tool.ever_loaded)
    s=ns.step(h);vv=p['velocity'][:ns.n]+h*(p['force'][:ns.n]+s['force'][:ns.n])/(2*RHO*V)
    xp_error=float(abs(pred-s['x']).max());vv_error=float(abs(vv-s['velocity'][:ns.n]).max())
    end_error=float(abs(s['force_stage_velocity']-s['velocity']).max())
    uc,fc,_=C(ref).evaluate(s['x'],s['d']);z=ns.tool.evaluate(s['x'],s['velocity'],s['d'],s['time'])
    port=float(abs(s['force']-s['material_force']-fc-z['force']).max())
    work=ns.tool.work(p,s)
    res=float(work.sum()-(s['scalars'][0]-p['scalars'][0]+s['K']-p['K']+s['scalars'][1]-p['scalars'][1]+s['Uc']-p['Uc']))
    assert xp_error<2e-18 and vv_error<1e-12 and end_error<1e-11 and port<1e-15
    assert np.array_equal(s['d'],p['d']) and np.array_equal(s['kappa']>=p['kappa'],np.ones(len(s['d']),bool))
    assert work[1]<0 and abs(res)/max(abs(work.sum()),1e-30)<.01
    assert not s['external_force'][ns.n:].any()
    np.savez_compressed(O/'checks'/mode/'one_step.npz',**s)
    results.append(dict(mode=mode,native_steps=1,scope='Single new callback/integrator step with nonzero velocity at HOLD; no fracture calibration or main trial',prediction_error_m=xp_error,VV_velocity_error=vv_error,implicit_velocity_error=end_error,force_port_error_N=port,actual_drive_work_J=work[0],actual_drag_work_J=work[1],energy_residual_J=res,relative_energy_error=abs(res)/abs(work.sum()),controller=ctrl))
    ns.close()
dump(O/'KEY_CHECKS.json',dict(status='PASS',native_version=version,cases=results,new_native_micro_steps=2,main_dynamics_runs=0,wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
print('NATIVE_NEW_PATH_PASS',results,flush=True)
