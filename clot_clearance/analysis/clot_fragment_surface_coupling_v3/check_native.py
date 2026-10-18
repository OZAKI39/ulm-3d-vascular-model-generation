"""Minimal new-surface native check: 2 microsteps, one restarted continuation.

Starts from the real V2 end-state restart. This is NOT a formal protocol or a
material/load search; it exercises the new active L2 coupling and sidecar.
"""
from common import *
from native import native,predict
from fragment_recording import FragmentFluidRecorder
import time

def main():
    start=time.perf_counter();a=O/'checks/continuous';ns,C=native(a,'biphasic',V2_O/'cases/biphasic/final.restart')
    ref=ns.ref;model=energy_model(ref);contact=C(ref);initial=ns.dynamics()
    for key in ['d','kappa','x','atom_ids']:
        assert np.array_equal(initial[key],load_npz(V2_O/'cases/biphasic/keys/final.npz')[key]),key
    z=ns.tool.evaluate(initial['x'],initial['velocity'],initial['d'],initial['time']);assert z['active'][ns.tool.L:].any();ns.tool.commit(z)
    accepted=ns.tool._accepted;history=(initial['kappa'].copy(),initial['d'].copy(),float(initial['scalars'][1]));ever=ns.tool.ever_loaded.copy()
    # Predictor is read-only; static force evaluations are noncommitting.
    xt=initial['x'].copy();xt[:ns.n,0]+=1e-6
    ns.tool.geometry(xt,initial['d']);ns.static(xt);ns.static(initial['x']);ns.tool.rollback_trial();p=ns.state()
    assert np.array_equal(ns.field(2),history[0]) and np.array_equal(ns.field(3),history[1]) and ns.field(8)[1]==history[2]
    assert ns.tool._accepted is accepted and np.array_equal(ever,ns.tool.ever_loaded)
    results=[];saved=None
    for j in range(2):
        h=MIN_H;xp=p['x'].copy();xp[:ns.n]+=h*(p['velocity'][:ns.n]+h*p['force'][:ns.n]/(2*RHO*V))
        kp,dp=predict(model,ref,xp,p['kappa'],ns.field(13),ns.field(14));s=ns.step(h)
        assert np.max(abs(s['x']-xp))<2e-18 and np.max(abs(s['d']-dp))<1e-12
        assert (s['d']>=p['d']).all() and (s['kappa']>=p['kappa']).all()
        z=ns.tool.evaluate(s['x'],s['velocity'],s['d'],s['time']);ns.tool.commit(z)
        recorder=FragmentFluidRecorder(ref,ns.tool);recorder.advance(p,s)
        graph=module('graph_for_recorder_check',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
        lab,sizes,attached=graph.topology(ref,s['d']);recorder.set_state(z,s['x'],s['velocity'],lab)
        records=[recorder.component(int(k),np.flatnonzero(lab==k)) for k in np.flatnonzero(~attached)]
        assert all(q['fluid_power_W']<=1e-17 for q in records)
        work=ns.tool.work(p,s);nf=(s['d']==1)&(p['d']<1);b=ref['bonds'][nf]
        J=.5*contact.c*np.sum(np.minimum(np.linalg.norm(s['x'][b[:,1]]-s['x'][b[:,0]],axis=1)-DX,0)**2)
        delta=s['scalars'][0]-p['scalars'][0]+s['K']-p['K']+s['scalars'][1]-p['scalars'][1]+s['Uc']-p['Uc']-J
        res=work.sum()-delta;rel=abs(res)/max(abs(work).sum(),abs(delta),1e-30)
        assert rel<.01 and work[1]<0 and abs(work[0])<1e-20,(res,rel)
        assert np.max(abs(z['force']-s['external_force']))<1e-15
        results.append(dict(step=s['step'],dt_s=h,relative_increment_energy_residual=rel,residual_J=res,drive_J=work[0],drag_J=work[1],L2_active_patches=int(z['active'][ns.tool.L:].sum())))
        if j==0:
            saved=s;ns.l.command(f'write_restart {a/"audit.restart"}');ns.tool.save_checkpoint(a/'surface_restart.npz',s)
        p=s
    continuous=s;on=ns.field(13);off=ns.field(14);ns.close()
    rs,_=native(O/'checks/restarted','biphasic',a/'audit.restart');s=rs.dynamics();rs.tool.restore_checkpoint(a/'surface_restart.npz',s)
    assert np.array_equal(rs.field(13),on) and np.array_equal(rs.field(14),off)
    for key in ['x','velocity','d','kappa','atom_ids']:assert np.array_equal(s[key],saved[key]),key
    r=rs.step(MIN_H);errors={key:float(np.max(abs(r[key]-continuous[key]))) for key in ['x','velocity','force','d','kappa','scalars']}
    # Instrumentation counters and wall time are intentionally not restart invariants.
    errors['physical_scalars']=float(np.max(abs(r['scalars'][:6]-continuous['scalars'][:6])));del errors['scalars']
    assert errors['x']<1e-18 and errors['velocity']<1e-12 and errors['force']<1e-14 and errors['d']<1e-12 and errors['kappa']<1e-12 and errors['physical_scalars']<1e-17,errors
    rs.close()
    dump(O/'NATIVE_SURFACE_HISTORY_RESTART_AUDIT.json',dict(status='PASS',microsteps=3,formal_protocols=0,source='real V2 Biphasic final.restart',surface_force_microsteps=results,restart_errors=errors,accepted_only_commit=True,predictor_and_static_trial_no_history_commit=True,original_ID_exact=True,wall_s=time.perf_counter()-start))
    dump(O/'KEY_CHECKS.json',dict(status='PASS',surface='SURFACE_GEOMETRY_AND_MAPPING_AUDIT.json',native='NATIVE_SURFACE_HISTORY_RESTART_AUDIT.json',formal_runs_started=0))
    print('NATIVE_SURFACE_CHECK_PASS',results,errors,flush=True)

if __name__=='__main__':main()
