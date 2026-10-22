"""Complete last accepted segments, display-contract checks and immutable hashes."""
import hashlib
import json
from pathlib import Path
import numpy as np
from diagnose import HERE, POP, OriginalRHS, save
from diagnostic_core import segment_cells


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    base=POP.parent
    inputs=json.loads((HERE/'INPUT_SHA256.json').read_text())
    historical=json.loads((HERE/'HISTORICAL_SOURCE_SHA.json').read_text())
    differences=[]
    for name,expected in historical['files'].items():
        path=POP/'source/particle_3d/src/particle_3d'/name
        digest=sha(path);inputs[str(path.relative_to(base))]=digest
        if digest!=expected:differences.append(name)
    assert differences==['population_v1/run.py'],differences
    rhs=OriginalRHS()
    for pid in [1189,1329,642,822]:
        result=json.loads((HERE/f'result_{pid}.json').read_text())
        with np.load(HERE/f'original_{pid}.npz') as z:
            t=z['physical_time_s'];x=z['x_m']
            spans,hits=segment_cells(rhs.env.field,x[-2],x[-1])
            segment=len(t)-2
            assert result['interfaces']['all_segment_partitions'][-1]['segment']==segment-1
            result['interfaces']['all_segment_partitions'].append(dict(segment=segment,
                cells=[v[2] for v in spans],fractions=[[v[0],v[1]] for v in spans]))
            for hit in hits:
                hit.update(segment=segment,time_s=float(t[-2]+hit['fraction']*(t[-1]-t[-2])),kind='TETRA',
                           role='Final recorded interval; complete traversal, no extra interface experiment')
                result['interfaces']['all_tetra_crossings'].append(hit)
            result['interfaces']['recorded_intervals_enumerated']=len(t)-1
            assert len(result['interfaces']['all_segment_partitions'])==len(t)-1
        save(f'result_{pid}.json',result)
    display=POP/'candidate_0/views_0to1p2s_large_axes/az240_el25p0/data/animation_samples.npz'
    render=display.parent/'render_manifest.json'
    old=POP/'analysis/realization_0000_contactfix/visualization_data.npz'
    for path in [display,render,old,POP/'scripts/render_view05_flow_axes.py']:
        inputs[str(path.relative_to(base))]=sha(path)
    curves={};qc={}
    with np.load(display) as z:
        for pid in [642,822]:
            j=int(np.flatnonzero(z['particle_id']==pid)[0])
            event=next(e for e in json.loads((HERE/'FROZEN_PLAN.json').read_text())['events'] if e['particle_id']==pid)
            ft=z['frame_time_s'];mask=(ft>=event['window_s'][0])&(ft<=event['window_s'][1])
            with np.load(POP/f'analysis/realization_0000_contactfix/tracks/mb_{pid:06d}.npz') as tr:
                t=tr['physical_time_s'];v=tr['velocity_m_s'];speed=np.linalg.norm(v[1:],axis=1)
                live=(ft>=t[0])&(ft<t[-1]);mask&=live
                integral=np.r_[0.,np.cumsum(speed*np.diff(t))]
                lo=np.maximum(t[0],ft[mask]-.05);hi=np.minimum(t[-1],ft[mask]+.05)
                means=(np.interp(hi,t,integral)-np.interp(lo,t,integral))/(hi-lo)
            stored=z['mean_speed_m_s'][mask,j]
            assert np.max(abs(means-stored))<1e-12
            shown=z['displayed'][mask,j];prior=bool(z['prior_exclusion'][j])
            assert np.array_equal(shown,(means>=.002)&~np.full(len(means),prior))
            curves[f'time_{pid}']=ft[mask];curves[f'rolling_speed_{pid}']=stored;curves[f'displayed_{pid}']=shown
            qc[str(pid)]=dict(max_mean_speed_reconstruction_error_m_s=float(np.max(abs(means-stored))),
                permanently_hidden=prior,frame_count_in_window=int(mask.sum()),threshold_mm_s=2,
                averaging_window_offsets_s=[-.05,.05])
        latest_dt=float(np.diff(z['frame_time_s']).mean())
    with np.load(old) as z:
        for e in json.loads((HERE/'FROZEN_PLAN.json').read_text())['events']:
            pid=e['particle_id'];j=int(np.flatnonzero(z['particle_id']==pid)[0]);t=z['frame_time_s'];x=z['frame_x_m'][:,j]
            valid=np.isfinite(x[:-1]).all(1)&np.isfinite(x[1:]).all(1)
            keep=valid&(t[:-1]>=e['window_s'][0])&(t[:-1]<=e['window_s'][1])
            curves[f'old_frame_time_{pid}']=t[:-1][keep]
            curves[f'old_frame_displacement_speed_{pid}']=np.linalg.norm(np.diff(x,axis=0)[keep],axis=1)/np.diff(t)[keep]
        old_dt=float(np.diff(t).mean())
    np.savez_compressed(HERE/'display_windows.npz',**curves)
    save('DISPLAY_CONTRACT_CHECK.json',dict(latest=qc,latest_physical_frame_dt_s=latest_dt,
        old_physical_frame_dt_s=old_dt,role='Read-only display postprocessing; no state or integration changes',
        future_lookahead_s=.05,mask_does_not_remove_physical_particle=True))
    # Production inputs are rehashed after all RHS/reintegration stages.
    for name,expected in inputs.items():assert sha(base/name)==expected,name
    save('INPUT_SHA256.json',inputs)
    save('FINAL_INPUT_INTEGRITY.json',dict(status='PASS',input_count=len(inputs),
        original_scientific_source_files_checked=len(historical['files']),
        changed_historical_entrypoint_not_imported=differences,scientific_source_or_original_data_modified=False))
    print('FINAL_INPUT_INTEGRITY_PASS',len(inputs),'historical_sources',len(historical['files']),flush=True)


if __name__=='__main__':main()
