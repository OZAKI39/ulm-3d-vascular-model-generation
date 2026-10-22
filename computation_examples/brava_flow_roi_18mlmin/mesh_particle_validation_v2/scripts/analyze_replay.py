"""Read every accepted state; retain held-velocity semantics and exact visit sequence."""
from pathlib import Path
import argparse,json,sys,hashlib,gzip
import numpy as np
import torch
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v.item())+'\n')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,required=True);ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();s=a.base/'mesh_particle_validation_v2';run=a.run
    summary=json.loads((run/'SUMMARY.json').read_text());identity=json.loads((run/'IDENTITY.json').read_text());name=identity['candidate']
    sys.path.insert(0,str(a.base/'population_v1/source/particle_3d/src'))
    from particle_3d.population_v1.environment import load_environment
    from particle_3d.population_v1.gpu_wall import GPUWall
    env=load_environment(identity['environment']);gpuwall=GPUWall(env.wall);audit=np.load(s/'04_field_audit'/name/'FIELD_ARRAYS.npz');eta=audit['eta'];div=audit['div'];out=run/'analysis';out.mkdir(exist_ok=False);rows=[]
    for folder in sorted(run.glob('particle_*')):
        manifest=json.loads((folder/'MANIFEST.json').read_text());parts=[]
        for c in manifest['chunks']:
            assert hashlib.sha256((folder/c['path']).read_bytes()).hexdigest()==c['file_sha256'],'Accepted-state file hash mismatch'
            z=np.load(folder/c['path']);counts=np.diff(z['state_offsets']);d={k:z[k] for k in z.files if k not in ['state_offsets','state_time_s']};d['time_s']=np.repeat(z['state_time_s'],counts);parts.append(d)
        data={k:np.concatenate([d[k] for d in parts]) for k in parts[0]};time=data['time_s'];keep=np.r_[time[1:]!=time[:-1],True];data={k:v[keep] for k,v in data.items()};time=data['time_s'];x=data['x'];ci=data['tetra_id'].copy();assert np.all(np.diff(time)>0)
        # Birth rows carry tetra=-1 by the unchanged recorder. Locate them using
        # the same authoritative predicate; a genuinely outside exit remains -1.
        for j in np.flatnonzero(ci<0):ci[j]=env.field.locate(x[j])[0]
        valid=ci>=0;fem=np.full_like(x,np.nan);jj=np.flatnonzero(valid);weights=env.field.geometry.weights(x[jj],ci[jj])
        assert np.all(env.field.geometry.contains(weights,ci[jj])),'Saved endpoint outside recorded tetra'
        fem[jj]=np.einsum('ni,nij->nj',weights,env.field.velocity_nodes_m_s[env.field.tetra[ci[jj]]])
        bubble=np.linalg.norm(data['v'],axis=1);femspeed=np.linalg.norm(fem,axis=1);gap=data['wall_gap'];dtime=np.diff(time);ds=np.linalg.norm(np.diff(x,axis=0),axis=1)
        exposure_eta=np.where(valid,eta[np.maximum(ci,0)],np.nan);exposure_div=np.where(valid,div[np.maximum(ci,0)],np.nan)
        # Event definition is fixed across meshes, using bubble speed and a
        # 50 ms lookback. It is a diagnostic definition, not an acceleration law.
        onset=[];last_event=-np.inf
        for j in np.flatnonzero(bubble>=.01):
            k=np.searchsorted(time,time[j]-.05);prior=bubble[k:j]
            if len(prior) and prior.min()<.002 and time[j]-last_event>.05:
                onset.append(dict(time_s=float(time[j]),position_m=x[j].tolist(),previous_min_m_s=float(prior.min()),speed_m_s=float(bubble[j]),tetra_id=int(ci[j])));last_event=time[j]
        prior=bubble[1:-1];rel=np.divide(abs(np.diff(bubble[1:])),prior,out=np.full_like(prior,np.nan),where=prior>0);terminal=manifest['terminal'];pid=terminal['particle_id']
        wallparts=[];contactparts=[];nearestparts=[]
        for begin in range(0,len(x),10000):
            end=min(begin+10000,len(x));wg=gpuwall.query(torch.as_tensor(x[begin:end],device='cuda'),torch.as_tensor(data['radius'][begin:end],device='cuda'))
            wallparts.append(wg[1].cpu().numpy());contactparts.append(wg[4].cpu().numpy());nearestparts.append(wg[0].cpu().numpy())
        wallgeom=np.concatenate(wallparts);wallfeature=(wallgeom[:,10:13]>128*np.finfo(float).eps).sum(1);contactcount=np.concatenate(contactparts);wall_triangle=np.concatenate(nearestparts)
        change=np.flatnonzero(wallfeature[1:]!=wallfeature[:-1])+1;feature_transitions=[dict(time_s=float(time[k]),position_m=x[k].tolist(),from_feature=int(wallfeature[k-1]),to_feature=int(wallfeature[k]),gap_m=float(gap[k]),wall_triangle=int(wall_triangle[k])) for k in change]
        assert np.max(abs(wallgeom[:,0]-gap))<1e-10,'Recorded/recomputed accepted endpoint wall gap disagree'
        transitions=[];previous=None
        with gzip.open(folder/'interaction_audit.jsonl.gz','rt') as f:
            for line in f:
                ev=json.loads(line)
                if ev.get('event')=='REJECTED_TRIAL':continue
                # Contact topology is available in original CPU component logs;
                # retain those records rather than invent a FACE/EDGE label from speed.
                for comp in ev.get('cpu_components',[]) or []:
                    if not isinstance(comp,dict):continue
                    ids=comp.get('contact_ids');tag=json.dumps(ids,sort_keys=True)
                    if ids and tag!=previous:transitions.append(dict(time_s=ev['physical_time_s'],contact_ids=ids));previous=tag
        row=dict(particle_id=pid,terminal=terminal,accepted_state_count=len(time),residence_s=terminal['residence_time_s'],path_length_m=float(ds.sum()),endpoint_m=x[-1],
           max_bubble_speed_m_s=float(bubble.max()),max_FEM_speed_m_s=float(np.nanmax(femspeed)),min_gap_m=float(gap.min()),max_per_step_displacement_m=float(ds.max(initial=0.)),
           max_relative_speed_change_excluding_birth=float(rel[np.isfinite(rel)].max(initial=0.)),zero_previous_speed_interval_count=int((prior==0).sum()),zero_to_positive_speed_count=int(((prior==0)&(bubble[2:]>0)).sum()),outside_FEM_endpoint_count=int((~valid).sum()),
           high_eta_ge_0p5_duration_s=float(dtime[exposure_eta[:-1]>=.5].sum()),high_eta_ge_0p1_duration_s=float(dtime[exposure_eta[:-1]>=.1].sum()),
           integrated_abs_div=float(np.nansum(dtime*abs(exposure_div[:-1]))),slow_fast_events=onset,contact_transitions=transitions,wall_feature_transitions=feature_transitions,
           wall_feature_definition='1 vertex, 2 edge, 3 face; nearest finite triangle barycentrics with unchanged CUDA 128*eps predicate; geometry feature alone is not active contact',
           velocity_semantics='v[i] held on (t[i-1],t[i]); birth velocity is initial free flow; field sampled separately',exposure_semantics='Accepted interval assigned to start cell; not exact within-step cell residence')
        dump(out/('particle_%04d.json'%pid),row);rows.append(row)
        np.savez_compressed(out/('particle_%04d.npz'%pid),time_s=time,x=x,v=data['v'],q=data['q'],fem_velocity=fem,fem_speed=femspeed,bubble_speed=bubble,gap=gap,tetra_id=ci,eta=exposure_eta,div=exposure_div,
            visited_cell_sequence=ci[np.r_[True,ci[1:]!=ci[:-1]]],wall_feature=wallfeature,contact_count=contactcount,wall_triangle=wall_triangle)
    dump(out/'SUMMARY.json',dict(candidate=name,dt_s=identity['dt_s'],counts=summary['counts'],trajectories=rows,run_identity_sha256=hashlib.sha256((run/'IDENTITY.json').read_bytes()).hexdigest(),
        analysis_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),accepted_state_file_hashes_verified=True,
        slow_fast_definition=dict(lookback_s=.05,slow_below_m_s=.002,fast_at_least_m_s=.01,refractory_s=.05,diagnostic_screen_not_physical_classification=True),
        historical_comparison_caveat='Historical population has MB interactions and a shared event/refinement clock. New M0 uses independent singleton clocks; differences from historical outcomes cannot be attributed only to MB interaction or mesh.'))
    print('REPLAY_ANALYZED',name,len(rows))
if __name__=='__main__':main()
