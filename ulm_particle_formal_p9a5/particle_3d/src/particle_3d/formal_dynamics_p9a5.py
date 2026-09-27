"""P9-A.5 scheduling/observation adapter around the unchanged P9-A.1 integrator.

The same stepper lives across 3/6/12 s. A private iterator extends its time budget;
The user-selected nominal dt is 1.0 ms. Force, contact, refinement, state update
and inlet-generation algorithms remain unchanged.
"""
from pathlib import Path
from dataclasses import replace
from collections import deque
from copy import deepcopy
import gzip,hashlib,json,math,os,resource,socket,time,traceback
import numpy as np
from .formal_cohort_p9a5 import (DT,HORIZONS,REL,canonical,content_sha,digest,write_new,
                                HorizonSchedule,completion_matches,terminal_status)

ENV = None
IDENTITY = None


def safe(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def nonfinite(value):
    if isinstance(value,dict):return sum(nonfinite(v) for v in value.values())
    if isinstance(value,(list,tuple)):return sum(nonfinite(v) for v in value)
    return int(not math.isfinite(value)) if isinstance(value,(float,np.floating)) else 0


def leaves(c):
    if c['proof']=='CONTINUOUS_ORIGINAL_SUPPORT_PLANE_MINUS_H_LOWER':return [c]
    if c['proof']!='UNION_OF_ORIGINAL_CONTINUOUS_HANDOFF_CERTIFICATES':raise ValueError('Unknown continuous certificate')
    parts=c['subcertificates']
    if not (c['proof_partition_only'] and c['held_velocity_path_unchanged'] and parts and parts[0]['t0_fraction']==0 and parts[-1]['t1_fraction']==1):
        raise ValueError('Invalid continuous certificate partition')
    if any(a['t1_fraction']!=b['t0_fraction'] for a,b in zip(parts,parts[1:])) or any(s['t0_fraction']>=s['t1_fraction'] for s in parts):
        raise ValueError('Noncontiguous continuous certificate partition')
    return [v for s in parts for v in leaves(s['certificate'])]


class AuditSink:
    """Stream full evidence to gzip; retain only the final support window in RAM."""
    def __init__(self,path):
        self.raw=Path(path).open('xb');self.stream=gzip.GzipFile(fileobj=self.raw,mode='wb',filename='',mtime=0,compresslevel=3)
        self.last=deque(maxlen=32);self.recent=deque(maxlen=5);self.accepted=0;self.rejected=0
        self.certificates=0;self.certificate_violations=0;self.nonfinite=0;self.handoffs=0
        self.contacts=0;self.contact_duration=0.;self.previous_contact=False;self.max_contacts=0
        self.rejections={};self.science=hashlib.sha256()
    def append(self,row):
        data=(json.dumps(row,default=safe,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
        self.stream.write(data);self.science.update(data);self.recent.append(row)
        if not row['accepted']:
            self.rejected+=1;reason=str(row.get('error'));self.rejections[reason]=self.rejections.get(reason,0)+1;return
        self.accepted+=1;self.nonfinite+=nonfinite(row)
        self.handoffs+=int(row.get('handoff_event') is not None)
        for c in row['continuous_certificates']:
            for leaf in leaves(c):
                self.certificates+=1;self.certificate_violations+=int(leaf['minimum_g_nf_bound_m'] < -leaf['roundoff_m'])
        if 'solver' in row:
            self.last.append(row);count=row['solver'].get('contact_count',0);on=count>0
            self.contacts+=int(on and not self.previous_contact);self.previous_contact=on
            self.contact_duration+=row['accepted_dt_s'] if on else 0.
            self.max_contacts=max(self.max_contacts,count)
    def close(self):
        self.stream.close();self.raw.flush();os.fsync(self.raw.fileno());self.raw.close()
    def supported(self):
        return len(self.last)==32 and all(r['solver'].get('contact_redundancy',{}).get('rank_after')==3
            and np.max(np.abs(np.asarray(r['velocity'][:3])))<=r['solver'].get('contact_kkt',{}).get('velocity_budget_m_s',0)
            and np.all(np.asarray(r['solver'].get('multipliers',[-1]))>=0) for r in self.last)


def array_sha(a):
    a=np.ascontiguousarray(a,dtype='<f8')
    return hashlib.sha256(canonical(dict(shape=list(a.shape),dtype='<f8'))+a.tobytes()).hexdigest()


def audit_result(event,meta,samples,sink,schedule,folder):
    from .particle82a_geometry import InletGeometryAudit,lower_gap
    from .validation_boundary import ValidationBoundaryClassifier
    tol=ENV.wall.roundoff_m;radius=event['radius_m'];n=len(samples)
    if not n:raise ValueError('Missing initial and terminal saved states')
    gaps=np.asarray([ENV.wall.nearest_center_triangle(p)[1]-radius for p in samples[:,1:4]])
    normal=InletGeometryAudit(ENV).normal;inlet=ValidationBoundaryClassifier({'INLET':ENV.boundaries['INLET']})
    escapes=[]
    for i,(x,y) in enumerate(zip(samples[:-1,1:4],samples[1:,1:4])):
        if (y-x)@normal<0 and inlet.first_event(x,y) is not None:escapes.append(i)
    supported=bool(sink.supported());status=terminal_status(meta,supported,HORIZONS[-1])
    outlet='O'+str(int(meta['exit_outlet'].split('_')[-1])) if meta.get('completed') else None
    hit=ENV.classifier.first_event(samples[-2,1:4],samples[-1,1:4]) if meta.get('completed') and n>=2 else None
    classification_ok=not meta.get('completed') or hit is not None and hit.role==meta['exit_outlet']
    dt=np.diff(samples[:,0]);speed=np.linalg.norm(samples[:,4:7],axis=1)
    for p in sorted((folder/'horizons').glob('*.npz')) if (folder/'horizons').exists() else []:
        prior=np.load(p)['samples']
        if not np.array_equal(prior,samples[:len(prior)]):raise ValueError('Extended trajectory changed its old prefix')
    metrics=dict(particle_id=event['particle_id'],source_event_id=event['source_event_id'],diameter_um=event['diameter_um'],
        radius_m=radius,birth_time_s=event['birth_time_s'],status=status,outlet=outlet,end_reason=meta['end_reason'],
        failure_detail=meta.get('failure_detail'),trajectory_age_s=float(samples[-1,0]),transit_time_s=meta.get('residence_time_s'),
        residence_time_s=float(samples[-1,0]),path_length_m=float(np.linalg.norm(np.diff(samples[:,1:4],axis=0),axis=1).sum()),
        minimum_wall_gap_m=float(gaps.min()),minimum_saved_wall_gap_m=float(samples[:,14].min()),minimum_g_nf_m=float(samples[:,15].min()),
        nearwall_exposure_s=float(dt[gaps[:-1]/radius<=.1].sum()),nearwall_definition='real_wall_gap/radius <= 0.1 at preceding accepted state',
        contact_count=sink.contacts,contact_duration_s=sink.contact_duration,handoff_count=sink.handoffs,max_simultaneous_contacts=sink.max_contacts,
        maximum_speed_m_s=float(speed.max()),representative_speed_m_s=float(np.median(speed)),initial_position_m=samples[0,1:4].tolist(),
        final_position_m=samples[-1,1:4].tolist(),initial_q=samples[0,10:14].tolist(),sample_count=n,accepted_steps=meta.get('accepted_steps',0),
        rejected_trials=meta.get('rejected_trials',0),provider_calls=meta.get('provider_calls',0),maximum_horizon_used_s=schedule.used,
        horizon_extensions=schedule.extensions,horizon_reached_count=len(schedule.extensions)+int(status=='LONG_RESIDENCE_CENSORED'),
        stationary_supported=supported and status=='SUPPORTED_STATIONARY',
        penetration_count=int((gaps < -tol).sum())+int((samples[:,14]<-tol).sum()),
        handoff_violation_count=int((gaps-lower_gap(radius)<-tol).sum())+int((samples[:,15]<-tol).sum())+sink.certificate_violations,
        inlet_escape_count=len(escapes),nan_inf_count=int((~np.isfinite(samples)).sum())+sink.nonfinite,
        unclassified_corruption_count=int(not classification_ok),continuous_certificates=sink.certificates,
        trajectory_file_sha256=digest(folder/'trajectory.npz'),array_scientific_sha256=array_sha(samples),
        audit_scientific_sha256=sink.science.hexdigest(),wall_seconds=meta['wall_seconds'])
    support=dict(particle_id=event['particle_id'],status=status,contact_support_verified=supported,
        last_32_accepted_contact_solves=list(sink.last),last_observed_trials=list(sink.recent),
        rejected_trial_reasons=sink.rejections,outward_inlet_segments=escapes,terminal_outlet_reclassified=classification_ok,
        last_valid_state=samples[-1].tolist(),scope='Observed unchanged P9-A.1 solver and original real-wall geometry; support is not a physiological trapping claim')
    write_new(folder/'support.json',json.loads(json.dumps(support,default=safe,allow_nan=False)))
    write_new(folder/'metrics.json',metrics)
    return metrics


def job(spec):
    event,folder=spec;folder=Path(folder)
    if completion_matches(folder,IDENTITY,event):
        row=json.loads((folder/'metrics.json').read_text());return dict(row,reused=True)
    if folder.exists():
        archive=folder.parent/'interrupted_or_mismatched';archive.mkdir(exist_ok=True)
        os.rename(folder,archive/(folder.name+'_'+str(time.time_ns())))
    folder.mkdir(parents=True)
    from .particle9a_motion import Particle9AStepper
    from .particle9a1_audit import instrument
    from .particle82a_integration import integrate_admitted
    from .particle6_stepper import bind_query_dependency
    created=[];sink=AuditSink(folder/'audit.jsonl.gz');start=time.time();cpu=time.process_time()
    def checkpoint(old,new):
        stepper=created[0];samples=np.asarray(stepper.samples,dtype=float)
        if abs(samples[-1,0]-old)>1e-12 or stepper.boundary_event!='ACTIVE':raise ValueError('Invalid horizon extension state')
        h=folder/'horizons';h.mkdir(exist_ok=True)
        p=h/f'age_{old:g}s.npz';np.savez_compressed(p,samples=samples)
        write_new(h/f'age_{old:g}s.json',dict(status='PHYSICAL_RESIDENCE_HORIZON_REACHED',old_horizon_s=old,new_horizon_s=new,
            sample_count=len(samples),prefix_array_sha256=array_sha(samples),file_sha256=digest(p),same_live_stepper=True,last_valid_state=samples[-1].tolist()))
    schedule=HorizonSchedule(checkpoint=checkpoint)
    def factory(particles,*args,**kwargs):
        p=particles[0];sample=ENV.field.sample(p.position)
        if not sample.inside_lumen:raise ValueError('Birth is outside NEW flow')
        p=replace(p,velocity=sample.velocity_m_s.copy(),omega=(.5*sample.vorticity_s_inv).copy())
        stepper=Particle9AStepper([p],*args,gradient_provider=lambda x:ENV.field.sample(x).velocity_gradient_s_inv,**kwargs)
        created.append(stepper);return instrument(stepper,enabled=True,rows=sink)
    def iteration(*args):
        if args!=(1,round(HORIZONS[-1]/DT)+1):raise ValueError('Unexpected frozen integration loop')
        return schedule.steps()
    # The integrator's exact loop/state guards are retained. Only the authorized
    # horizon and proportional total computation allowance are supplied privately.
    integration=bind_query_dependency(integrate_admitted,dict(environment=lambda:ENV,SavedTrajectoryStepper=factory,
        HORIZON=HORIZONS[-1],MAX_PROVIDER_CALLS=round(16000*HORIZONS[-1]/1.5),range=iteration))
    try:
        before=content_sha(event);meta=integration(deepcopy(event),dt=DT,output=folder/'raw')
        if content_sha(event)!=before:raise ValueError('Frozen event mutated')
        sink.close()
        raw=folder/'raw/trajectories'/f"mb_{event['particle_id']:06d}"
        os.replace(raw.with_suffix('.npz'),folder/'trajectory.npz')
        samples=np.load(folder/'trajectory.npz')['samples']
        meta.update(model='UNCHANGED_P9A1_P65',formal_identity=IDENTITY,samples_path='trajectory.npz',
                    initial_horizon_s=HORIZONS[0],maximum_horizon_used_s=schedule.used,horizon_extensions=schedule.extensions,
                    same_live_stepper_across_horizons=True,planar_statistics=created[0].planar_statistics)
        write_new(folder/'trajectory.json',meta)
        raw.with_suffix('.json').unlink();raw.parent.rmdir();raw.parent.parent.rmdir()
        metrics=audit_result(event,meta,samples,sink,schedule,folder)
        scientific=content_sha({k:v for k,v in metrics.items() if k not in {'wall_seconds','trajectory_file_sha256'}})
        receipt=dict(hostname=socket.gethostname(),pid=os.getpid(),started_unix_s=start,ended_unix_s=time.time(),
                     cpu_seconds=time.process_time()-cpu,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                     scientific_result_sha256=scientific)
        write_new(folder/'receipt.json',receipt)
        files={str(p.relative_to(folder)):digest(p) for p in sorted(folder.rglob('*')) if p.is_file()}
        write_new(folder/'COMPLETE.json',dict(identity=IDENTITY,event_sha256=content_sha(event),files=files,
                  scientific_result_sha256=scientific,terminal_status=metrics['status']))
        return dict(metrics,reused=False)
    except BaseException:
        try:sink.close()
        except (OSError,ValueError):pass
        (folder/'error.log').write_text(traceback.format_exc())
        if created:
            np.savez_compressed(folder/'last_valid_samples.npz',samples=np.asarray(created[0].samples,dtype=float))
        write_new(folder/'FAILED_EXECUTION.json',dict(particle_id=event['particle_id'],event_sha256=content_sha(event),identity=IDENTITY,
                  status='SOLVER_FAILURE',error_log_sha256=digest(folder/'error.log'),no_completed_marker=True))
        raise
