"""Streaming mixed-population validation with exact common plug translation.

The short open control section is a bookkeeping stress fixture, not a vessel
passage or suspension-physics result. Every original finite shape is retained.
A shared constant velocity preserves all pair separations between births.
"""
from pathlib import Path
from dataclasses import replace
from collections import Counter
import csv,gzip,math,time
import numpy as np
from .particle7_cases import channel,SONOVUE
from .injection_population import PopulationSource,InjectionScheduler,LinearProfile,C_MB
from .injection_admission import FiniteSizeAdmission,event_shape
from .particle3_cases import write_json,write_rows


def run_long(folder,mb_target=1000,*,seed=2026092107):
    folder=Path(folder); folder.mkdir(parents=True,exist_ok=True)
    width=200e-6; length=1e-8; sampler,wall,classifier,velocity=channel(width,length)
    horizon=(mb_target+.25)/(C_MB*sampler.Q_m3_s)
    source=PopulationSource(SONOVUE,seed); scheduler=InjectionScheduler(source,LinearProfile([0],[sampler.Q_m3_s]))
    admission=FiniteSizeAdmission(sampler,source,velocity=velocity)
    # Exact support planes below certify the WHOLE translated path against all
    # synthetic side-wall triangles, including shapes straddling both open caps.
    transit=length/velocity[2]; active={}; counts=Counter(); exit_counts=Counter(); timeline=[]; rejections=Counter()
    started=time.monotonic(); previous_id=0; samples=[]; min_wall_bound=math.inf; candidate_count=0; narrow_pair_count=0
    columns=['particle_id','species','scheduled_time_s','admitted_time_s','exit_time_s','final_status','outlet_role','volume_m3',
        'D_um','V_fL','a_m','b_m','c_m','r','jeffery_lambda','mb_diameter_um','qw','qx','qy','qz','x_m','y_m','z_m','shape_mode',
        'sample_seed','sample_sequence','candidate_id','mb_uniform','admission_draws','wall_path_lower_bound_m']
    ledger=folder/'11_long_lifecycle.csv.gz'
    with gzip.open(ledger,'wt',newline='',compresslevel=1) as stream:
        writer=csv.DictWriter(stream,fieldnames=columns); writer.writeheader()
        while scheduler.next_time()<=horizon:
            e=scheduler.pop(); t=e['scheduled_time_s']; tag=e['particle_id']
            if tag!=previous_id+1: raise AssertionError('Nonmonotonic or reused stable ID')
            previous_id=tag
            for i in list(active):
                if active[i]['exit_time']<=t:
                    exit_counts[active[i]['species']]+=1; del active[i]
            # Center-distance broad phase in the exact common-translation frame.
            rng=source.rng[e['species']+'_POSITION']; p=None; last=None
            for draw in range(512):
                candidate_count+=1; point,tid=sampler.sample(rng); point=point[0]
                candidate_shape=event_shape(e,point)
                from .particle_shapes import Sphere
                support=np.array([candidate_shape.support(n)[j]-point[j] for j,n in enumerate(np.eye(3)[:2])])
                lower=admission.policy.lower_handoff_gap(candidate_shape.radius_m)['h_lower_m'] if isinstance(candidate_shape,Sphere) else 0.
                wall_margin=float(np.min(width/2-np.abs(point[:2])-support))
                if wall_margin<lower+wall.roundoff_m:
                    rejections['FULL_PLUG_WALL_PATH_UNCERTIFIED']+=1
                    continue
                radius=candidate_shape.bounding_radius_m
                near={}
                for i,old in active.items():
                    x=old['shape'].center_m+velocity*(t-old['birth_time'])
                    if np.linalg.norm(point-x)<=radius+old['shape'].bounding_radius_m+1e-8:
                        near[i]=old['shape'].moved(x); narrow_pair_count+=1
                p,status,detail=admission.check(e,point,near)
                detail['wall_gap_lower_bound_m']=wall_margin
                if p is not None: break
                rejections[status]+=1
                last=dict(particle_id=tag,status=status,position_m=point.tolist(),**detail)
            if p is None: raise AssertionError(dict(status='LONG_VALIDATION_BACKLOG',event=e,last_candidate=last))
            # Original P1/P3 triangle classifier supplies the first open-outlet crossing.
            hit=classifier.first_event(p.position,p.position+2*length*np.array([0.,0.,1.]))
            if hit is None or hit.role!='OUTLET_01': raise AssertionError('Outlet intersection missing')
            lifetime=hit.segment_fraction*2*transit; exit_time=t+lifetime
            active[tag]=dict(shape=p.shape(),birth_time=t,exit_time=exit_time,species=e['species'],volume_m3=e['volume_m3'])
            counts[e['species']]+=1
            if 'wall_gap_lower_bound_m' in detail: min_wall_bound=min(min_wall_bound,detail['wall_gap_lower_bound_m'])
            g=e.get('geometry'); prov=e['provenance']; row=dict(particle_id=tag,species=e['species'],scheduled_time_s=t,admitted_time_s=t,
                exit_time_s=exit_time if exit_time<=horizon else '',final_status='EXITED' if exit_time<=horizon else 'ACTIVE',
                outlet_role=hit.role if exit_time<=horizon else '',volume_m3=e['volume_m3'],D_um=prov.get('D_um',''),V_fL=prov.get('V_fL',''),
                a_m=g['a_m'] if g else '',b_m=g['b_m'] if g else '',c_m=g['c_m'] if g else '',
                r=g['c_m']/g['a_m'] if g else '',jeffery_lambda=((g['c_m']/g['a_m'])**2-1)/((g['c_m']/g['a_m'])**2+1) if g else '',
                mb_diameter_um=e.get('diameter_um',''),qw=e['q'][0],qx=e['q'][1],qy=e['q'][2],qz=e['q'][3],
                x_m=p.position[0],y_m=p.position[1],z_m=p.position[2],shape_mode=p.shape().mode,
                sample_seed=prov.get('seed',''),sample_sequence=prov['sequence'],candidate_id=prov.get('candidate_id',''),mb_uniform=prov.get('uniform',''),admission_draws=draw+1,wall_path_lower_bound_m=wall_margin)
            writer.writerow(row)
            if tag<=5000 or e['species']=='MB': samples.append(row)
            if tag%1000==0:
                timeline.append(_row(scheduler,t,active,counts,width**2*length))
            if tag%10000==0: print('long',tag,'MB',counts['MB'],'elapsed_s',round(time.monotonic()-started,1),flush=True)
        for i in list(active):
            if active[i]['exit_time']<=horizon: exit_counts[active[i]['species']]+=1; del active[i]
        final=_row(scheduler,horizon,active,counts,width**2*length); timeline.append(final)
    residual=scheduler.rbc_clock.cumulative(horizon)-scheduler.rbc_volume
    if not 0<=scheduler.mb_clock.cumulative(horizon)-scheduler.mb_count<1: raise AssertionError('MB accumulator error')
    if not -1e-25<=residual<scheduler.next_rbc['volume_m3']: raise AssertionError('RBC residual error')
    result=dict(status='PASS',role='VALIDATION_POPULATION_COMMON_PLUG_TRANSLATION',population_parameters_unchanged=True,
        Q_m3_s=sampler.Q_m3_s,horizon_s=horizon,width_m=width,length_m=length,velocity_m_s=velocity.tolist(),
        shape_outside_open_caps_allowed=True,section_role='SHORT_OPEN_CONTROL_SECTION_FOR_LIFECYCLE_STRESS_NOT_RBC_PASSAGE',
        geometry_sampling='ORIGINAL_P2_AND_SONOVUE_UNMODIFIED',ledger=str(ledger.name),ledger_complete=True,
        final=final,next_rbc_volume_m3=scheduler.next_rbc['volume_m3'],rbc_volume_residual_m3=residual,
        candidate_draws=candidate_count,narrow_pair_calls=narrow_pair_count,rejections=dict(rejections),
        duplicate_id_count=0,birth_overlap_count=0,birth_wall_violation_count=0,birth_nearfield_handoff_violation_count=0,
        nonoverlap_certificate='EXACT_P4_ADMISSION_PLUS_COMMON_TRANSLATION_PRESERVES_PAIR_SEPARATIONS_AND_EXACT_SUPPORT_PLANES_CERTIFY_ALL_SIDE_WALLS',
        minimum_saved_wall_gap_lower_bound_m=min_wall_bound,
        lammps_scope='SEPARATE_DYNAMIC_BINARY_RESTART_VALIDATION; THIS_STRESS_RUN_USES_STANDALONE_ACTIVE_STATE',
        elapsed_s=time.monotonic()-started)
    write_json(folder/'11_long_summary.json',result); write_rows(folder/'11_long_balance.csv',timeline)
    write_rows(folder/'12_long_distribution_display_sample.csv',samples)
    write_json(folder/'11_long_rng_final.json',source.state())
    return result


def _row(scheduler,t,active,counts,control_volume):
    out=dict(time_s=t,mb_expected=scheduler.mb_clock.cumulative(t),rbc_target_m3=scheduler.rbc_clock.cumulative(t),
        scheduled_rbc_volume_m3=scheduler.rbc_volume,admitted_rbc_volume_m3=scheduler.rbc_volume,
        pending_rbc_volume_m3=0.,pending_mb_count=0,pending_rbc_count=0,next_rbc_volume_m3=scheduler.next_rbc['volume_m3'])
    for sp in ['MB','RBC']:
        key=sp.lower(); values=[v for v in active.values() if v['species']==sp]
        out.update({f'scheduled_{key}_count':counts[sp],f'admitted_{key}_count':counts[sp],f'active_{key}_count':len(values),f'exited_{key}_count':counts[sp]-len(values)})
        if sp=='RBC':
            vol=math.fsum(v['volume_m3'] for v in values); out.update(active_rbc_volume_m3=vol,exited_rbc_volume_m3=scheduler.rbc_volume-vol,observed_tube_hct=vol/control_volume,control_volume_m3=control_volume)
    return out
