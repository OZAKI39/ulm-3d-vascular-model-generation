"""Independent reconstruction of every accepted surface constraint; C++ normals are observations only."""
from pathlib import Path
import csv,collections,numpy as np
from surface_reference import SurfaceReference,project_velocity
S=Path(__file__).resolve().parents[1];_topology=None

def check_wall_records(D,rows,g):
    global _topology
    if _topology is None:_topology=SurfaceReference.production(S)
    ref=_topology;cfg=dict(line.split() for line in (D/'case.cfg').read_text().splitlines() if line and not line.startswith('#'));factor=float(cfg.get('tie_factor',1));theta=float(cfg.get('smooth_dihedral_threshold_deg',15))
    wall=rows(D/'wall_constraint_timeseries.csv');cor=rows(D/'GEOMETRIC_PROJECTION_EVENTS.csv');used={};offset={};maxnormal=0.;maxvelocity=0.;active=0;stats={};hist=collections.Counter();clusterhist=collections.Counter();multi=0
    with (D/'WALL_SURFACE_CONTACT_EVENTS.csv').open() as f:
        contact={ (int(a['step']),int(a['stage']),int(a['particle_id'])):a for a in csv.DictReader(f) if a['accepted']=='1'}
    assert len(contact)==len(wall)
    for r in wall:
        f=lambda k:float(r[k]);key=(int(r['step']),int(r['stage']),int(r['particle_id']));c=contact[key];x=np.array([float(c['query_'+k]) for k in 'xyz']);q=ref.query(x,factor,theta)
        assert len(q['candidates'])==int(c['candidate_triangle_count']) and len(q['clusters'])==int(c['surface_cluster_count']),(key,'CLUSTER_COUNT')
        assert q['status']==c['surface_mode'],(key,q['status'],c['surface_mode'])
        assert [a['id'] for a in q['candidates']]==list(map(int,c['candidate_ids'].split('|')))
        assert [a['id'] for a in q['clusters']]==list(map(int,c['cluster_ids'].split('|')))
        N=np.array([a['normal'] for a in q['clusters']]);engineN=np.array([list(map(float,a.split(';'))) for a in c['normal_set'].split('|')]);ne=float(np.max(np.abs(N-engineN)));maxnormal=max(maxnormal,ne);assert ne<2e-8,(key,'PSEUDONORMAL',ne)
        observed_mesh_normal=np.array([float(c['pseudonormal_'+k]) for k in 'xyz']);assert np.linalg.norm(observed_mesh_normal-q['clusters'][0]['surface_normal'])<2e-8
        raw=np.array([f('v'+k+'_raw') for k in 'xyz']);u=np.array([f('v'+k+'_used') for k in 'xyz']);on=int(r['wall_constraint_active']);rawsafe=int(c['raw_safe'])
        if rawsafe:expected=raw;assert c['constraint_mode']=='NONE'
        else:
            objective,expected,subset,lam=project_velocity(raw,N);assert np.min(N@u)>=-2e-11*max(np.linalg.norm(raw),1e-12)
            assert abs(.5*np.sum((u-raw)**2)-objective)<2e-18*max(np.linalg.norm(raw),1e-12)
            logged=list(map(int,c['active_indices'].split('|'))) if c['active_indices'] else []
            for j in logged:assert abs(N[j]@u)<2e-11*max(np.linalg.norm(raw),1e-12)
        ue=np.linalg.norm(u-expected);maxvelocity=max(maxvelocity,float(ue));assert ue<2e-11*max(np.linalg.norm(raw),1e-12),(key,'QP_MINIMALITY',ue)
        assert int(c['safe'])==1 and rawsafe==int(r['raw_segment_safe']);assert all(r['omega_raw_'+k]==r['omega_used_'+k] for k in 'xyz'),'FAIL_ROTATION_MODIFIED'
        assert on==int(np.linalg.norm(u-raw)>0);assert key not in used;used[key]=u;hist[c['constraint_mode']]+=1;clusterhist[len(N)]+=1;multi+=int(len(N)>1 and not rawsafe)
        stat=stats.setdefault(key[2],dict(particle_id=key[2],wall_constraint_event_count=0,total_time_under_constraint=0.,maximum_removed_normal_speed=0.,removed_speed_sum=0.,total_geometric_projection_distance=0.,maximum_single_projection_distance=0.,wall_stall_count=0))
        if on:active+=1;stat['wall_constraint_event_count']+=1;stat['total_time_under_constraint']+=float(c['dt']);stat['maximum_removed_normal_speed']=max(stat['maximum_removed_normal_speed'],float(np.linalg.norm(u-raw)));stat['removed_speed_sum']+=float(np.linalg.norm(u-raw))
    for r in cor:
        f=lambda k:float(r[k]);before=np.array([f('before_'+k) for k in 'xyz']);after=np.array([f('after_'+k) for k in 'xyz']);normal=np.array([f('n'+k) for k in 'xyz']);distance=f('correction_distance');q=ref.query(before,factor,theta)
        assert len(q['clusters'])==1,'MULTISURFACE_POSITION_CORRECTION_FORBIDDEN';n=before-q['clusters'][0]['point'];n/=np.linalg.norm(n)
        assert 0<distance<=2.5e-11 and np.linalg.norm(normal-n)<2e-8
        assert np.linalg.norm(after-before-distance*n)<4e-20
        if int(r['accepted']):
            key=(int(r['step']),int(r['stage']),int(r['particle_id']));assert key not in offset;offset[key]=after-before
            assert contact[key]['surface_cluster_count']=='1'
            stat=stats[key[2]];stat['total_geometric_projection_distance']+=distance;stat['maximum_single_projection_distance']=max(stat['maximum_single_projection_distance'],distance)
    for s in stats.values():s['mean_removed_normal_speed']=s.pop('removed_speed_sum')/max(s['wall_constraint_event_count'],1)
    return used,offset,dict(wall_stage_records=len(wall),active_velocity_projections=active,position_corrections_all_attempts=len(cor),position_corrections_accepted=len(offset),max_independent_normal_error=maxnormal,max_independent_velocity_error_m_s=maxvelocity,true_multi_surface_active=multi,constraint_modes=dict(hist),surface_cluster_counts=dict(clusterhist),particles=list(stats.values()))
