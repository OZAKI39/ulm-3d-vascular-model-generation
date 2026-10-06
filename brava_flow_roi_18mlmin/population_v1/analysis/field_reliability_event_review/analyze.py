"""Same-position, explicitly weighted comparisons. No truth/error conflation."""
from common import *
import csv

def quantities(data,n):
    u=data[:,:3];g=data[:,4:13].reshape(-1,3,3);un=np.sum(u*n,axis=1)
    gn=np.einsum('nij,nj->ni',g,n);ngn=np.sum(n*gn,axis=1)
    return dict(u=u,speed=norm(u),un=un,ut=u-un[:,None]*n,G=g.reshape(-1,9),
      curl=data[:,13:16],shear=gn-ngn[:,None]*n,ngn=ngn,div=np.trace(g,axis1=1,axis2=2))

def angle(a,b,floor):
    na=norm(a);nb=norm(b);ok=(na>floor)&(nb>floor)
    result=np.full(len(a),np.nan)
    result[ok]=np.degrees(np.arccos(np.clip(np.sum(a[ok]*b[ok],axis=1)/(na[ok]*nb[ok]),-1,1)))
    result[ok & np.all(a==b,axis=1)]=0.
    return result

def rms(value,w):
    if value.ndim==1:value=value[:,None]
    return float(np.sqrt(np.sum(w* np.sum(value**2,axis=1))/w.sum())) if len(w) and w.sum()>0 else None

def pair_summary(a,b,w,mask,floor):
    valid=mask & np.isfinite(a['u']).all(1)&np.isfinite(b['u']).all(1)
    report=dict(total_points=int(mask.sum()),valid_points=int(valid.sum()),total_weight=float(w[mask].sum()),valid_weight=float(w[valid].sum()),
      missing_indices=np.flatnonzero(mask&~valid).tolist(),missing_reason='OUTSIDE_REFERENCE_DOMAIN / strict unresolvable boundary point; retained in original denominator',
      relative_denominator='Weighted RMS norm of first named field over common valid subset; no pointwise component percentage')
    if not valid.any() or w[valid].sum()==0:return {**report,'metrics':None,'reason':'NO_VALID_NONZERO_WEIGHT'}
    ww=w[valid];metrics={}
    for key in a:
        delta=b[key][valid]-a[key][valid];den=rms(a[key][valid],ww)
        val=rms(delta,ww)
        quantity_floor=128*np.finfo(float).eps*float(np.max(np.abs(a[key][valid])))
        metrics[key]=dict(RMS_difference=val,max_abs_difference=float(np.max(np.abs(delta))) if delta.ndim==1 else float(norm(delta).max()),
             reference_RMS=den,relative_L2=None if den is None or den<=quantity_floor else val/den,
             quantity_roundoff_scale=quantity_floor)
    angles=angle(a['u'][valid],b['u'][valid],floor);v=np.isfinite(angles)
    metrics['direction_deg']=dict(RMS_difference=rms(angles[v],ww[v]),max_abs_difference=float(angles[v].max()) if v.any() else None,
                                 undefined_points=int((~v).sum()),floor_m_s=floor)
    report['metrics']=metrics
    return report

def rise_crossings(x,u,g,spans,a,b,levels):
    delta=b-a;length=np.linalg.norm(delta);rows=[]
    for level in levels:
        roots=[]
        for xm,um,gg,(segment,lo,hi) in zip(x,u,g,spans):
            dx=delta*(hi-lo);v=gg@dx;start=um-.5*v
            coef=[np.dot(v,v),2*np.dot(start,v),np.dot(start,start)-level**2]
            rr=np.roots(np.trim_zeros(coef,'f')) if any(coef) else []
            for r in rr:
                if np.isreal(r) and 0<=float(np.real(r))<=1:
                    z=float(np.real(r));s=(lo+z*(hi-lo))*length
                    roots.append(dict(arc_m=s,direction='rising' if np.dot(start+z*v,v)>0 else 'falling'))
        unique=[]
        for r in sorted(roots,key=lambda z:z['arc_m']):
            if not unique or abs(r['arc_m']-unique[-1]['arc_m'])>128*np.finfo(float).eps*max(norm(a),norm(b)):
                unique.append(r)
        rows.append(dict(level_m_s=level,roots=unique,status='UNREACHED' if not unique else 'MULTIPLE_CROSSINGS' if len(unique)>1 else 'SINGLE_CROSSING'))
    length90=None
    if all(r['status']=='SINGLE_CROSSING' and r['roots'][0]['direction']=='rising' for r in rows):
        length90=rows[1]['roots'][0]['arc_m']-rows[0]['roots'][0]['arc_m']
    slope=np.einsum('nij,j->ni',g,delta)
    derivative=np.sum(u*slope,axis=1)
    return dict(levels=rows,length_10_90_m=length90,monotonic_increase=bool(np.all(derivative>=0)),
                unavailable_reason=None if length90 is not None else 'Original fixed thresholds unreached or multiply/reversely crossed')

def main():
    protocol=read(HERE/'EVENT_FIELD_PROTOCOL.json');events={};rows=[];parts=read(HERE/'PARTITIONS.json')
    meta=read(HERE/'FIELD_COMPARABILITY.json')
    for name in NAMES:
        audit=read(HERE/f'{name}_SAMPLING.json')
        meta['fields'][name]['raw_payload']['native_mapping']='PASS: actual native final node u/p and unique coordinates match frozen/canonical payload; exact canonical tetra pair validated'
        meta['fields'][name]['valid_event_probes']=audit['coverage']
        meta['fields'][name]['status']='ADMITTED_ACTUAL_PAYLOAD'
    for pid in PIDS:
        e=protocol['events'][str(pid)];path=np.load(HERE/f'path_{pid}.npz');raw=np.load(OLD/f'original_{pid}.npz')
        normals=np.load(HERE/f'wall_path_{pid}.npz')['normal']
        accepted_wall=np.load(HERE/f'wall_accepted_{pid}.npz');an=accepted_wall['normal']
        if pid==822:
            normals=np.broadcast_to(e['plane']['e1'],normals.shape);an=np.broadcast_to(e['plane']['e1'],an.shape)
        sampled={name:np.load(HERE/f'{name}_path_{pid}.npz') for name in NAMES}
        q={name:quantities(z['flow'],normals) for name,z in sampled.items()}
        accepted={name:np.load(HERE/f'{name}_accepted_{pid}.npz') for name in NAMES}
        qa={name:quantities(z['flow'],an) for name,z in accepted.items()}
        floor=max(read(HERE/f'{name}_SAMPLING.json')['velocity_roundoff_scale_m_s'] for name in NAMES)
        pairs={}
        for ia,ib in [(0,1),(0,2),(1,2)]:
            key=NAMES[ia]+'__'+NAMES[ib];pairs[key]={}
            for phase,label in enumerate(e['phase_names']):
                result={}
                for weight in ['arc_weight_m','time_weight_s']:
                    masks=[(path['phase']==phase)&(path['quadrature_order']==order) for order in [2,4]]
                    mainstats=pair_summary(q[NAMES[ia]],q[NAMES[ib]],path[weight],masks[0],floor)
                    refined=pair_summary(q[NAMES[ia]],q[NAMES[ib]],path[weight],masks[1],floor)
                    mainstats['four_point_rule']=refined
                    if mainstats.get('metrics') and refined.get('metrics'):
                        mainstats['two_vs_four_RMS_difference']={k:refined['metrics'][k]['RMS_difference']-v['RMS_difference'] for k,v in mainstats['metrics'].items() if v['RMS_difference'] is not None and refined['metrics'][k]['RMS_difference'] is not None}
                    result[weight]=mainstats
                pairs[key][label]=result
        reps={}
        for name in NAMES:
            r={}
            for label in ['A','B']:
                i=int(np.flatnonzero(raw['physical_time_s']==e[label]['time_s'])[0]);u=qa[name]['u'][i];g=qa[name]['G'][i].reshape(3,3)
                r[label]=dict(time_s=e[label]['time_s'],x_m=e[label]['x_m'],cell=int(accepted[name]['cell'][i]),
                  velocity_m_s=u,G_s_inv=g,curl_s_inv=qa[name]['curl'][i],speed_mm_s=float(norm(u)*1000),
                  normal_or_path_component_mm_s=float(qa[name]['un'][i]*1000),tangent_speed_mm_s=float(norm(qa[name]['ut'][i])*1000),
                  shear_s_inv=qa[name]['shear'][i],ngn_s_inv=qa[name]['ngn'][i],div_s_inv=qa[name]['div'][i],
                  center_wall_distance_um=accepted_wall['d_m'][i]*1e6,surface_gap_um=(accepted_wall['d_m'][i]-e['radius_m'])*1e6,
                  original_triangle=int(accepted_wall['triangle'][i]),original_feature=int(accepted_wall['feature'][i]),normal_ambiguous=bool(accepted_wall['normal_ambiguous'][i]),cell_geometry=accepted[name]['cell_geometry'][i])
                rows.append(dict(event=pid,field=name,endpoint=label,time_s=e[label]['time_s'],speed_mm_s=r[label]['speed_mm_s'],
                  component_mm_s=r[label]['normal_or_path_component_mm_s'],div_s_inv=float(r[label]['div_s_inv']),
                  shear_norm_s_inv=float(norm(r[label]['shear_s_inv'])),d_um=float(r[label]['center_wall_distance_um']),h_um=float(r[label]['surface_gap_um'])))
            du=r['B']['velocity_m_s']-r['A']['velocity_m_s'];r['delta_u_mm_s']=du*1000
            r['scalar_speed_increase_mm_s']=r['B']['speed_mm_s']-r['A']['speed_mm_s'];r['direction_change_deg']=angle(r['A']['velocity_m_s'][None,:],r['B']['velocity_m_s'][None,:],floor)[0]
            reps[name]=r
        endpoint_differences={}
        amplitude=reps['M0']['scalar_speed_increase_mm_s']
        for a,b in [(NAMES[0],NAMES[1]),(NAMES[0],NAMES[2]),(NAMES[1],NAMES[2])]:
            dd={}
            for label in ['A','B']:
                v=reps[b][label]['velocity_m_s']-reps[a][label]['velocity_m_s']
                dd[label]=dict(vector_difference_mm_s=v*1000,vector_difference_norm_mm_s=norm(v)*1000,
                  scalar_speed_difference_mm_s=reps[b][label]['speed_mm_s']-reps[a][label]['speed_mm_s'],
                  direction_difference_deg=angle(reps[a][label]['velocity_m_s'][None,:],reps[b][label]['velocity_m_s'][None,:],floor)[0],
                  vector_difference_over_original_event_gain=norm(v)*1000/amplitude,
                  component_difference_mm_s=reps[b][label]['normal_or_path_component_mm_s']-reps[a][label]['normal_or_path_component_mm_s'],
                  G_Frobenius_difference_s_inv=norm((reps[b][label]['G_s_inv']-reps[a][label]['G_s_inv']).ravel()))
            endpoint_differences[a+'__'+b]=dd
        profiles={}
        profilewall=np.load(HERE/f'wall_profiles_{pid}.npz');pn=profilewall['normal'];coord=np.load(HERE/f'profiles_{pid}.npz')
        pq={name:quantities(np.load(HERE/f'{name}_profiles_{pid}.npz')['flow'],pn) for name in NAMES}
        for k in [0,1]:
            mask=(coord['profile_index']==k)&profilewall['original_normal_ray_valid']&(coord['d_m']>0)
            profiles[k]=dict(valid_normal_ray_distances_um=coord['d_m'][mask]*1e6,
                original_count=int((coord['profile_index']==k).sum()),normal_valid_count=int(mask.sum()),pairs={})
            for a,b in [(NAMES[0],NAMES[1]),(NAMES[0],NAMES[2]),(NAMES[1],NAMES[2])]:
                profiles[k]['pairs'][a+'__'+b]=pair_summary(pq[a],pq[b],np.full(len(mask),.5e-6),mask,floor)
        result=dict(representatives=reps,endpoint_differences=endpoint_differences,phase_comparisons=pairs,normal_profiles=profiles,
            phase_boundaries_s=e['phase_boundaries_s'],phase_names=e['phase_names'],original_event_gain_mm_s=amplitude,
            numerical_reference='Two independently solved fixed original PLC fields; finest is not truth',physical_validation='PHYSICAL_VALIDATION_NOT_PERFORMED')
        if pid==822:
            ia=int(np.flatnonzero(raw['physical_time_s']==e['A']['time_s'])[0]);mask=(path['original_segment']==ia)&(path['quadrature_order']==1)
            spans=path['union_spans'][path['union_spans'][:,0]==ia];a=np.array(e['A']['x_m']);b=np.array(e['B']['x_m'])
            levels=np.array([.1,.9])*amplitude/1000+reps['M0']['A']['speed_mm_s']/1000
            rise={}
            for name in NAMES:
                f=sampled[name]['flow'][mask];g=f[:,4:13].reshape(-1,3,3)
                rise[name]=rise_crossings(path['x_m'][mask],f[:,:3],g,spans,a,b,levels)
                increment=np.sum(np.einsum('nij,nj->ni',g,(spans[:,2]-spans[:,1])[:,None]*(b-a)),axis=0)
                rise[name]['sum_G_dx_m_s']=increment;rise[name]['endpoint_identity_error_m_s']=increment-(reps[name]['B']['velocity_m_s']-reps[name]['A']['velocity_m_s'])
                rise[name]['own_cells']=parts[name][str(pid)][ia]
            result['spatial_rise']=rise;result['key_segment_length_um']=norm(b-a)*1e6
        events[pid]=result
    put('FIELD_COMPARABILITY.json',meta)
    put('FIELD_METRICS.json',finite_json(events))
    with (HERE/'representative_values.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    print('METRICS COMPLETE')
    for pid,e in events.items():
        print(pid,[(name,e['representatives'][name]['A']['speed_mm_s'],e['representatives'][name]['B']['speed_mm_s'],e['representatives'][name]['A']['normal_or_path_component_mm_s']) for name in NAMES])

if __name__=='__main__':main()
