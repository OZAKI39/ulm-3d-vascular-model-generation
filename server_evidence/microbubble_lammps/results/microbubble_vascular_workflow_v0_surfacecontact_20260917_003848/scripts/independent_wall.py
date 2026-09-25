"""NumPy/VTK geometric oracle; no engine code or engine normal is trusted."""
import numpy as np,vtk
from numpy.linalg import norm

def segment_distance(g,a,b,required):
    da=g.distance(a)[0];db=g.distance(b)[0]
    if min(da,db)-norm(b-a)/2>=required:return min(da,db)-norm(b-a)/2
    ids=vtk.vtkIdList();lo=np.minimum(a,b)-required;hi=np.maximum(a,b)+required
    g.locator.FindCellsWithinBounds(np.column_stack([lo,hi]).ravel(),ids)
    ix=np.array([ids.GetId(i) for i in range(ids.GetNumberOfIds())],int)
    if not len(ix):return required
    if not hasattr(g,'_wall_triangles'):g._wall_triangles=g.g['points'][g.g['faces'][g.g['classes']<=2]]
    triangles=g._wall_triangles[ix]
    # endpoints to triangles are supplied by the independent VTK locator.
    best=min(da,db)
    u=b-a;A=u@u
    for k in range(3):
        p=triangles[:,k];q=triangles[:,(k+1)%3];v=q-p;w=a-p
        C=np.einsum('ij,ij->i',v,v);B=v@u;D=w@u;E=np.einsum('ij,ij->i',v,w)
        for x in [a,b]:
            t=np.clip(np.einsum('ij,ij->i',x-p,v)/np.maximum(C,1e-300),0,1)
            best=min(best,float(norm(x-p-t[:,None]*v,axis=1).min()))
        for x in [p,q]:
            s=np.clip((x-a)@u/max(A,1e-300),0,1)
            best=min(best,float(norm(a+s[:,None]*u-x,axis=1).min()))
        det=A*C-B*B;valid=det>0
        s=np.divide(B*E-C*D,det,out=np.zeros_like(det),where=valid);t=np.divide(A*E-B*D,det,out=np.zeros_like(det),where=valid)
        valid &= (s>=0)&(s<=1)&(t>=0)&(t<=1)
        if valid.any():best=min(best,float(norm(w[valid]+s[valid,None]*u-t[valid,None]*v[valid],axis=1).min()))
    # Segment intersection with triangle interior: independent barycentric solve.
    v=triangles[:,1]-triangles[:,0];w=triangles[:,2]-triangles[:,0];n=np.cross(v,w);den=n@u
    f=np.divide(np.einsum('ij,ij->i',triangles[:,0]-a,n),den,out=np.full(len(den),np.inf),where=den!=0)
    valid=(f>=0)&(f<=1)
    if valid.any():
        v=v[valid];w=w[valid];h=a+f[valid,None]*u-triangles[valid,0]
        vv=np.einsum('ij,ij->i',v,v);vw=np.einsum('ij,ij->i',v,w);ww=np.einsum('ij,ij->i',w,w);hv=np.einsum('ij,ij->i',h,v);hw=np.einsum('ij,ij->i',h,w);det=vv*ww-vw*vw
        c1=np.divide(ww*hv-vw*hw,det,out=np.full(len(det),np.inf),where=det>0);c2=np.divide(vv*hw-vw*hv,det,out=np.full(len(det),np.inf),where=det>0)
        if ((c1>=0)&(c2>=0)&(c1+c2<=1)).any():best=0.
    return best

def check_wall_records(D,rows,g):
    wall=rows(D/'wall_constraint_timeseries.csv');cor=rows(D/'GEOMETRIC_PROJECTION_EVENTS.csv');used={};offset={};max_t=0;max_n=0;active=0
    stats={}
    for r in wall:
        f=lambda k:float(r[k]);x=np.array([f('x_eval'),f('y_eval'),f('z_eval')]);d,p,_=g.distance(x);n=(x-p)/d
        engine_n=np.array([f('n'+k) for k in 'xyz']);assert norm(n-engine_n)<2e-7,('STOP_WALL_NORMAL_AMBIGUITY',n,engine_n)
        raw=np.array([f('v'+k+'_raw') for k in 'xyz']);u=np.array([f('v'+k+'_used') for k in 'xyz']);vn=raw@n;on=int(r['wall_constraint_active'])
        expected=raw-min(raw@engine_n,0)*engine_n if on else raw
        assert norm(u-expected)<2e-18
        delta=u-raw;terr=norm(delta-(delta@n)*n);max_t=max(max_t,terr);assert terr<2e-11*max(norm(raw),1e-12)
        if on:assert vn<1e-14 and int(r['raw_segment_safe'])==0 and u@n>=-2e-11*max(norm(raw),1e-12);active+=1
        if raw@engine_n>=0:assert np.array_equal(u,raw)
        assert all(r['omega_raw_'+k]==r['omega_used_'+k] for k in 'xyz'),'FAIL_ROTATION_UNEXPECTEDLY_MODIFIED'
        key=(int(r['step']),int(r['stage']),int(r['particle_id']));assert key not in used;used[key]=u
        stat=stats.setdefault(key[2],dict(particle_id=key[2],wall_constraint_event_count=0,total_time_under_constraint=0.,maximum_removed_normal_speed=0.,removed_speed_sum=0.,total_geometric_projection_distance=0.,maximum_single_projection_distance=0.,wall_stall_count=0))
        if on:
            stat['wall_constraint_event_count']+=1;stat['total_time_under_constraint']+=f('dt')/2;stat['maximum_removed_normal_speed']=max(stat['maximum_removed_normal_speed'],f('removed_normal_speed'));stat['removed_speed_sum']+=f('removed_normal_speed')
    for r in cor:
        f=lambda k:float(r[k]);before=np.array([f('before_'+k) for k in 'xyz']);after=np.array([f('after_'+k) for k in 'xyz']);d,p,_=g.distance(before);n=(before-p)/d;distance=f('correction_distance')
        assert 0<distance<=2.5e-11,'STOP_POSITION_PROJECTION_TOO_LARGE'
        assert norm(after-before-distance*n)<3e-20
        assert abs((g.distance(after)[0]-d)-(f('gap_after')-f('gap_before')))<2e-14
        if int(r['accepted']):
            key=(int(r['step']),int(r['stage']),int(r['particle_id']));assert key not in offset;offset[key]=after-before
            stat=stats[key[2]];stat['total_geometric_projection_distance']+=distance;stat['maximum_single_projection_distance']=max(stat['maximum_single_projection_distance'],distance)
    for s in stats.values():s['mean_removed_normal_speed']=s.pop('removed_speed_sum')/max(s['wall_constraint_event_count'],1)
    return used,offset,dict(wall_stage_records=len(wall),active_velocity_projections=active,position_corrections_all_attempts=len(cor),position_corrections_accepted=len(offset),max_tangential_vector_error_m_s=max_t,particles=list(stats.values()))
