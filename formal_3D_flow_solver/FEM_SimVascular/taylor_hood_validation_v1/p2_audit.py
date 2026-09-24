"""True P2 conservation audit; never reduces ten-node velocity to P1."""
from common import *
from p2 import *
from p2_measure import P2Measurements
from particle_3d.flowfield_conservation_diagnosis import face_topology,boundary_owners,upstream_cells,clip_polygon,tetra_plane_polygon
from scipy.spatial import ConvexHull
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import math,time,argparse

ENV=None

def weighted_quantiles(v,w,quantiles):
    v=np.asarray(v).ravel();w=np.broadcast_to(w,np.shape(v)) if np.ndim(w)==0 else np.asarray(w).ravel()
    order=np.argsort(v);v=v[order];w=w[order];return np.interp(np.array(quantiles)*w.sum(),np.cumsum(w)-w/2,v)

def divergence_nodes(points,tetra,u):
    result=np.empty((len(tetra),4))
    for i in range(0,len(tetra),10000):
        cells=tetra[i:i+10000];g=corner_gradients(points[cells[:,:4]])
        dn=gradient_at(np.eye(4)[None,:,:],g[:,None,:,:])
        result[i:i+len(cells)]=np.einsum('tqai,tai->tq',dn,u[cells])
    return result

def clipped_centroid_volume(x,center,normal):
    d=(x-center)@normal
    if np.all(d<=0):return abs(np.linalg.det(x[1:]-x[0]))/6,x.mean(axis=0)
    if np.all(d>=0):return 0.,x.mean(axis=0)
    vertices=list(x[d<=0])
    for i,j in EDGES:
        if d[i]*d[j]<0:vertices.append(x[i]+d[i]/(d[i]-d[j])*(x[j]-x[i]))
    vertices=np.asarray(vertices);scale=np.ptp(x,axis=0).max();z=(vertices-x[0])/scale;origin=z.mean(axis=0)
    hull=ConvexHull(z);tri=z[hull.simplices]
    vol=np.abs(np.einsum('ij,ij->i',tri[:,0]-origin,np.cross(tri[:,1]-origin,tri[:,2]-origin)))/6
    centroid=np.einsum('i,ij->j',vol,(tri.sum(axis=1)+origin)/4)/vol.sum()
    return float(vol.sum()*scale**3),centroid*scale+x[0]

def polygon_p2_flux(polygon,normal,corner,values):
    if len(polygon)<3:return 0.,0.
    triangles=np.stack([np.broadcast_to(polygon[0],polygon[1:-1].shape),polygon[1:-1],polygon[2:]],axis=1)
    xyz=np.einsum('qa,tai->tqi',TRI_Q,triangles)
    lam=barycentric(xyz,np.broadcast_to(corner,(*xyz.shape[:-1],4,3)))
    uq=np.einsum('tqa,ai->tqi',basis(lam),values)
    area=np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)/2
    return float(area@((uq@normal)@TRI_W)),float(area.sum())

def slab_job(row):
    points,tetra,u,divcorners,volumes,top,boundaries=ENV
    center=np.array([float(row['center_m_'+k]) for k in 'xyz']);normal=np.array([float(row['normal_'+k]) for k in 'xyz'])
    corners=tetra[:,:4];mask,d=upstream_cells(points,corners,top,center,normal,boundaries['INLET'][1])
    maxd=d[corners].max(axis=1);full=np.flatnonzero(mask&(maxd<=0));partial=np.flatnonzero(mask&(maxd>0))
    terms=list(volumes[full]*divcorners[full].mean(axis=1));volume=math.fsum(volumes[full]);qs=[];section_area=0.
    for cell in partial:
        x=points[corners[cell]];v,c=clipped_centroid_volume(x,center,normal)
        terms.append(v*float(barycentric(c,x)@divcorners[cell]));volume+=v
        polygon,_=tetra_plane_polygon(x,np.zeros((4,3)),center,normal)
        q,area=polygon_p2_flux(polygon,normal,x,u[tetra[cell]]);qs.append(q);section_area+=area
    surfaces={}
    for name,(tri,owners) in boundaries.items():
        qq=[];areas=[]
        for face,cell in zip(tri,owners):
            if not mask[cell]:continue
            x=points[face[:3]];av=np.cross(x[1]-x[0],x[2]-x[0]);n=av/np.linalg.norm(av)
            if n@(x.mean(axis=0)-points[corners[cell]].mean(axis=0))<0:n=-n
            polygon,_=clip_polygon(x,np.zeros((3,3)),center,normal)
            q,area=polygon_p2_flux(polygon,n,points[corners[cell]],u[tetra[cell]]);qq.append(q);areas.append(area)
        surfaces[name]=dict(Q_m3_s=math.fsum(qq),area_m2=math.fsum(areas))
    assert sum(v['area_m2'] for k,v in surfaces.items() if k.startswith('OUTLET'))==0
    Qin=-surfaces['INLET']['Q_m3_s'];Qw=surfaces['WALL']['Q_m3_s'];Qs=math.fsum(qs);D=math.fsum(terms)
    assert abs(Qs+Qw-Qin-D)/Qin<1e-10,'Independent P2 Gauss diagnostic failed'
    return dict(section_s_m=float(row['arclength_m']),Q_in_m3_s=Qin,Q_section_m3_s=Qs,Q_wall_m3_s=Qw,
        integral_div_m3_s=D,negative_integral_div_m3_s=-D,loss_minus_wall_m3_s=Qin-Qs-Qw,
        closure_relative=(Qs+Qw-Qin-D)/Qin,section_relative_error=Qs/Qin-1,
        clipped_volume_m3=volume,full_tetra_count=len(full),partial_tetra_count=len(partial))

def main():
    global ENV
    parser=argparse.ArgumentParser();parser.add_argument('field');args=parser.parse_args();start=time.time()
    m=P2Measurements(CASE);u,p=m.read(args.field);global_state=m.measure(u,p);points=m.points;tetra=m.tetra;corners=tetra[:,:4]
    assert global_state['velocity_finite'] and global_state['pressure_finite']
    a=baseline_arrays();geom=np.load(REPORT/'data/p1_geometry_diagnostics.npz');p1=json.loads((REPORT/'data/p1_network_bc_local_conservation_baseline.json').read_text())
    root=json.loads((OLD/'data/root_topology.json').read_text());q,w=tetra_rule15();N=basis(q);vol=m.volumes
    dc=divergence_nodes(points,tetra,u);dq=dc@q.T;weights=vol[:,None]*w
    means=dq@w;meansquares=(dq*dq)@w
    masks={'domain':np.ones(len(tetra),bool),'root':geom['root_mask'],'wall_adjacent':geom['wall_adjacent'],
           'junction':(~geom['root_mask'])&(np.linalg.norm(geom['centroids']-np.array(root['junction_m']),axis=1)<5e-6)}
    divstats={}
    for name,mask in masks.items():
        v=dq[mask];ww=weights[mask];quant=weighted_quantiles(v,ww,[.01,.05,.5,.95,.99])
        L2=math.sqrt(float(vol[mask]@meansquares[mask]))
        divstats[name]=dict(quadrature_samples=int(v.size),volume_m3=float(vol[mask].sum()),min=float(v.min()),max=float(v.max()),
            volume_weighted_mean=float(np.sum(v*ww)/ww.sum()),volume_weighted_RMS=L2/math.sqrt(vol[mask].sum()),L2_div_m_3over2_s_inv=L2,
            **dict(zip(['P01','P05','P50','P95','P99'],map(float,quant))))
    csvwrite(REPORT/'data/p2_tetra_divergence.csv',(dict(tetra_id=i,volume_m3=vol[i],element_mean_div_s_inv=means[i],element_RMS_div_s_inv=math.sqrt(meansquares[i]),
        min_quad_div_s_inv=dq[i].min(),max_quad_div_s_inv=dq[i].max()) for i in range(len(tetra))))
    np.savez_compressed(REPORT/'data/p2_geometry_diagnostics.npz',divergence_at_corners=dc,divergence_quadrature=dq,quadrature_weights=w,volumes=vol)
    shared=np.load(REPORT/'data/shared_section_geometry.npz');candidate=json.loads((REPORT/'data/shared_section_candidates.json').read_text());sections=[]
    for i,row in enumerate(candidate['all']):
        lo,hi=shared['offsets'][i:i+2];tri=shared['triangles_xyz'][lo:hi];parents=shared['parent_tetra'][lo:hi]
        xyz=np.einsum('qa,tai->tqi',TRI_Q,tri);lam=barycentric(xyz,points[corners[parents]][:,None,:,:])
        uq=np.einsum('tqa,tai->tqi',basis(lam),u[tetra[parents]]);normal=shared['normals'][i]
        area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
        Q=math.fsum(area*((uq@normal)@TRI_W));Qin=global_state['Q_in_m3_s']
        sections.append(dict(candidate_id=int(row['candidate_id']),refinement=int(row['refinement']),section_s_m=float(row['arclength_m']),
            original_geometry_pass=row['geometry_pass']=='True',Q_signed_m3_s=Q,Q_over_Qin=Q/Qin,relative_error=Q/Qin-1,area_m2=math.fsum(area)))
    csvwrite(REPORT/'data/p2_root_section_flux.csv',sections)
    top=face_topology(corners);boundaries={name:(tri,boundary_owners(top,tri[:,:3])) for name,(tri,av) in m.boundary.items()}
    ENV=(points,tetra,u,dc,vol,top,boundaries)
    with ProcessPoolExecutor(max_workers=6,mp_context=mp.get_context('fork')) as pool:slabs=list(pool.map(slab_job,candidate['selected_six']))
    csvwrite(REPORT/'data/p2_slab_gauss_balance.csv',slabs)
    difference=np.empty_like(dq);speed1=np.empty_like(dq);pressure_difference=np.empty_like(dq);pressure1=np.empty_like(dq)
    offset1=float(vol@(a['pressure_pa'][corners].mean(axis=1))/vol.sum());offset2=float(vol@(p[corners].mean(axis=1))/vol.sum())
    for i in range(0,len(tetra),10000):
        cells=tetra[i:i+10000];old=np.einsum('qa,tai->tqi',q,a['velocity_m_s'][cells[:,:4]]);new=np.einsum('qa,tai->tqi',N,u[cells])
        difference[i:i+len(cells)]=np.linalg.norm(new-old,axis=-1);speed1[i:i+len(cells)]=np.linalg.norm(old,axis=-1)
        pressure1[i:i+len(cells)]=a['pressure_pa'][cells[:,:4]]@q.T-offset1
        pressure_difference[i:i+len(cells)]=(p[cells[:,:4]]@q.T-offset2)-pressure1[i:i+len(cells)]
    vd=[]
    for name,mask in masks.items():
        rel=np.sqrt(np.sum(weights[mask]*difference[mask]**2)/np.sum(weights[mask]*speed1[mask]**2))
        quant=weighted_quantiles(difference[mask],weights[mask],[.5,.9,.99])
        vd.append(dict(region=name,relative_L2_difference=float(rel),**dict(zip(['P50_difference_m_s','P90_difference_m_s','P99_difference_m_s'],map(float,quant)))))
    csvwrite(REPORT/'data/p1_vs_p2_velocity_difference.csv',vd);dump(REPORT/'data/p1_vs_p2_velocity_difference.json',dict(regions=vd,definition='Common 15 physical quadrature points per original tetra; volume-weighted percentiles and integrated norms'))
    pd=dict(relative_pressure_pattern_L2=float(np.sqrt(np.sum(weights*pressure_difference**2)/np.sum(weights*pressure1**2))),
        gauge='Subtract each field volume-weighted mean, on the same physical domain',P1_removed_mean_pa=offset1,P2_removed_mean_pa=offset2,
        corner_difference_range_pa=[float((p[:m.ncorner]-a['pressure_pa']-offset2+offset1).min()),float((p[:m.ncorner]-a['pressure_pa']-offset2+offset1).max())])
    dump(REPORT/'data/p1_vs_p2_pressure_difference.json',pd)
    walltri,av=m.boundary['WALL'];wallnorm=av/np.linalg.norm(av,axis=1)[:,None];wv=np.einsum('qa,tai->tqi',tri_basis(TRI_Q),u[walltri])
    vn=np.einsum('tqi,ti->tq',wv,wallnorm);wall=dict(max_velocity_m_s=float(np.linalg.norm(u[m.wall_nodes],axis=1).max()),
        max_normal_velocity_m_s=float(abs(vn).max()),signed_flux_m3_s=float(np.linalg.norm(av,axis=1)@(vn@TRI_W)),
        absolute_flux_quadrature_m3_s=float(np.linalg.norm(av,axis=1)@(abs(vn)@TRI_W)),wall_nodes=len(m.wall_nodes),
        wall_midpoint_nodes=int(np.sum(m.wall_nodes>=m.ncorner)))
    legal=np.array([v['relative_error'] for v in sections if v['original_geometry_pass']]);rms=float(np.sqrt(np.mean(legal**2)));maximum=float(abs(legal).max())
    result=dict(global_flow=global_state,wall=wall,divergence=divstats,root_section_max_error=maximum,root_section_RMS_error=rms,
        section_count=len(sections),six_slabs=slabs,div_rms_improvement_factor=p1['divergence']['domain']['volume_weighted_RMS']/divstats['domain']['volume_weighted_RMS'],
        section_max_improvement_factor=p1['root_section_max_error']/maximum,section_RMS_improvement_factor=p1['root_section_RMS_error']/rms,
        velocity_difference=vd,pressure_difference=pd,elapsed_seconds=time.time()-start,source_native_VTU=str(args.field),
        source_native_VTU_sha256=sha(args.field),geometry='Exact double input, verified native output point/connectivity mapping',velocity='Native Float64 on all 10 nodes',
        quantile_definition='P2 volume-times-quadrature weighted distribution; P1 original baseline JSON retains its documented element quantiles')
    dump(REPORT/'data/p2_local_conservation_audit.json',result)
    bins=np.linspace(0,root['arclength_m'][-1],21);rows=[]
    for left,right in zip(bins[:-1],bins[1:]):
        mask=geom['root_mask']&(geom['root_arclength']>=left)&(geom['root_arclength']<right)
        if not mask.any():continue
        for label,val,weight in [('P1',geom['divergence'][mask],vol[mask]),('P2',dq[mask],weights[mask])]:
            quant=weighted_quantiles(val,weight,[.01,.1,.9,.99]);rows.append(dict(field=label,s_left_m=left,s_right_m=right,
                mean=float(np.sum(val*weight)/weight.sum()),**dict(zip(['P01','P10','P90','P99'],map(float,quant)))))
    csvwrite(REPORT/'data/root_divergence_bins.csv',rows)
    print(json.dumps({k:result[k] for k in ('root_section_max_error','root_section_RMS_error','div_rms_improvement_factor','section_max_improvement_factor','elapsed_seconds')}),flush=True)

if __name__=='__main__':main()
