"""Re-audit the exact network-pressure P1 field, without running particles."""
from common import *
import os, time, math
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import pyvista as pv
from particle_3d.flowfield_conservation_diagnosis import (p1_divergence,face_topology,boundary_owners,slab_balance,nearest_root_arclength,stats)
from particle_3d.interior_section import cut_tetrahedra,root_component,p1_flux

ENV=None

def section_job(row):
    a,grid,div,vol,top,boundaries,qin=ENV
    center=np.array([float(row['center_m_'+i]) for i in 'xyz'])
    normal=np.array([float(row['normal_'+i]) for i in 'xyz'])
    sec,ncomp,contains=root_component(cut_tetrahedra(grid,center,normal),center,1e-16)
    assert contains
    tri=sec.faces.reshape(-1,4)[:,1:]
    result=p1_flux(sec,normal)
    out=dict(candidate_id=int(row['candidate_id']),refinement=int(row['refinement']),section_s_m=float(row['arclength_m']),
        original_geometry_pass=row['geometry_pass']=='True',Q_signed_m3_s=result['signed_Q_m3_s'],Q_positive_m3_s=result['positive_Q_m3_s'],
        Q_over_Qin=result['signed_Q_m3_s']/qin,relative_error=result['signed_Q_m3_s']/qin-1,area_m2=result['area_m2'],
        root_contains_center=bool(contains),global_plane_components=int(ncomp))
    return out,np.asarray(sec.points,float)[tri],np.asarray(sec.cell_data['parent_tetra_id'],int)

def slab_job(row):
    a,grid,div,vol,top,boundaries,qin=ENV
    c=np.array([float(row['center_m_'+i]) for i in 'xyz']);n=np.array([float(row['normal_'+i]) for i in 'xyz'])
    d,w=slab_balance(a['points_m'],a['tetra'],a['velocity_m_s'],div,vol,top,boundaries,c,n)
    Qwall=d['surface_fluxes']['WALL']['Q_m3_s'];Qi=-d['surface_fluxes']['INLET']['Q_m3_s'];Qs=d['Q_section_numpy_m3_s'];D=d['volume_integral_divergence_m3_s']
    assert sum(x['area_m2'] for k,x in d['surface_fluxes'].items() if k.startswith('OUTLET'))==0
    assert abs(Qi/qin-1)<1e-10
    assert abs(Qs+Qwall-Qi-D)/qin<1e-10
    out=dict(section_s_m=float(row['arclength_m']),Q_in_m3_s=Qi,Q_section_m3_s=Qs,Q_wall_m3_s=Qwall,
        integral_div_m3_s=D,negative_integral_div_m3_s=-D,loss_minus_wall_m3_s=Qi-Qs-Qwall,
        closure_relative=(Qs+Qwall-Qi-D)/qin,section_relative_error=Qs/qin-1,clipped_volume_m3=d['clipped_volume_m3'])
    return out

def main():
    global ENV
    started=time.monotonic();a=baseline_arrays();p,t,u=a['points_m'],a['tetra'],a['velocity_m_s']
    grid=pv.read(BASE/'frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu')
    assert np.array_equal(grid['Velocity'],u) and np.array_equal(grid.points,p)
    grid.points=p.copy();grid.cell_data['parent_tetra_id']=np.arange(len(t))
    div,vol=p1_divergence(p,t,u);top=face_topology(t);boundaries={};flux={};wall={}
    centers=p[t].mean(axis=1)
    for path in (BASE/'SV_MESH/mesh-surfaces').glob('*.vtp'):
        surface=pv.read(path);ids=np.asarray(surface['GlobalNodeID'],int)-1
        tri=ids[surface.faces.reshape(-1,4)[:,1:]];owner=boundary_owners(top,tri);boundaries[path.stem]=(tri,owner)
        x=p[tri];av=np.cross(x[:,1]-x[:,0],x[:,2]-x[:,0])/2
        av[np.einsum('ij,ij->i',av,x.mean(axis=1)-centers[owner])<0]*=-1
        q=np.einsum('ij,ij->i',av,u[tri].mean(axis=1));flux[path.stem]=math.fsum(q)
        if path.stem=='WALL':
            normal=av/np.linalg.norm(av,axis=1)[:,None]
            wall=dict(max_velocity_m_s=float(np.linalg.norm(u[np.unique(tri)],axis=1).max()),
                max_normal_velocity_m_s=float(abs(np.einsum('tni,ti->tn',u[tri],normal)).max()),
                signed_flux_m3_s=math.fsum(q),total_absolute_triangle_flux_m3_s=math.fsum(abs(q)))
    qin=-flux['INLET'];outs={k:flux[k] for k in ('OUTLET_01','OUTLET_02','OUTLET_03')}
    root=json.loads((OLD/'data/root_topology.json').read_text())
    arc,dist=nearest_root_arclength(centers,np.array(root['points_m']))
    daughter=np.full(len(t),np.inf)
    for x,y in np.array(root['daughter_segments_m']):
        e=y-x;s=np.clip((centers-x)@e/(e@e),0,1)
        daughter=np.minimum(daughter,np.linalg.norm(centers-x-s[:,None]*e,axis=1))
    rootmask=(dist<=daughter)&(arc<root['arclength_m'][-1]);wallmask=np.zeros(len(t),bool);wallmask[boundaries['WALL'][1]]=True
    np.savez_compressed(REPORT/'data/p1_geometry_diagnostics.npz',divergence=div,volumes=vol,centroids=centers,
        root_arclength=arc,root_mask=rootmask,wall_adjacent=wallmask,nearest_root_distance=dist)
    with (REPORT/'data/p1_tetra_divergence.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['tetra_id','x_m','y_m','z_m','volume_m3','div_s_inv','root_arclength_m','root_region','wall_adjacent'])
        w.writerows((i,*centers[i],vol[i],div[i],arc[i],int(rootmask[i]),int(wallmask[i])) for i in range(len(t)))
    records={k:stats(div[m],vol[m]) for k,m in {'domain':np.ones(len(t),bool),'root':rootmask,'wall_adjacent':wallmask,'junction':(~rootmask)&(np.linalg.norm(centers-np.array(root['junction_m']),axis=1)<5e-6)}.items()}
    for v in records.values():
        if v['count']:v['L2_div_m_3over2_s_inv']=v['volume_weighted_RMS']*np.sqrt(v['volume_m3'])
    candidates=list(csv.DictReader((OLD/'data/interior_section_search.csv').open()))
    chosen=json.loads((OLD_DIAG/'data/selected_candidates.json').read_text())['candidates']
    ENV=(a,grid,div,vol,top,boundaries,qin)
    print('Exact P1 divergence computed; auditing',len(candidates),'sections and',len(chosen),'slabs',flush=True)
    with ProcessPoolExecutor(max_workers=6,mp_context=mp.get_context('fork')) as pool:
        sections=list(pool.map(section_job,candidates));slabs=list(pool.map(slab_job,chosen))
    csvwrite(REPORT/'data/p1_root_section_flux.csv',[r[0] for r in sections]);csvwrite(REPORT/'data/p1_slab_gauss_balance.csv',slabs)
    offsets=np.r_[0,np.cumsum([len(r[2]) for r in sections])]
    np.savez_compressed(REPORT/'data/shared_section_geometry.npz',offsets=offsets,
        triangles_xyz=np.concatenate([r[1] for r in sections]),parent_tetra=np.concatenate([r[2] for r in sections]),
        normals=np.array([[float(r['normal_'+k]) for k in 'xyz'] for r in candidates]),
        centers=np.array([[float(r['center_m_'+k]) for k in 'xyz'] for r in candidates]))
    dump(REPORT/'data/shared_section_candidates.json',dict(all=candidates,selected_six=chosen,source_sha256=sha(OLD/'data/interior_section_search.csv')))
    legal=np.array([r[0]['relative_error'] for r in sections if r[0]['original_geometry_pass']])
    data=dict(status='PASS',case=str(BASE),source_field_sha256=sha(BASE/'frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu'),
        Q_in_m3_s=qin,outlet_Q_m3_s=outs,outlet_fractions={k:q/qin for k,q in outs.items()},
        global_signed_flux_residual_m3_s=math.fsum(flux.values()),global_mass_error=abs(math.fsum(flux.values()))/qin,
        wall=wall,divergence=records,candidate_section_count=len(candidates),legal_root_section_count=len(legal),
        root_section_max_error=float(abs(legal).max()),root_section_RMS_error=float(np.sqrt(np.mean(legal**2))),
        six_slabs=slabs,maximum_slab_closure_relative=max(abs(x['closure_relative']) for x in slabs),
        root_region_definition='nearest root centerline versus daughter segments, before junction',
        section_gate='Only original geometry eligibility; no arbitrary internal flux tolerance',
        process_workers=6,elapsed_seconds=time.monotonic()-started)
    dump(REPORT/'data/p1_network_bc_local_conservation_baseline.json',data)
    print(json.dumps({k:data[k] for k in ['status','Q_in_m3_s','global_mass_error','root_section_max_error','root_section_RMS_error','elapsed_seconds']}),flush=True)

if __name__=='__main__':main()
