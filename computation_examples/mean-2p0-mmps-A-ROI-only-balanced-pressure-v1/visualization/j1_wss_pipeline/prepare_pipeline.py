"""Read-only derivation of the five J1 WSS visualization stages, in SI units."""
from pathlib import Path
import csv, hashlib, importlib.util, json
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
J1=np.array([92.,49.,111.])
RADIUS=11.

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def dump(p,obj):
    Path(p).write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=lambda a:a.tolist() if isinstance(a,np.ndarray) else a.item())+'\n')

def core_module():
    spec=importlib.util.spec_from_file_location('wss_core',HERE/'audit/wss_core.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core);return core

def spaced_indices(points,spacing):
    # Deterministic surface sampling. Selection does not depend on vector magnitude.
    tree=cKDTree(points);available=np.ones(len(points),bool);chosen=[]
    for i in np.lexsort(points.T):
        if available[i]:
            chosen.append(i);available[tree.query_ball_point(points[i],spacing)]=False
    return np.array(chosen,dtype=int)

def stats(v,area,region,quantity,unit):
    order=np.argsort(v);cdf=(np.cumsum(area[order])-.5*area[order])/area.sum()
    p=np.interp([.05,.5,.95],cdf,v[order])
    return dict(region=region,quantity=quantity,unit=unit,facets=len(v),area_um2=area.sum()*1e12,
        min=v.min(),mean=np.average(v,weights=area),P5=p[0],P50=p[1],P95=p[2],max=v.max())

def main():
    if (HERE/'COMPUTE_VALIDATION.json').exists():raise FileExistsError('Preserve completed derived data')
    sources=['input_data/frozen_flow/flow_arrays_si.npz','input_data/frozen_flow/steady_flow_mean_2p0_mmps.vtu',
        'input_data/field_diagnostics/COMPUTE_VALIDATION.json','input_data/field_diagnostics/data/wall_wss_si.vtp',
        'render_visualization.py','render_surface_fields.py']
    sources += [f'input_data/solver_mesh/mesh-surfaces/{r}.vtp' for r in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']]
    lock={s:sha(ROOT/s) for s in sources};dump(HERE/'audit/source_lock.json',lock)
    c=json.loads((ROOT/sources[2]).read_text());assert c['all_pass']
    assert sha(HERE/'audit/wss_core.py')==c['core_sha256']
    assert lock[sources[0]]==c['flow_arrays_sha256'] and lock[sources[1]]==c['source_field_sha256']
    core=core_module();a=np.load(ROOT/sources[0]);x=a['points_m'];tet=a['tetra'];u=a['velocity_m_s']
    tag=a['facet_tags'];tri=a['boundary_triangles'][tag==1];wall_ids=np.flatnonzero(tag==1)
    g=core.p1_gradients(x,tet,u);owner=core.boundary_owners(tet,tri)
    centers,area,n=core.wall_geometry(x,tet,tri,owner);gw=g[owner];mu=c['material']['mu_Pa_s']
    traction=mu*np.einsum('nij,nj->ni',gw+gw.swapaxes(1,2),n)
    normal_signed=np.einsum('ij,ij->i',traction,n)
    tau=core.tangential_traction(gw,n,mu);wss=np.linalg.norm(tau,axis=1)
    gn=np.einsum('nij,nj->ni',gw,n);gnmag=np.linalg.norm(gn,axis=1)
    # Independent projection using explicit 3x3 projectors.
    tau2=np.einsum('nij,nj->ni',np.eye(3)-n[:,:,None]*n[:,None,:],traction)
    original=pv.read(ROOT/sources[3])
    assert np.array_equal(original.cell_data['Parent_tetra_zero_based'],owner)
    assert np.array_equal(original.cell_data['Global_boundary_facet_zero_based'],wall_ids)
    assert np.allclose(n,original.cell_data['Outward_normal'],rtol=0,atol=1e-14)
    assert np.max(abs(wss-original.cell_data['WSS_raw_Pa']))<1e-12
    assert np.max(abs(tau-original.cell_data['Tangential_viscous_traction_Pa']))<1e-12
    assert np.max(abs(tau-tau2))<1e-12
    assert np.max(abs(np.linalg.norm(n,axis=1)-1))<1e-14
    assert np.max(np.linalg.norm(u[np.unique(tri)],axis=1))==0
    full,ids=core.surface(x,tri)
    values={'VelocityGradient_s_inv':gw.reshape(-1,9),'Outward_normal':n,
        'ViscousTraction_Pa':traction,'NormalViscousTraction_signed_Pa':normal_signed,
        'WSSVector_Pa':tau,'FluidOnWallShear_Pa':-tau,'WSS_raw_Pa':wss,
        'NormalVelocityDerivative_s_inv':gn,'NormalVelocityDerivativeMagnitude_s_inv':gnmag,
        'Area_m2':area,'Parent_tetra_zero_based':owner,'Global_boundary_facet_zero_based':wall_ids,
        'Boundary_tag':np.ones(len(tri),dtype=np.int32),
        'TractionAngleFromTangent_deg':np.degrees(np.arctan2(abs(normal_signed),wss))}
    for name,v in values.items():full.cell_data[name]=v
    # Vertex normals are area-weighted incident full-wall facet normals. They
    # illustrate nodal orientation only; production traction uses facet normals.
    node_n=np.zeros_like(x)
    for k in range(3):np.add.at(node_n,tri[:,k],area[:,None]*n)
    node_n[ids]/=np.linalg.norm(node_n[ids],axis=1)[:,None]
    full.point_data['NodeOutwardNormal']=node_n[ids]
    display,_=core.nodal_average(tri,wss,area,len(x))
    assert np.array_equal(original.point_data['GlobalNodeID_zero_based'],ids)
    assert np.max(abs(display[ids]-original.point_data['WSS_display_Pa']))<1e-12
    full.point_data['WSS_display_Pa']=original.point_data['WSS_display_Pa'].copy()
    ng,_=core.nodal_average(tri,gnmag,area,len(x))
    full.point_data['NormalVelocityDerivative_display_s_inv']=ng[ids]
    sel=np.linalg.norm(centers*1e6-J1,axis=1)<RADIUS
    # Whole existing triangles, with no clipping-generated faces or capping.
    local=full.extract_cells(sel).extract_surface(algorithm='dataset_surface')
    assert local.n_cells==sel.sum()
    local.save(HERE/'data/J1_wss_pipeline_si.vtp')
    lc=local.cell_centers().points*1e6;ln=local.cell_data['Outward_normal']
    fi=spaced_indices(lc,1.25)
    ni=spaced_indices(local.points*1e6,1.25)
    # Identical face samples and identical linear arrow scale for steps 3/4/5.
    np.savez_compressed(HERE/'data/glyph_samples.npz',face_indices=fi,node_indices=ni)
    local_tau=local.cell_data['WSSVector_Pa'];local_t=local.cell_data['ViscousTraction_Pa']
    local_max_t=float(np.linalg.norm(local_t,axis=1).max())
    arrow_scale=1.15/local_max_t # micrometres per Pa, no normalization floor.
    # P1 no-slip identity: with three zero wall-node values G is rank one;
    # n.t_mu = 2 mu div(u), and the chosen wall-on-fluid tau opposes interior u_t.
    near_u=u[tet[owner]].mean(axis=1);near_ut=near_u-(near_u*n).sum(1)[:,None]*n
    dots=np.einsum('ij,ij->i',tau,near_ut)
    normal_identity_error=np.max(abs(normal_signed-2*mu*np.trace(gw,axis1=1,axis2=2)))
    valid=(np.linalg.norm(near_ut,axis=1)>1e-14)&(wss>1e-12)
    cos=dots[valid]/(wss[valid]*np.linalg.norm(near_ut[valid],axis=1))
    assert normal_identity_error<1e-10 and np.max(cos)<-.999999999
    # Affine simple-shear sign and projection sanity check in the production core.
    xx=np.array([[0,0,0],[2,0,0],[0,3,0],[0,0,4]])*1e-6
    gg=np.array([[0,0,1000],[0,0,0],[0,0,0]],float)
    recovered=core.p1_gradients(xx,np.array([[0,1,2,3]]),xx@gg.T)
    analytic=core.tangential_traction(recovered,np.array([[0.,0.,-1.]]),mu)[0]
    expected=np.array([-mu*1000,0,0]);assert np.max(abs(analytic-expected))<1e-12
    # Annotation coordinates from previously established graph paths. No WSS
    # values or velocity from old cases are reused.
    rows=list(csv.DictReader((HERE/'audit/reference_paths_um.csv').open()))
    pathdata={}
    for name in ['INLET_to_J1','J1_to_O2','J1_to_J2']:
        rr=[r for r in rows if r['path']==name]
        ss=np.array([float(r['s_um']) for r in rr]);xyz=np.array([[float(r[k]) for k in ['x_um','y_um','z_um']] for r in rr])
        assert np.linalg.norm((xyz[-1] if name=='INLET_to_J1' else xyz[0])-J1)<1e-10
        s=ss[-1]-8.5 if name=='INLET_to_J1' else 8.5
        target=np.array([np.interp(s,ss,xyz[:,k]) for k in range(3)])
        assert np.linalg.norm(target-J1)<RADIUS
        dist,face=cKDTree(lc).query(target)
        pathdata[name]=dict(centerline_target_um=target,wall_arrow_target_um=lc[face],
            local_face_id=int(face),distance_to_wall_um=float(dist))
    ports=[('J1',J1),('Toward inlet',pathdata['INLET_to_J1']['wall_arrow_target_um']),
        ('Toward O2',pathdata['J1_to_O2']['wall_arrow_target_um']),
        ('Toward J2',pathdata['J1_to_J2']['wall_arrow_target_um'])]
    # Record all actual cap coordinates for an unambiguous anatomical identity.
    caps={}
    for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        cap=pv.read(ROOT/f'input_data/solver_mesh/mesh-surfaces/{role}.vtp')
        caps[role]=(np.asarray(cap.center)*1e6).tolist()
    capdiff={}
    for path,role,first in [('INLET_to_J1','INLET',True),('J1_to_O2','OUTLET_02',False)]:
        rr=[r for r in rows if r['path']==path];r=rr[0 if first else -1]
        old=np.array([float(r[k]) for k in ['x_um','y_um','z_um']])
        capdiff[role]=float(np.linalg.norm(old-caps[role]));assert capdiff[role]<.05
    dump(HERE/'data/region_and_annotations.json',dict(J1_um=J1,display_radius_um=RADIUS,
        region_rule='Original wall triangle centroid within 11 um of J1; open crop edges are display cuts, not outlets',
        paths=pathdata,annotations=ports,current_cap_centers_um=caps,graph_cap_discrepancy_um=capdiff))
    summaries=[]
    for rad in [5.,RADIUS]:
        mask=np.linalg.norm(centers*1e6-J1,axis=1)<rad
        for name,v,unit in [('WSS_raw',wss,'Pa'),('ViscousTractionMagnitude',np.linalg.norm(traction,axis=1),'Pa'),
            ('NormalViscousTraction_abs',abs(normal_signed),'Pa'),('NormalVelocityDerivativeMagnitude',gnmag,'s^-1'),
            ('TractionAngleFromTangent',values['TractionAngleFromTangent_deg'],'degree')]:
            summaries.append(stats(v[mask],area[mask],f'J1_r{rad:g}um',name,unit))
    with (HERE/'data/J1_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summaries[0]));w.writeheader();w.writerows(summaries)
    with (HERE/'data/J1_facet_values.csv').open('w',newline='') as f:
        names=['global_boundary_facet','parent_tetra','x_um','y_um','z_um','area_um2','nx','ny','nz',
            'Gxx_s_inv','Gxy_s_inv','Gxz_s_inv','Gyx_s_inv','Gyy_s_inv','Gyz_s_inv','Gzx_s_inv','Gzy_s_inv','Gzz_s_inv',
            't_mu_x_Pa','t_mu_y_Pa','t_mu_z_Pa','normal_traction_Pa','tau_x_Pa','tau_y_Pa','tau_z_Pa','WSS_raw_Pa']
        w=csv.writer(f);w.writerow(names)
        for i in np.flatnonzero(sel):w.writerow([wall_ids[i],owner[i],*(centers[i]*1e6),area[i]*1e12,*n[i],*gw[i].ravel(),*traction[i],normal_signed[i],*tau[i],wss[i]])
    result=dict(all_pass=True,case=c['case'],source_field_sha256=c['source_field_sha256'],source_lock=lock,
        mu_Pa_s=mu,model='Current Newtonian rigid-wall no-slip FEM field',J1_um=J1,display_radius_um=RADIUS,
        wall_facets=len(tri),local_facets=int(sel.sum()),local_nodes=local.n_points,face_glyph_count=len(fi),node_glyph_count=len(ni),
        glyph_spacing_um=1.25,normal_glyph_length_um=.9,traction_glyph_scale_um_per_Pa=arrow_scale,
        glyph_offset_um=.075,gradient_core_sha256=sha(HERE/'audit/wss_core.py'),
        WSS_reproduction_max_error_Pa=float(np.max(abs(wss-original.cell_data['WSS_raw_Pa']))),
        projection_crosscheck_max_error_Pa=float(np.max(abs(tau-tau2))),
        tangency_max_error_Pa=float(abs(np.einsum('ij,ij->i',tau,n)).max()),
        normal_identity_max_error_Pa=float(normal_identity_error),wall_speed_max_m_s=0.,
        near_wall_velocity_cosine_range=[float(cos.min()),float(cos.max())],
        affine_simple_shear=dict(expected_Pa=expected,computed_Pa=analytic,absolute_error_Pa=float(np.max(abs(analytic-expected)))),
        source_files_unchanged=True,CFD_calls=0,scalar_clipping=False,geometry_smoothing=False,
        scalar_display='Existing full-wall area-weighted nodal WSS magnitude; raw face counterpart retained',
        gradient_display='Full-wall area-weighted nodal mean of norm(G n), not a tensor component or WSS',
        traction_vector_convention='outward-fluid normal; wall-on-fluid tangential traction; -tau is fluid-on-wall',
        statistics=summaries,files={str(p.relative_to(HERE)):sha(p) for p in sorted((HERE/'data').glob('*'))})
    for name,h in lock.items():assert sha(ROOT/name)==h
    dump(HERE/'COMPUTE_VALIDATION.json',result)
    print(json.dumps({k:result[k] for k in ['local_facets','local_nodes','face_glyph_count','node_glyph_count','tangency_max_error_Pa','normal_identity_max_error_Pa','statistics']},default=lambda x:x.item()),flush=True)

if __name__=='__main__':main()
