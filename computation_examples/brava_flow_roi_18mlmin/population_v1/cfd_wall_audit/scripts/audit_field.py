"""Read-only canonical P1 whole-domain audit. No cached gradients or CFD execution."""
from pathlib import Path
import argparse,hashlib,json,itertools,sys,time
import numpy as np
import pandas as pd
import pyvista as pv
import vtk
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def js(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item())+'\n')
def ratio(a,b):return np.divide(a,b,out=np.zeros_like(a,dtype=float),where=b!=0)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--environment',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--commit',required=True);args=ap.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True);(out/'validation').mkdir(exist_ok=True)
    spec=json.loads(args.environment.read_text());paths=dict(mesh=spec['mesh'],flow=spec['flow'],**spec['boundaries'])
    hashes={k:sha(v) for k,v in paths.items()};assert hashes==spec['sha256']
    m,f,w=(pv.read(paths[k]) for k in ['mesh','flow','WALL']);N=m.n_cells
    assert np.all(m.celltypes==10) and np.all(f.celltypes==10)
    c=np.asarray(m.cells).reshape(-1,5)[:,1:];fc=np.asarray(f.cells).reshape(-1,5)[:,1:]
    xyz=np.asarray(m.points,dtype=float);u=np.asarray(f['Velocity']);p=np.asarray(f['Pressure'])
    assert np.array_equal(xyz,f.points) and np.array_equal(np.sort(c,axis=1),np.sort(fc,axis=1))
    assert np.array_equal(m['GlobalNodeID'],np.arange(1,len(xyz)+1))
    assert np.array_equal(m['GlobalElementID'],np.arange(1,N+1))
    for name in spec['boundaries']:
        b=pv.read(paths[name]);assert np.array_equal(b.points,xyz[np.asarray(b['GlobalNodeID'])-1])
    assert np.isfinite(u).all() and np.isfinite(p).all()
    wg=np.asarray(w['GlobalNodeID'],dtype=np.int64)-1;wc=wg[np.asarray(w.faces).reshape(-1,4)[:,1:]]
    identity=dict(input_paths=paths,sha256=hashes,source_code_commit=args.commit,audit_script_sha256=sha(__file__),node_count=len(xyz),tetra_count=N,wall_triangle_count=len(wc),units=dict(x='m',u='m/s',p='Pa',mu='Pa.s',gradient='1/s',volume='m3',flux='m3/s'),mu_Pa_s=spec['mu_Pa_s'],numpy=np.__version__,pyvista=pv.__version__,vtk=vtk.vtkVersion.GetVTKVersion(),mapping=dict(canonical_GlobalNodeID_contiguous=True,surface_GlobalNodeID_coordinate_equality=True,flow_has_GlobalNodeID='GlobalNodeID' in f.point_data,flow_coordinate_order_exact=True,canonical_cell_vertex_sets_exact=True,flow_local_tetra_vertex_permutations=int(np.any(c!=fc,axis=1).sum()),limitation='Frozen flow lacks GlobalNodeID; identity independently established by exact ordered coordinates and cellwise vertex sets, not by an unavailable field ID.'))
    js(out/'AUDIT_IDENTITY.json',identity)
    js(out/'AUDIT_CONFIG.json',dict(W0='owns official WALL face',W1='face neighbors of W0 excluding W0',W2='face neighbors of W1 excluding W0,W1',distance='exact Euclidean point-to-finite-WALL-triangle, vtkStaticCellLocator, centroid',h_cell='maximum edge',bands=[0,.5,1,2,4,'inf'],eta_div='abs(trace(G))/norm(G); zero only when denominator exactly zero',tiers=[0,.01,.1,.5,'inf'],severity='SCREENING_ONLY_NOT_PHYSICAL_PASS_FAIL',normal='inward; W0 official own face, other cells nearest finite wall triangle; cells with several WALL faces use largest-area face for cell summaries',strain='(G+G.T)/2',clusters='face adjacency within W0/W1/W2; eta>=0.5 or abs(un)/speed>=0.5; separate families',Q_in_m3_s=3e-7,fast_interior_screen_m_s=.1))
    print('Identity verified; recomputing all P1 tensors',flush=True)
    X=xyz[c];U=u[c];P=p[c];cent=X.mean(1);uc=U.mean(1);pc=P.mean(1)
    mat=X[:,1:]-X[:,:1];G=np.linalg.solve(mat,U[:,1:]-U[:,:1]).transpose(0,2,1)
    gp=np.linalg.solve(mat,(P[:,1:]-P[:,:1])[...,None])[...,0]
    det=np.linalg.det(mat);vol=np.abs(det)/6;assert np.all(vol>0)
    div=np.trace(G,axis1=1,axis2=2);gn=np.linalg.norm(G,axis=(1,2));eta=ratio(abs(div),gn)
    strain=(G+G.transpose(0,2,1))/2;curl=np.stack([G[:,2,1]-G[:,1,2],G[:,0,2]-G[:,2,0],G[:,1,0]-G[:,0,1]],axis=1)
    edges=np.stack([np.linalg.norm(X[:,a]-X[:,b],axis=1) for a,b in itertools.combinations(range(4),2)],axis=1)
    faces=np.array([[1,2,3],[0,2,3],[0,1,3],[0,1,2]])
    av=[];q=[]
    for k,idx in enumerate(faces):
        v=np.cross(X[:,idx[1]]-X[:,idx[0]],X[:,idx[2]]-X[:,idx[0]])/2
        v*=np.where(np.einsum('ij,ij->i',v,X[:,k]-X[:,idx].mean(1))>0,-1,1)[:,None]
        av.append(v);q.append(np.einsum('ij,ij->i',v,U[:,idx].mean(1)))
    av=np.stack(av,axis=1);area=np.linalg.norm(av,axis=2);fn=av/area[:,:,None];q=np.stack(q,axis=1)
    qnet=q.sum(1);err=qnet-div*vol
    hmin=3*vol/area.max(1);hmax=edges.max(1);deq=(6*vol/np.pi)**(1/3)
    shape=12*(3*vol)**(2/3)/np.sum(edges**2,axis=1)
    corner=[]
    for k in range(4):
        dx=X[:,[j for j in range(4) if j!=k]]-X[:,k,None]
        corner.append(np.sqrt(2)*abs(np.linalg.det(dx))/np.prod(np.linalg.norm(dx,axis=2),axis=1))
    sj=np.min(corner,axis=0)
    angles=np.stack([np.arccos(np.clip(-np.einsum('ij,ij->i',fn[:,a],fn[:,b]),-1,1)) for a,b in itertools.combinations(range(4),2)],axis=1)*180/np.pi
    print('Building exact canonical face adjacency',flush=True)
    cf=np.sort(c[:,faces].reshape(-1,3),axis=1)
    keys=np.ascontiguousarray(cf).view(np.dtype((np.void,cf.dtype.itemsize*3))).ravel()
    uniq,first,inv,cnt=np.unique(keys,return_index=True,return_inverse=True,return_counts=True);assert cnt.max()==2
    order=np.argsort(inv,kind='stable');start=np.r_[0,np.cumsum(cnt)[:-1]]
    both=np.flatnonzero(cnt==2);aa=order[start[both]]//4;bb=order[start[both]+1]//4
    adj=np.full((N,4),-1,dtype=np.int64);aaidx=order[start[both]];bbidx=order[start[both]+1]
    adj.ravel()[aaidx]=bbidx//4;adj.ravel()[bbidx]=aaidx//4
    wk=np.ascontiguousarray(np.sort(wc,axis=1)).view(keys.dtype).ravel();wi=np.searchsorted(uniq,wk)
    assert np.array_equal(uniq[wi],wk) and np.all(cnt[wi]==1)
    owner=first[wi]//4;ownface=first[wi]%4;opp=c[owner,ownface]
    assert np.array_equal(np.asarray(w.cell_data['GlobalElementID'])-1,owner)
    layer=np.full(N,3,dtype=np.int8);layer[np.unique(owner)]=0
    for k in [1,2]:
        nb=adj[layer==k-1].ravel();nb=nb[nb>=0];layer[nb[layer[nb]==3]]=k
    normal=-fn[owner,ownface];wa=area[owner,ownface];wfcenter=xyz[wc].mean(1)
    # Cell scalar for W0 takes largest-area owned WALL face; retain all faces separately.
    choose=np.full(N,-1,dtype=int)
    for i in np.argsort(wa):choose[owner[i]]=i
    print('Exact finite-WALL centroid distances for all tetra',flush=True)
    locator=vtk.vtkStaticCellLocator();locator.SetDataSet(w);locator.BuildLocator()
    nearest=np.empty(N,dtype=int);distance=np.empty(N);closest=np.empty((N,3));point=[0.,0.,0.];cid=vtk.reference(0);sub=vtk.reference(0);d2=vtk.reference(0.)
    for i in range(N):
        locator.FindClosestPoint(cent[i],point,cid,sub,d2);nearest[i]=int(cid);distance[i]=float(d2)**.5;closest[i]=point
    normals=normal[nearest].copy();isw=layer==0;normals[isw]=normal[choose[isw]]
    un=np.einsum('ij,ij->i',uc,normals);speed=np.linalg.norm(uc,axis=1);nf=ratio(abs(un),speed)
    bands=np.searchsorted([.5,1,2,4],distance/hmax,side='right');labels=['0-0.5','0.5-1','1-2','2-4','>=4']
    metrics=dict(tetra_id=np.arange(N),GlobalElementID=np.arange(N)+1,layer=layer,distance_band=bands,volume_m3=vol,centroid_x_m=cent[:,0],centroid_y_m=cent[:,1],centroid_z_m=cent[:,2],div_s_inv=div,abs_div_s_inv=abs(div),G_norm_s_inv=gn,eta_div=eta,speed_m_s=speed,pressure_Pa=pc,grad_p_norm_Pa_m=np.linalg.norm(gp,axis=1),strain_norm_s_inv=np.linalg.norm(strain,axis=(1,2)),curl_norm_s_inv=np.linalg.norm(curl,axis=1),u_normal_m_s=un,abs_u_normal_m_s=abs(un),normal_fraction=nf,toward_wall=(un<0),d_wall_m=distance,d_over_hmax=distance/hmax,h_max_m=hmax,h_min_altitude_m=hmin,equivalent_diameter_m=deq,aspect_hmax_hmin=hmax/hmin,mean_ratio_quality=shape,scaled_jacobian_abs_regular_normalized=sj,signed_det_m3=det,min_dihedral_deg=angles.min(1),Q_net_m3_s=qnet,Q_div_m3_s=div*vol,flux_identity_error_m3_s=err,all_four_u_zero=np.all(U==0,axis=(1,2)))
    frame=pd.DataFrame(metrics);frame.to_csv(out/'all_cells.csv',index=False)
    near=layer<=2;wide=frame.loc[near].copy()
    for k in range(4):
        wide[f'node{k}_GlobalNodeID']=c[near,k]+1
        for a,axis in enumerate('xyz'):
            wide[f'node{k}_{axis}_m']=X[near,k,a];wide[f'node{k}_u{axis}_m_s']=U[near,k,a]
        wide[f'node{k}_p_Pa']=P[near,k]
    for k,(a,b) in enumerate(itertools.combinations(range(4),2)):wide[f'edge_{a}{b}_m']=edges[near,k]
    wide.to_csv(out/'wall_adjacent_cells.csv',index=False)
    np.savez_compressed(out/'recomputed_arrays.npz',xyz=xyz,velocity=u,pressure=p,connectivity=c,G=G,strain=strain,curl=curl,grad_p=gp,adjacency=adj,face_flux=q,layer=layer,eta=eta,div=div,volume=vol,centroid=cent,owner=owner,wall_connectivity=wc,wall_normal_in=normal,wall_area=wa,nearest_wall_triangle=nearest,nearest_wall_point=closest,wall_distance=distance)
    wn=np.unique(wc);ws=np.linalg.norm(u[wn],axis=1)
    nod=pd.DataFrame(dict(GlobalNodeID=wn+1,x_m=xyz[wn,0],y_m=xyz[wn,1],z_m=xyz[wn,2],ux_m_s=u[wn,0],uy_m_s=u[wn,1],uz_m_s=u[wn,2],speed_m_s=ws));nod.to_csv(out/'wall_nodes.csv',index=False)
    nonzero=[]
    for nid in wn[ws!=0]:nonzero.append(dict(GlobalNodeID=int(nid+1),wall_triangles=np.flatnonzero(np.any(wc==nid,axis=1)),tetra=np.flatnonzero(np.any(c==nid,axis=1)),position_m=xyz[nid],velocity_m_s=u[nid]))
    js(out/'validation/nonzero_wall_nodes.json',nonzero)
    wsmax=np.linalg.norm(u[wc],axis=2).max(1);uo=u[opp];uo_n=np.einsum('ij,ij->i',uo,normal);uo_s=np.linalg.norm(uo,axis=1)
    wf=pd.DataFrame(dict(wall_triangle_id=np.arange(len(wc)),tetra_id=owner,opposite_GlobalNodeID=opp+1,area_m2=wa,wall_face_node_speed_max_m_s=wsmax,interior_vertex_speed_m_s=uo_s,interior_vertex_height_m=3*vol[owner]/wa,interior_vertex_normal_m_s=uo_n,interior_vertex_tangent_m_s=np.linalg.norm(uo-uo_n[:,None]*normal,axis=1),eta_div=eta[owner],div_s_inv=div[owner],hmax_m=hmax[owner]))
    for j in range(3):wf[f'wall_node{j}_speed_m_s']=np.linalg.norm(u[wc[:,j]],axis=1)
    # Normal changes at real shared WALL edges, face-edge adjacency only.
    wedges=np.sort(wc[:,[[0,1],[1,2],[0,2]]].reshape(-1,2),axis=1)
    ek=np.ascontiguousarray(wedges).view(np.dtype((np.void,16))).ravel();_,ef,ei,ec=np.unique(ek,return_index=True,return_inverse=True,return_counts=True)
    eo=np.argsort(ei,kind='stable');es=np.r_[0,np.cumsum(ec)[:-1]];eb=np.flatnonzero(ec==2);ea=eo[es[eb]]//3;eb2=eo[es[eb]+1]//3
    turn=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',normal[ea],normal[eb2]),-1,1)));maxturn=np.zeros(len(wc));np.maximum.at(maxturn,ea,turn);np.maximum.at(maxturn,eb2,turn);wf['max_adjacent_wall_normal_angle_deg']=maxturn
    traction=spec['mu_Pa_s']*np.einsum('nij,nj->ni',G[owner]+G[owner].transpose(0,2,1),-normal)
    tang=traction-np.einsum('ij,ij->i',traction,normal)[:,None]*normal;wss=np.linalg.norm(tang,axis=1);wf['recomputed_WSS_Pa']=wss
    oldpath=Path(spec['flow']).parent.parent/'postprocess/wall_wss_si.vtp';old=pv.read(oldpath)
    og=np.asarray(old['GlobalNodeID_zero_based'],dtype=np.int64);oc=og[np.asarray(old.faces).reshape(-1,4)[:,1:]]
    # Match triangles by canonical global vertex sets, not cell row assumption.
    ok=np.ascontiguousarray(np.sort(oc,axis=1)).view(keys.dtype).ravel();oi=np.argsort(ok);mapidx=oi[np.searchsorted(ok[oi],wk)];assert np.array_equal(ok[mapidx],wk)
    assert np.array_equal(np.asarray(old['Parent_tetra_zero_based'])[mapidx],owner)
    diff=wss-np.asarray(old['WSS_raw_Pa'])[mapidx];wf['existing_WSS_Pa']=np.asarray(old['WSS_raw_Pa'])[mapidx];wf['WSS_difference_Pa']=diff;wf.to_csv(out/'wall_face_structure.csv',index=False)
    sample=[]
    for name,frac in [('tetra_centroid',.25),('wall_to_opposite_10pct',.1),('wall_to_opposite_25pct',.25),('wall_to_opposite_50pct',.5),('wall_to_opposite_75pct',.75)]:
        sx=(1-frac)*wfcenter+frac*xyz[opp];su=(1-frac)*u[wc].mean(1)+frac*uo;sn=np.einsum('ij,ij->i',su,normal);ss=np.linalg.norm(su,axis=1)
        sample.append(pd.DataFrame(dict(wall_triangle_id=np.arange(len(wc)),tetra_id=owner,sample=name,fraction=frac,x_m=sx[:,0],y_m=sx[:,1],z_m=sx[:,2],u_normal_m_s=sn,u_tangent_m_s=np.linalg.norm(su-sn[:,None]*normal,axis=1),speed_m_s=ss,normal_fraction=ratio(abs(sn),ss),toward_wall=sn<0)))
    pd.concat(sample,ignore_index=True).to_csv(out/'wall_sample_metrics.csv',index=False)
    print('Statistics and connected components',flush=True)
    groups={f'W{k}':layer==k for k in range(3)};groups.update({f'd_over_hmax_{label}':bands==k for k,label in enumerate(labels)})
    stats=[];fracs=[];corr=[]
    for name,mask in groups.items():
        for met in ['abs_div_s_inv','eta_div','G_norm_s_inv','speed_m_s','abs_u_normal_m_s','normal_fraction','h_max_m','h_min_altitude_m','mean_ratio_quality','min_dihedral_deg','pressure_Pa','grad_p_norm_Pa_m']:
            values=frame.loc[mask,met].values
            if len(values):stats.append(dict(region=name,metric=met,**dict(zip(['min','median','P90','P95','P99','max'],np.quantile(values,[0,.5,.9,.95,.99,1])))))
        for met in ['eta_div','normal_fraction']:
            for lo,hi in [(0,.01),(.01,.1),(.1,.5),(.5,np.inf),(.1,np.inf)]:
                sel=mask&(frame[met].values>=lo)&(frame[met].values<hi)
                fracs.append(dict(region=name,metric=met,lower=lo,upper=hi,total_count=int(mask.sum()),flag_count=int(sel.sum()),count_fraction=float(sel.sum()/mask.sum()) if mask.any() else 0,total_volume_m3=vol[mask].sum(),flag_volume_m3=vol[sel].sum(),volume_fraction=vol[sel].sum()/vol[mask].sum() if mask.any() else 0))
        for met in ['h_max_m','mean_ratio_quality','min_dihedral_deg','normal_fraction']:
            if mask.sum()>1:corr.append(dict(region=name,x='eta_div',y=met,spearman=float(spearmanr(eta[mask],frame.loc[mask,met]).statistic)))
    pd.DataFrame(stats).to_csv(out/'layer_distributions.csv',index=False);pd.DataFrame(fracs).to_csv(out/'layer_tiers.csv',index=False);pd.DataFrame(corr).to_csv(out/'correlations.csv',index=False)
    clusterrows=[];clmembers=[];wallmembers=[];cllabels={}
    def components(flag,family):
        ids=np.flatnonzero(flag);index=np.full(N,-1,dtype=int);index[ids]=np.arange(len(ids));mask=flag[aa]&flag[bb]
        graph=coo_matrix((np.ones(mask.sum()),(index[aa[mask]],index[bb[mask]])),shape=(len(ids),len(ids))).tocsr()
        n,lab=connected_components(graph,directed=False);global_lab=np.full(N,-1,dtype=int);global_lab[ids]=lab
        size=np.bincount(lab);order=np.argsort(-size);rank=np.empty(n,dtype=int);rank[order]=np.arange(n);global_lab[ids]=rank[lab];cllabels[family]=global_lab
        cellorder=np.argsort(global_lab[ids],kind='stable');sortedids=ids[cellorder];starts=np.r_[0,np.cumsum(size[order])]
        walllab=global_lab[owner]
        for k in range(n):
            ci=sortedids[starts[k]:starts[k+1]];wi2=np.flatnonzero(walllab==k);bounds=np.stack([X[ci].min(axis=(0,1)),X[ci].max(axis=(0,1))])
            clusterrows.append(dict(family=family,cluster_id=k,cell_count=len(ci),volume_m3=vol[ci].sum(),wall_area_m2=wa[wi2].sum(),bbox_min_x_m=bounds[0,0],bbox_min_y_m=bounds[0,1],bbox_min_z_m=bounds[0,2],bbox_max_x_m=bounds[1,0],bbox_max_y_m=bounds[1,1],bbox_max_z_m=bounds[1,2],eta_max=eta[ci].max(),eta_median=np.median(eta[ci]),normal_fraction_max=nf[ci].max(),normal_fraction_median=np.median(nf[ci])))
        clmembers.append(pd.DataFrame(dict(family=family,tetra_id=ids,cluster_id=global_lab[ids])))
        wx=walllab>=0;wallmembers.append(pd.DataFrame(dict(family=family,wall_triangle_id=np.flatnonzero(wx),tetra_id=owner[wx],cluster_id=walllab[wx])))
        return global_lab,n
    for family,flag in [('eta_ge_0p5',near&(eta>=.5)),('normal_ge_0p5',near&(nf>=.5)),('toward_wall_normal_ge_0p5',near&(nf>=.5)&(un<0)),('all_zero_velocity',np.all(U==0,axis=(1,2)))]:components(flag,family)
    pd.DataFrame(clusterrows).to_csv(out/'anomaly_clusters.csv',index=False);pd.concat(clmembers).to_csv(out/'cluster_cells.csv',index=False);pd.concat(wallmembers).to_csv(out/'cluster_wall_triangles.csv',index=False)
    np.savez_compressed(out/'cluster_labels.npz',**cllabels)
    severe=near&((eta>=.5)|((nf>=.5)&(un<0)));frame[severe].sort_values('eta_div',ascending=False).to_csv(out/'worst_cells.csv',index=False)
    sc=np.flatnonzero(severe);src=np.repeat(sc,4);nb=adj[sc].ravel();valid=nb>=0
    ndf=frame.iloc[nb[valid]].copy();ndf.insert(0,'flagged_tetra_id',src[valid]);ndf.to_csv(out/'flagged_cell_neighbors.csv',index=False)
    pressure_edges=pd.DataFrame(dict(tetra_a=aa,tetra_b=bb,centroid_pressure_difference_Pa=abs(pc[aa]-pc[bb]),centroid_distance_m=np.linalg.norm(cent[aa]-cent[bb],axis=1),gradient_jump_Pa_m=np.linalg.norm(gp[aa]-gp[bb],axis=1)))
    pressure_edges['scaled_centroid_gradient_Pa_m']=pressure_edges.centroid_pressure_difference_Pa/pressure_edges.centroid_distance_m
    pressure_edges.nlargest(2000,'gradient_jump_Pa_m').to_csv(out/'pressure_neighbor_extremes.csv',index=False)
    for tid in [388769,138584]:
        nb=adj[tid];nb=nb[nb>=0];wallids=np.flatnonzero(owner==tid);zeroid=int(cllabels['all_zero_velocity'][tid])
        js(out/f'known_regression_{tid}.json',dict(tetra_id=tid,GlobalNodeIDs=c[tid]+1,xyz_m=X[tid],velocity_m_s=U[tid],pressure_Pa=P[tid],G_s_inv=G[tid],grad_p_Pa_m=gp[tid],metrics=frame.iloc[tid].to_dict(),neighbors=frame.iloc[nb].to_dict('records'),wall_triangles=wallids,wall_structure=wf.iloc[wallids].to_dict('records'),zero_component_id=zeroid,zero_component_cell_count=int(np.sum(cllabels['all_zero_velocity']==zeroid)) if zeroid>=0 else 0,all_four_nodes_on_WALL=bool(np.isin(c[tid],wn).all())))
    summary=dict(wall_nodes=dict(count=len(wn),max_m_s=ws.max(),median_m_s=np.median(ws),P95_m_s=np.quantile(ws,.95),P99_m_s=np.quantile(ws,.99),nonzero_count=int((ws!=0).sum())),layer_counts={f'W{k}':int((layer==k).sum()) for k in range(3)},flux_identity=dict(max_absolute_error_m3_s=max(abs(err)),max_error_scaled_by_sum_abs_face_flux=float(np.max(ratio(abs(err),np.sum(abs(q),axis=1))))),global_sum_Qnet_m3_s=qnet.sum(),WSS_comparison=dict(existing_sha256=sha(oldpath),max_absolute_difference_Pa=max(abs(diff)),rms_difference_Pa=float(np.sqrt(np.mean(diff**2)))),zero_cells=int(np.all(U==0,axis=(1,2)).sum()),wall_face_zero3_count=int((wsmax==0).sum()),wall_face_zero3_opposite_gt_0p1_m_s_count=int(((wsmax==0)&(uo_s>.1)).sum()),W0_toward_normal_ge_0p5_count=int((isw&(un<0)&(nf>=.5)).sum()),W0_normal_ge_0p5_count=int((isw&(nf>=.5)).sum()),wall_area_m2=wa.sum(),W0_eta_ge_0p5_wall_area_fraction=float(wa[eta[owner]>=.5].sum()/wa.sum()),pressure_neighbor_gradient_jump_quantiles_Pa_m=np.quantile(pressure_edges.gradient_jump_Pa_m,[.5,.95,.99,1]),regression_388769_pass=bool(abs(div[388769]+827.4291634662534)<1e-8 and abs(eta[388769]-.9948873431226156)<1e-10),minimum_volume_m3=vol.min(),det_positive_count=int((det>0).sum()),det_negative_count=int((det<0).sum()))
    js(out/'SUMMARY.json',summary)
    # Output derived-only VTK, no mutation of any authoritative input.
    for key in ['eta_div','div_s_inv','normal_fraction','layer','d_over_hmax','mean_ratio_quality']:m.cell_data[key]=metrics[key]
    m.extract_cells(near).save(out/'wall_layers_audit.vtu')
    w.cell_data['eta_div']=eta[owner];w.cell_data['normal_fraction']=ratio(abs(uo_n),uo_s);w.cell_data['WSS_recomputed_Pa']=wss;w.cell_data['canonical_owner_tetra']=owner;w.save(out/'wall_audit.vtp')
    assert hashes=={k:sha(v) for k,v in paths.items()};js(out/'validation/input_immutability.json',dict(before=hashes,after=hashes,unchanged=True))
    print((out/'SUMMARY.json').read_text(),flush=True)
if __name__=='__main__':main()
