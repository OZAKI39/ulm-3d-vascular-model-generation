"""Read-only project audit. All outputs and recovered code stay inside wss_audit.
Run with FEM_SimVascular/.venv/bin/python -B scripts/audit.py.
Production WSS functions are imported byte-for-byte from the recorded Git object.
"""
from pathlib import Path
import sys, json, csv, hashlib, shutil, subprocess, importlib.util, platform
import xml.etree.ElementTree as ET
import numpy as np
import scipy, scipy.spatial, scipy.sparse, scipy.sparse.csgraph
import pyvista as pv

A = Path(__file__).resolve().parents[1]
B = Path('/home/lzy/projects'); C = B/'formal_3D_flow_solver/FEM_SimVascular'
F = B/'ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular'
H = F/'flow_cases/mean-2p0-mmps-A-H0-pressure-v1'; O = F/'flow_cases/mean-2p0-mmps'
V = C/'rotate_visualization'; N = B/'ulm_3D_vascular'
ND = N/'reports/a_network_1d0d_boundary_v2_idealized/data'
sys.dont_write_bytecode = True
sys.path.insert(0, str(C/'solver_support/src/flow_solver_support'))
sys.path.insert(0, str(C/'mesh_generate/src'))
spec = importlib.util.spec_from_file_location('production_wss', A/'evidence/source/compute_field_diagnostics.py')
prod = importlib.util.module_from_spec(spec); spec.loader.exec_module(prod)
from vascular_validation.postprocess import SolutionMeasurements
from vascular_validation.mesh_diagnostics import tetra_sicn

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean(x):
    if isinstance(x, np.ndarray): return clean(x.tolist())
    if isinstance(x, np.generic): return clean(x.item())
    if isinstance(x, dict): return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x, (tuple,list)): return [clean(v) for v in x]
    return x
def dump(name, x): (A/name).write_text(json.dumps(clean(x), ensure_ascii=False, indent=2, allow_nan=False)+'\n')
def csvout(name, rows):
    with (A/name).open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(clean(rows))
def qtile(v,w,qs=(.01,.05,.5,.95,.99)):
    i=np.argsort(v); return np.interp(qs, (np.cumsum(w[i])-.5*w[i])/w.sum(),v[i])
def stats(v,w):
    return dict(min_Pa=v.min(),max_Pa=v.max(),mean_Pa=np.average(v,weights=w),
        **{f'p{int(q*100):02d}_Pa':x for q,x in zip((.01,.05,.5,.95,.99),qtile(v,w))},
        area_below_1Pa_pct=100*w[v<1].sum()/w.sum(), area_below_2Pa_pct=100*w[v<2].sum()/w.sum(),
        area_above_30Pa_pct=100*w[v>30].sum()/w.sum())

def main():
    protected={}
    def record(p, dest=None):
        p=Path(p); protected[str(p)]=sha(p)
        if dest:
            d=A/dest; d.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,d)
    code={
        'prepare_surface_data.py':V/'prepare_surface_data.py','render_surface_fields.py':V/'render_surface_fields.py',
        'render_visualization.py':V/'render_visualization.py','flow_parser.py':C/'solver_support/src/flow_solver_support/flow_parser.py',
        'postprocess.py':C/'mesh_generate/src/vascular_validation/postprocess.py',
        'mesh_diagnostics.py':C/'mesh_generate/src/vascular_validation/mesh_diagnostics.py',
        'generate_tetra_mesh.py':C/'mesh_generate/scripts/generate_tetra_mesh.py',
        'prepare_surface.py':C/'mesh_generate/scripts/prepare_surface.py',
        'import_surface_model.py':C/'mesh_generate/scripts/import_surface_model.py',
        'network_hydraulic_resistance.py':N/'network_1d0d/hydraulic_resistance.py',
        'network_solver.py':N/'network_1d0d/network_solver.py',
        'idealized_h0.py':N/'network_1d0d/idealized_h0.py',
        'extension_transfer.py':N/'network_1d0d/extension_transfer.py',
        'boundary_conditions.py':N/'network_1d0d/boundary_conditions.py',
        'fem_h0_case.py':N/'network_1d0d/fem_h0_case.py',
        'solve_a_h0_fem_remote.py':N/'scripts/solve_a_h0_fem_remote.py',
        'audit_mesh.py':C/'mesh_generate/scripts/audit_mesh.py',
        'geometry.py':C/'mesh_generate/src/vascular_validation/geometry.py',
        'swc_graph.py':N/'network_1d0d/swc_graph.py',
        'case_audit.py':N/'network_1d0d/case_audit.py'}
    for name,p in code.items(): record(p,'evidence/source/'+name)
    for name in ['fluid.cpp','post.cpp','set_bc.cpp','nn.cpp']:
        p=C/'external/flow_solver_source/Code/Source/solver'/name
        if p.exists():record(p,'evidence/source/'+name)
    for name in ['set_bc.cpp','fluid.cpp','post.cpp']:
        record(F/'vendor/svMultiPhysics_stage_q/Code/Source/solver'/name,'evidence/source/solver_frozen/'+name)
    record(V/'BUILD_CONTEXT.json','evidence/BUILD_CONTEXT.json')
    record(V/'results/surface_fields/figures/wss_overview_4k.png')
    record(C/'outputs/mesh_and_flow/mesh_generation/primary/generation.json','evidence/mesh/generation.json')
    for case,label in [(H,'H0'),(O,'zero_outlets')]:
        for path in ['frozen_flow/flow_arrays_si.npz','frozen_flow/manifest.json','run/solver.xml','policy.json','run/PETSC_OPTIONS.txt']:
            if (case/path).exists():record(case/path,f'inputs/{label}/'+Path(path).name)
    for p in [V/'input_data/field_diagnostics/data/wall_wss_si.vtp',V/'input_data/field_diagnostics/COMPUTE_VALIDATION.json',V/'input_data/BUILD_CONTEXT.json',V/'results/surface_fields/render_manifest.json']:
        if p.exists():record(p,'inputs/display/'+p.name)
    record(H/'SV_MESH/mesh_arrays.npz','inputs/mesh_arrays.npz')
    for name in ['roi_fixed_pressure_bc_H0.json','current_FEM_extension_sections.json','3d_case_config_diff.json','3d_physics_validation.json','analysis_A_H0_graph_si.npz']:
        record(ND/name,'inputs/network/'+name)
    portfile=N/'reports/a_network_1d0d_boundary_v1/data/roi_ports_in_a.json'
    record(portfile,'inputs/network/roi_ports_in_a.json')
    for name in ['mesh_policy.json','reference_physics.yaml','face_map.json']:
        record(C/'configs'/name,'inputs/'+name)
    record(H/'run/solver.log','evidence/solver.log')
    for p in sorted((C/'reports/mesh_and_flow').glob('*.json')): record(p,'evidence/mesh/'+p.name)
    build=subprocess.check_output(['git','show','fab4cf0eb3c8f8ad1ce40dfe995181c147cb9a9d:formal_3D_flow_solver/FEM_SimVascular/reports/sv1_3q/remote/svmp_reuse_build.json'],cwd=B/'ulm_flow_mean_2p0_mmps')
    (A/'evidence/solver_build_from_git.json').write_bytes(build)
    versions=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,pyvista=pv.__version__,vtk=pv.vtk_version_info,
        platform=platform.platform(),solver_git=subprocess.check_output(['git','describe','--always','--tags','--dirty'],cwd=C/'external/flow_solver_source',text=True).strip(),
        recovered_git_commit='fab4cf0eb3c8f8ad1ce40dfe995181c147cb9a9d',recovered_wss_sha=sha(A/'evidence/source/compute_field_diagnostics.py'))
    assert versions['recovered_wss_sha']=='efa12e218bcaf99ba07db7194d3b8da79251da709c0a54a97ad6f6179be6e438'
    dump('evidence/versions.json',versions)
    d=np.load(A/'inputs/H0/flow_arrays_si.npz'); old=np.load(A/'inputs/zero_outlets/flow_arrays_si.npz'); m=np.load(A/'inputs/mesh_arrays.npz')
    assert (H/'run/PETSC_OPTIONS.txt').read_bytes()==(O/'run/PETSC_OPTIONS.txt').read_bytes()
    assert (H/'policy.json').read_bytes()==(O/'policy.json').read_bytes()
    x,t,b,tags,u,p=[d[k] for k in ['points_m','tetra','boundary_triangles','facet_tags','velocity_m_s','pressure_pa']]
    tree=ET.parse(H/'run/solver.xml');mu=float(tree.find('.//Viscosity/Value').text);dt=float(tree.find('.//Time_step_size').text)
    for key in ['points_m','tetra','boundary_triangles','facet_tags']:
        assert np.array_equal(d[key],m[key]) and np.array_equal(d[key],old[key]),key
    mesh_comparisons={}
    for label,mp in [('current_mesh',C/'outputs/mesh_and_flow/solver_mesh/mesh_arrays.npz'),('viz_mesh',V/'input_data/solver_mesh/mesh_arrays.npz')]:
        if mp.exists():
            record(mp); mm=np.load(mp);mesh_comparisons[label]={k:np.array_equal(d[k],mm[k]) for k in ['points_m','tetra','boundary_triangles','facet_tags']}
    wallids=np.flatnonzero(tags==1);wall=b[wallids]; own=prod.boundary_owners(t,wall)
    center,area,n=prod.wall_geometry(x,t,wall,own); grad=prod.p1_gradients(x,t,u)
    tr=prod.tangential_traction(grad[own],n,mu);wss=np.linalg.norm(tr,axis=1)
    disp,weight=prod.nodal_average(wall,wss,area,len(x)); ws,ids=prod.surface(x,wall)
    original=pv.read(A/'inputs/display/wall_wss_si.vtp')
    idcheck=np.array_equal(original['GlobalNodeID_zero_based'],ids)
    assert idcheck and np.array_equal(original.points,x[ids])
    assert np.array_equal(original.faces,ws.faces)
    rawerr=float(np.max(abs(original['WSS_raw_Pa']-wss))); diserr=float(np.max(abs(original['WSS_display_Pa']-disp[ids])))
    assert rawerr==0 and diserr==0 and np.array_equal(original['Parent_tetra_zero_based'],own)
    assert np.array_equal(original['Global_boundary_facet_zero_based'],wallids)
    affine=np.concatenate([np.ones((len(t),4,1)),x[t]],axis=2)
    independent=np.linalg.solve(affine,u[t])[:,1:].swapaxes(1,2)
    gerr=float(abs(independent-grad).max())
    ws.cell_data['WSS_raw_Pa']=wss;ws.point_data['WSS_display_Pa']=disp[ids]
    ws.cell_data['parent_tetra_zero_based']=own;ws.cell_data['boundary_facet_zero_based']=wallids
    ws.save(A/'data/wall_audit.vtp')
    policy=json.loads((H/'policy.json').read_text())
    measure=SolutionMeasurements(A/'inputs/mesh_arrays.npz',policy['Q_target_m3_s'],policy['Umean_m_s'])
    flux=measure.measure(u,p);fluxold=measure.measure(old['velocity_m_s'],old['pressure_pa'])
    # VTU, NPZ, display bundle, frozen manifest and native checkpoint all cross-checked.
    correspondence={}
    for label,vtu in [('frozen',H/'frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu'),('run71',H/'run/1-procs/result_071.vtu'),('bundled',V/'input_data/frozen_flow/steady_flow_mean_2p0_mmps.vtu')]:
        record(vtu);grid=pv.read(vtu);uu,pp=measure.read(vtu)
        vt=grid.cells.reshape(-1,5)[:,1:]
        same_nodes=np.array_equal(np.sort(vt,axis=1),np.sort(t,axis=1))
        assert same_nodes, 'VTU element IDs do not correspond to NPZ tetrahedra'
        correspondence[label]=dict(sha256=sha(vtu),points_exact=np.array_equal(grid.points,x),tetra_local_order_exact=np.array_equal(vt,t),tetra_same_nodes_same_cell_order=same_nodes,
            local_order_note='VTU exchanges the first two vertices relative to positively oriented NPZ; same element IDs and node sets.',velocity_max_abs_m_s=abs(uu-u).max(),pressure_max_abs_Pa=abs(pp-p).max())
    # Native checkpoint Y_n precedes A_n; establish layout from comparison, not assumption.
    cp=H/'run/1-procs/stFile_071.bin';record(cp)
    values=np.frombuffer(cp.read_bytes(),dtype='<f8',offset=56).reshape(2,len(x),4)
    correspondence['checkpoint71']=dict(velocity_max_abs_m_s=abs(values[0,:,:3]-u).max(),pressure_max_abs_Pa=abs(values[0,:,3]-p).max())
    assert correspondence['checkpoint71']['velocity_max_abs_m_s']<1e-14
    assert correspondence['checkpoint71']['pressure_max_abs_Pa']<1e-7
    oldgrad=prod.p1_gradients(x,t,old['velocity_m_s']); oldw=np.linalg.norm(prod.tangential_traction(oldgrad[own],n,mu),axis=1)
    # Signed volume/Jacobian and documented signed inverse Frobenius condition number.
    det=np.linalg.det(x[t[:,1:]]-x[t[:,:1]]);vol=det/6
    sicn=tetra_sicn(x,t); hnear=3*vol[own]/area
    edge=np.linalg.norm(x[t[:,[0,0,0,1,1,2]]]-x[t[:,[1,2,3,2,3,3]]],axis=2)
    div=np.trace(grad,axis1=1,axis2=2); gnorm=np.linalg.norm(grad,axis=(1,2))
    wallflags=np.zeros(len(x),bool);wallflags[ids]=True
    allwall=wallflags[t].sum(axis=1)==4
    faces=t[:,[[0,1,2],[0,1,3],[0,2,3],[1,2,3]]].reshape(-1,3)
    _,_,count=np.unique(prod.face_keys(faces),return_index=True,return_counts=True)
    meshstats=dict(nodes=len(x),tetra=len(t),wall_facets=len(wall),boundary_facets=len(b),nonpositive_volume=int((vol<=0).sum()),min_volume_um3=vol.min()*1e18,
        quality_min=sicn.min(),quality_p01=np.quantile(sicn,.01),quality_p05=np.quantile(sicn,.05),quality_median=np.median(sicn),quality_recomputed_max_error=abs(sicn-m['min_sicn']).max(),
        wall_owner_quality_min=sicn[own].min(),edge_length_um_quantiles=np.quantile(edge,[0,.05,.5,.95,1])*1e6,
        nearwall_altitude_um_quantiles=qtile(hnear,area,(0,.05,.5,.95,1))*1e6,
        face_multiplicity_above2=int((count>2).sum()),exterior_count=int((count==1).sum()),
        tetra_all_four_nodes_on_wall=int(allwall.sum()),volume_rms_div_s_inv=np.sqrt(np.average(div**2,weights=vol)),
        rms_div_over_rms_grad=np.sqrt(np.dot(vol,div**2)/np.dot(vol,gnorm**2)),max_abs_div_s_inv=abs(div).max())
    # Exact ROI graph topology identifies the first junction downstream of the inlet.
    gd=np.load(A/'inputs/network/analysis_A_H0_graph_si.npz');edges=gd['edges'][gd['roi_internal_edge_mask']];xyz=gd['xyz_m'];gnodes=gd['ids']
    ports=json.loads((A/'inputs/network/roi_ports_in_a.json').read_text())['ports']
    adj={}
    for i,j in edges:adj.setdefault(int(i),[]).append(int(j));adj.setdefault(int(j),[]).append(int(i))
    inlet=int(np.flatnonzero(gnodes==-10001)[0]); parent={inlet:-1};dist={inlet:0.};queue=[inlet]
    for i in queue:
        for j in adj[i]:
            if j not in parent:parent[j]=i;dist[j]=dist[i]+np.linalg.norm(xyz[i]-xyz[j]);queue.append(j)
    junctions=sorted([i for i in adj if len(adj[i])>=3],key=lambda i:dist[i]); j1=junctions[0]
    regions={'all_wall':np.ones(len(wall),bool),'first_junction_r5um':np.linalg.norm(center-xyz[j1],axis=1)<5e-6}
    region_defs={'first_junction_r5um':dict(center_um=xyz[j1]*1e6,radius_um=5,graph_node_id=int(gnodes[j1]),path_from_real_inlet_um=dist[j1]*1e6),
        'junctions':[dict(node_id=int(gnodes[j]),xyz_um=xyz[j]*1e6,path_from_real_inlet_um=dist[j]*1e6) for j in junctions]}
    for port in ports:
        cp=np.array(port['fem_cap']['centroid_um'])*1e-6
        regions[port['name']+'_cap_r5um']=np.linalg.norm(center-cp,axis=1)<5e-6
        origin=np.array(port['real_cut_xyz_um'])*1e-6; normal=np.array(port['outward_normal']);s=(center-origin)@normal
        radial=np.linalg.norm(center-origin-s[:,None]*normal,axis=1)
        regions[port['name']+'_extension']= (s>0)&(s<port['fem_cap']['real_cut_to_cap_axial_um']*1e-6)&(radial<3*port['radius_um']*1e-6)
        region_defs[port['name']]=dict(cap_center_um=cp*1e6,real_cut_um=origin*1e6,outward_normal=normal,radius_um=port['fem_cap']['equivalent_radius_um'])
    regions['raw_below_2Pa']=wss<2;regions['raw_above_30Pa']=wss>30
    regionrows=[]
    for name,mask in regions.items():
        for case,v in [('H0',wss),('zero_outlets',oldw)]:
            regionrows.append(dict(case=case,region=name,facet_count=mask.sum(),area_um2=area[mask].sum()*1e12,**stats(v[mask],area[mask]),
                owner_sicn_min=sicn[own[mask]].min(),owner_sicn_median=np.median(sicn[own[mask]]),nearwall_altitude_mean_um=np.average(hnear[mask],weights=area[mask])*1e6))
    csvout('data/region_summary.csv',regionrows);dump('data/regions.json',region_defs)
    # Neighbor discontinuities on original triangles (no display processing).
    we=wall[:,[[0,1],[1,2],[2,0]]].reshape(-1,2);we=np.sort(we,axis=1)
    _,inverse,counts=np.unique(we,axis=0,return_inverse=True,return_counts=True)
    order=np.argsort(inverse,kind='stable');starts=np.r_[0,np.cumsum(counts)[:-1]];ind=starts[counts==2]
    left=order[ind]//3;right=order[ind+1]//3;jump=abs(wss[left]-wss[right]);ang=np.rad2deg(np.arccos(np.clip(np.sum(n[left]*n[right],axis=1),-1,1)))
    adjacencyrows=[]
    for ii in np.argsort(jump)[-30:][::-1]:
        l,r=left[ii],right[ii];pt=.5*(center[l]+center[r])
        adjacencyrows.append(dict(facet_a_zero_based=wallids[l],facet_b_zero_based=wallids[r],parent_a_zero_based=own[l],parent_b_zero_based=own[r],x_um=pt[0]*1e6,y_um=pt[1]*1e6,z_um=pt[2]*1e6,wss_a_Pa=wss[l],wss_b_Pa=wss[r],jump_Pa=jump[ii],normal_angle_deg=ang[ii],sicn_a=sicn[own[l]],sicn_b=sicn[own[r]]))
    csvout('data/largest_adjacent_jumps.csv',adjacencyrows)
    meshstats['adjacent_wss_jump_Pa_quantiles']=np.quantile(jump,[.5,.95,.99,1]);meshstats['adjacent_normal_angle_deg_quantiles']=np.quantile(ang,[.5,.95,.99,1])
    meshstats['adjacent_below2_above20_pairs']=int(((np.minimum(wss[left],wss[right])<2)&(np.maximum(wss[left],wss[right])>20)).sum())
    # All facets, explicit IDs and units. Display-at-centroid means linear interpolation of displayed nodal magnitude.
    rows=[]
    for i in range(len(wall)):
        rows.append(dict(case='H0_step71',facet_zero_based=wallids[i],parent_tetra_zero_based=own[i],x_um=center[i,0]*1e6,y_um=center[i,1]*1e6,z_um=center[i,2]*1e6,area_um2=area[i]*1e12,wss_raw_Pa=wss[i],display_centroid_Pa=disp[wall[i]].mean(),old_bc_wss_Pa=oldw[i],sicn=sicn[own[i]],nearwall_altitude_um=hnear[i]*1e6,owner_max_edge_um=edge[own[i]].max()*1e6,owner_div_s_inv=div[own[i]]))
    csvout('data/wall_facets.csv',rows)
    # O2 path from first junction to real cut, then artificial extension; nearest segment projection.
    outlet=int(np.flatnonzero(gnodes==-10003)[0]);path=[outlet]
    while path[-1]!=j1:path.append(parent[path[-1]])
    path=path[::-1];pathxyz=np.vstack([xyz[path],np.array(region_defs['O2']['cap_center_um'])*1e-6]);pathradius=np.r_[gd['radius_m'][path],region_defs['O2']['radius_um']*1e-6]
    seg=pathxyz[1:]-pathxyz[:-1];length=np.linalg.norm(seg,axis=1);arc=np.r_[0,np.cumsum(length)];sbest=np.zeros(len(center));rbest=np.zeros(len(center));dbest=np.full(len(center),np.inf)
    for i,(start,vec,ll) in enumerate(zip(pathxyz[:-1],seg,length)):
        f=np.clip((center-start)@vec/ll**2,0,1);distance=np.linalg.norm(center-start-f[:,None]*vec,axis=1);sel=distance<dbest
        sbest[sel]=arc[i]+f[sel]*ll;rbest[sel]=(1-f[sel])*pathradius[i]+f[sel]*pathradius[i+1];dbest[sel]=distance[sel]
    # Radius bound prevents remote branches from being assigned to O2; bins near the junction can contain blended geometry.
    pathmask=dbest<2.5*rbest;prof=[]
    for lo in np.arange(0,arc[-1]*1e6,1):
        ma=pathmask&(sbest*1e6>=lo)&(sbest*1e6<lo+1)
        if not ma.any():continue
        for name,v in [('H0',wss),('zero_outlets',oldw)]:
            prof.append(dict(case=name,s_mid_um=lo+.5,facet_count=ma.sum(),area_um2=area[ma].sum()*1e12,**stats(v[ma],area[ma]),radius_centerline_mean_um=np.mean(rbest[ma])*1e6,nearwall_altitude_mean_um=np.average(hnear[ma],weights=area[ma])*1e6))
    csvout('data/outlet2_path_profile.csv',prof)
    dump('data/outlet2_path.json',dict(node_ids=[int(gnodes[j]) for j in path],path_xyz_um=pathxyz*1e6,arc_um=arc*1e6,real_cut_arc_um=arc[-2]*1e6,filter='nearest path segment and distance < 2.5 interpolated SWC radius; first junction bins may overlap another branch'))
    # Reuse saved time series, no new CFD solve.
    times=[]
    for step in [10,20,30,40,50,60,70,71]:
        vp=H/f'run/1-procs/result_{step:03d}.vtu';record(vp)
        ut,pt=measure.read(vp);wt=np.linalg.norm(prod.tangential_traction(prod.p1_gradients(x,t,ut)[own],n,mu),axis=1)
        times.append(dict(case='H0_saved_time_series',step=step,time_s=step*dt,velocity_L2_relative_to71=measure.velocity_l2(ut-u)/measure.velocity_l2(u),wss_area_L2_relative_to71=np.sqrt(np.dot(area,(wt-wss)**2)/np.dot(area,wss**2)),wss_max_abs_delta_Pa=abs(wt-wss).max(),wss_area_mean_Pa=np.average(wt,weights=area)))
    csvout('data/time_convergence.csv',times)
    bcs=[]
    for case,trr,ff in [('H0',tree,flux),('zero_outlets',ET.parse(O/'run/solver.xml'),fluxold)]:
        for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
            bc=next(z for z in trr.findall('.//Add_BC') if z.attrib['name']==role)
            bcs.append(dict(case=case,boundary=role,type=bc.findtext('Type'),prescribed_value=float(bc.findtext('Value')),prescribed_unit='m3/s outward' if role=='INLET' else 'Pa',Q_outward_m3_s=ff['signed_outward_boundary_flows_m3_s'][role],flow_uL_min=abs(ff['signed_outward_boundary_flows_m3_s'][role])*6e10,flow_relative_to_inlet_pct=abs(ff['signed_outward_boundary_flows_m3_s'][role])/ff['Q_in_m3_s']*100,area_average_pressure_Pa=ff['area_average_pressure_pa'][role]))
    csvout('data/boundary_values.csv',bcs)
    # Canonical comparison excludes path spelling, but preserves every physical setting.
    ta=ET.parse(O/'run/solver.xml');tb=ET.parse(H/'run/solver.xml');xmlchanges=[]
    la=list(ta.getroot().iter());lb=list(tb.getroot().iter());assert len(la)==len(lb)
    for aa,bb in zip(la,lb):
        assert aa.tag==bb.tag and aa.attrib==bb.attrib
        if (aa.text or '').strip()!=(bb.text or '').strip():xmlchanges.append(dict(tag=aa.tag,old=(aa.text or '').strip(),new=(bb.text or '').strip()))
    dump('data/xml_changes.json',xmlchanges)
    hist=prod.parse_solver_log((H/'run/solver.log').read_text(),dt)
    solves=hist['linear_solves'];opts=(H/'run/PETSC_OPTIONS.txt').read_text().split();rtol=float(opts[opts.index('-ksp_rtol')+1]);atol=float(opts[opts.index('-ksp_atol')+1])
    residualrows=[]
    for r in solves:
        mm=r['petsc_monitor'];first,last=mm[0],mm[-1];tol=max(atol,rtol*first['true_residual_norm'])
        residualrows.append(dict(step=r['step'],nonlinear_iteration=r['nonlinear_iteration'],nonlinear_Ri_R0=r['nonlinear_Ri_over_R0'],linear_iterations=r['linear_iterations'],reason=r['petsc_reason']['reason'],final_true_residual=last['true_residual_norm'],effective_tolerance=tol,true_over_tolerance=last['true_residual_norm']/tol))
    csvout('data/residual_solves.csv',residualrows)
    residualsummary={k:hist[k] for k in ['unparsed_rows','failed_linear_solves','recovered_attempts','unassigned_petsc_monitor','nonlinear_failure_messages']}
    residualsummary.update(linear_solves=len(solves),max_true_over_tolerance=max(r['true_over_tolerance'] for r in residualrows),final_nonlinear_Ri_R0=residualrows[-1]['nonlinear_Ri_R0'])
    validation=dict(mu_Pa_s=mu,dt_s=dt,raw_reproduction_max_error_Pa=rawerr,display_reproduction_max_error_Pa=diserr,independent_all_tetra_gradient_max_error_s_inv=gerr,
        tangency_max_error_Pa=abs((tr*n).sum(1)).max(),normal_length_max_error=abs(np.linalg.norm(n,axis=1)-1).max(),
        wall_max_speed_m_s=np.linalg.norm(u[ids],axis=1).max(),display_unweighted_min_Pa=disp[ids].min(),display_unweighted_max_Pa=disp[ids].max(),
        display_zero_weight_wall_nodes=int((weight[ids]<=0).sum()),raw_nonfinite=int((~np.isfinite(wss)).sum()),raw_zero_count=int((wss==0).sum()),
        raw=stats(wss,area),mesh_correspondence=mesh_comparisons,flow_correspondence=correspondence,flux_H0=flux,flux_zero_outlets=fluxold,mesh=meshstats,residuals=residualsummary,
        old_new_wss_area_L2_delta_relative_to_old=np.sqrt(np.dot(area,(wss-oldw)**2)/np.dot(area,oldw**2)))
    dump('data/baseline_validation.json',validation)
    np.savez_compressed(A/'data/plot_arrays.npz',centers_um=center*1e6,area=area,wss=wss,oldw=oldw,display=disp[ids],sicn=sicn[own],h_um=hnear*1e6,own=own,ids=ids,wall=wall,first_junction_um=xyz[j1]*1e6,first_junction_mask=regions['first_junction_r5um'],path_s_um=sbest*1e6,pathmask=pathmask)
    dump('evidence/input_hashes.json',protected)
    changed=[p for p,h in protected.items() if sha(p)!=h];assert not changed
    dump('evidence/preservation_check.json',dict(checked=len(protected),changed=changed))
    print(json.dumps(clean(validation),indent=2),flush=True)

if __name__=='__main__':main()
