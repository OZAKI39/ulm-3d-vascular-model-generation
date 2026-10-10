"""Native event/graph evidence for local shedding; never invent a continuum crack."""
import json,time
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
from .continuous_handoff_audit import OUT,RESULTS,PRIOR,PRIOR_RESULTS,CASE,quant
from .config import dump
from .y_clot_morphology_v3 import load_geometry
from .y_mural_clot_geometry import locator
from .detached_fragment_metrics import native_events,graph_components,equivalent_diameter_um
from .failure_map import csv_out
import vtk

def distances_to_triangles(mesh,X):
    loc=locator(mesh);d=[];cp=[0.,0.,0.];cid=vtk.mutable(0);sub=vtk.mutable(0);dd=vtk.mutable(0.)
    for x in X:loc.FindClosestPoint(x,cp,cid,sub,dd);d.append(float(dd)**.5)
    return np.asarray(d)

def run():
    started=time.perf_counter();g=load_geometry();n=g['clot_count'];b=g['bonds'];typ=g['bond_types'];X=g['X'];ev,bi=native_events(CASE,g);tr=np.load(CASE/'accepted_trajectory.npz');alive=np.ones(len(b),bool);alive[bi]=False
    top=list(__import__('pd_clot.lammps_fragmentation.fragment_analysis',fromlist=['read_dump']).read_dump(CASE/'surviving_bonds.local'))[-1][2]
    actual=np.sort(top[:,1:3].astype(int)-1,axis=1);N=len(X);assert np.array_equal(np.sort(actual[:,0]*N+actual[:,1]),(b[:,0]*N+b[:,1])[alive])
    exposed=pv.read(PRIOR_RESULTS/'paraview/exposed_surface.vtp');attach=pv.read(PRIOR_RESULTS/'paraview/attachment_surface.vtp');mid=X[b].mean(1);bulk=typ==1;distance=distances_to_triangles(exposed,mid[bulk]);broken_bulk=~alive[bulk];bins=np.array([0,50,100,150,250,1000])*1e-6
    strata=[]
    for lo,hi in zip(bins[:-1],bins[1:]):
        sel=(distance>=lo)&(distance<hi);count=int(sel.sum());failed=int((sel&broken_bulk).sum());strata.append(dict(depth_min_um=lo*1e6,depth_max_um=hi*1e6,initial_bulk_bonds=count,broken_bulk_bonds=failed,broken_fraction_in_stratum=failed/count if count else None,fraction_of_all_broken_bulk=failed/broken_bulk.sum()))
    csv_out(OUT/'bulk_break_depth_strata.csv',strata)
    event_rows=[];dist_ev=distances_to_triangles(exposed,ev[:,7:10]);near_attach=distances_to_triangles(attach,ev[:,7:10])
    for i,e in enumerate(ev):
        event_rows.append(dict(native_step=int(e[1]),physical_time_s=float(e[1]*.5e-6),initial_bond_index=int(bi[i]),bond_type=int(typ[bi[i]]),atom_id_1=int(e[2]),atom_id_2=int(e[3]),event_x_m=e[4],event_y_m=e[5],event_z_m=e[6],reference_x_m=e[7],reference_y_m=e[8],reference_z_m=e[9],reference_distance_exposed_m=dist_ev[i],reference_distance_attachment_m=near_attach[i]))
    csv_out(OUT/'native_break_spatial_distribution.csv',event_rows)
    it=b[typ==2];at=b[alive&(typ==2)];initial_count=np.bincount(it[it<n],minlength=n);surviving=np.bincount(at[at<n],minlength=n);interface_nodes=np.flatnonzero(initial_count)
    inode_rows=[dict(atom_id=int(i+1),x_m=X[i,0],y_m=X[i,1],z_m=X[i,2],initial_interface_bonds=int(initial_count[i]),surviving_interface_bonds=int(surviving[i]),broken_fraction=float(1-surviving[i]/initial_count[i]),all_direct_interface_bonds_lost=bool(surviving[i]==0),still_attached_through_bulk=bool(tr['attached'][-1,i])) for i in interface_nodes]
    csv_out(OUT/'interface_attachment_node_map.csv',inode_rows)
    lab,cnt,att=graph_components(n,b,typ,alive);core=int(np.argmax(np.where(att,cnt,0)));coreids=np.flatnonzero(lab==core);coremask=np.zeros(n,bool);coremask[coreids]=True
    bulk_edges=b[typ==1];core_edges=b[alive&(typ==1)];core_edges=core_edges[coremask[core_edges[:,0]]&coremask[core_edges[:,1]]];cuts=[]
    for axis,name in enumerate('xyz'):
        mn,mx=X[coreids,axis].min(),X[coreids,axis].max()
        for coordinate in np.linspace(mn+.1*(mx-mn),mx-.1*(mx-mn),15):
            aa=X[bulk_edges,axis];cc=X[core_edges,axis];initial=int(((aa.min(1)<coordinate)&(aa.max(1)>=coordinate)).sum());retained=int(((cc.min(1)<coordinate)&(cc.max(1)>=coordinate)).sum());cuts.append(dict(axis=name,coordinate_m=float(coordinate),initial_bulk_crossing_bonds=initial,final_core_surviving_crossing_bonds=retained))
    csv_out(OUT/'core_axis_section_surviving_connectivity.csv',cuts)
    # All native steps, not only the10us output grid: exact-membership lifetimes.
    active=np.ones(len(b),bool);cursor=0;previous={};records={};nextid=1;series=[]
    for step in np.unique(ev[:,1]).astype(int):
        end=np.searchsorted(ev[:,1],step,side='right');active[bi[cursor:end]]=False;step_events=ev[cursor:end];cursor=end;labels,counts,attached=graph_components(n,b,typ,active);current={}
        for j in np.flatnonzero(~attached):
            ids=np.flatnonzero(labels==j);key=tuple((ids+1).tolist());uid=previous.get(key)
            if uid is None:
                uid=nextid;nextid+=1;trigger=step_events[np.isin(step_events[:,2],ids+1)|np.isin(step_events[:,3],ids+1)]
                records[uid]=dict(episode_id=uid,node_count=len(ids),member_atom_ids=list(key),reference_mass_kg=len(ids)*1.25e-10,reference_volume_m3=len(ids)*(50e-6)**3,reference_d_eq_um=float(equivalent_diameter_um(len(ids),50e-6)),first_detached_native_step=int(step),first_detached_native_time_s=float(step*.5e-6),last_native_event_step=int(step),end_exclusive_step=4660,right_censored=True,trigger_event_positions_m=trigger[:,4:7].tolist(),reference_COM_m=X[ids].mean(0).tolist(),first_saved_COM_m=None,first_saved_COM_velocity_m_s=None,first_saved_COM_time_s=None,location_caveat='Native trigger bond event position is exact; component COM at an unsaved native break step is unavailable and is not interpolated.')
                firstframe=int(np.searchsorted(tr['step'],step))
                if firstframe<len(tr['step']):
                    ll=tr['component_id'][firstframe];valid=len(np.unique(ll[ids]))==1 and (ll==ll[ids[0]]).sum()==len(ids) and not tr['attached'][firstframe,ids].any()
                    if valid:records[uid].update(first_saved_COM_m=tr['position_m'][firstframe,ids].mean(0).tolist(),first_saved_COM_velocity_m_s=tr['velocity_m_s'][firstframe,ids].mean(0).tolist(),first_saved_COM_time_s=float(tr['step'][firstframe]*.5e-6))
            records[uid]['last_native_event_step']=int(step);current[key]=uid
        for key,uid in previous.items():
            if key not in current:records[uid]['end_exclusive_step']=int(step);records[uid]['right_censored']=False
        previous=current;free=counts[~attached];series.append(dict(step=int(step),time_s=float(step*.5e-6),detached_components=int(len(free)),largest_detached_nodes=int(free.max(initial=0)),largest_attached_nodes=int(counts[attached].max(initial=0)),detached_nodes=int(free.sum())))
    for r in records.values():r['exact_detached_membership_duration_s']=(r['end_exclusive_step']-r['first_detached_native_step'])*.5e-6
    dump(OUT/'ALL_NATIVE_DETACHED_MEMBERSHIP_EPISODES.json',dict(records=list(records.values()),scope='Every distinct native break step; no inferred COM at unsaved steps',steps_checked=len(series)))
    csv_out(OUT/'detached_component_native_lifetimes.csv',[{**r,'member_atom_ids':json.dumps(r['member_atom_ids']),'trigger_event_positions_m':json.dumps(r['trigger_event_positions_m'])} for r in records.values()]);csv_out(OUT/'native_component_timeline.csv',series)
    core_iface=at[np.any(np.isin(at,coreids),axis=1)];degree=np.bincount(core_edges.ravel(),minlength=n)
    report=dict(CURRENT_FAILURE_MODE='LOCAL_DEBRIS_SHEDDING_WITH_ATTACHED_CORE',initial_bulk_bonds=int(bulk.sum()),broken_bulk_bonds=int(broken_bulk.sum()),initial_interface_bonds=int((typ==2).sum()),broken_interface_bonds=int(((typ==2)&~alive).sum()),bulk_broken_midpoints_depth_um=quant(distance[broken_bulk]*1e6),bulk_break_fraction_within150um=float((distance[broken_bulk]<=150e-6).mean()),bulk_break_fraction_within100um=float((distance[broken_bulk]<=100e-6).mean()),depth_strata=strata,interface_endpoint_nodes=len(interface_nodes),interface_nodes_losing_all_direct_interface_bonds=int((surviving[interface_nodes]==0).sum()),interface_fully_lost_but_still_bulk_attached=int(((surviving==0)&(initial_count>0)&tr['attached'][-1]).sum()),total_attached_nodes=int(cnt[att].sum()),attached_component_sizes=sorted(cnt[att].tolist(),reverse=True),retained_core_nodes=len(coreids),retained_core_mass_fraction=len(coreids)/n,retained_core_surviving_bulk_bonds=len(core_edges),retained_core_surviving_interface_bonds=len(core_iface),retained_core_anchor_nodes=int(len(np.unique(core_iface[core_iface>=n]))),retained_core_bulk_degree=quant(degree[coreids]),axis_sections_checked=len(cuts),minimum_surviving_cross_section_bonds=min(r['final_core_surviving_crossing_bonds'] for r in cuts),native_break_steps_checked=len(series),native_detached_membership_episodes=len(records),episodes_without_saved_COM=sum(r['first_saved_COM_m'] is None for r in records.values()),largest_ever_detached_nodes=max(r['largest_detached_nodes'] for r in series),largest_ever_detached_d_eq_um=float(equivalent_diameter_um(max(r['largest_detached_nodes'] for r in series),50e-6)),macro_fragment_200um='NOT_OBSERVED',macro_fragment_300um='NOT_OBSERVED',macroscopic_crack_surface='NOT_ESTABLISHED: bond events are not an explicit continuum crack surface. No34-node or114-node region becomes an independent graph component; all tested central axis sections retain bulk paths.',macro_interface_front='NOT_ESTABLISHED: spatial node/bond maps quantify local direct-interface depletion; no gross peel-off of the anchored core.',causality='Connectivity proves remaining paths, not a unique causal mechanism or a counterfactual load/material claim.',next_candidate='Separately labelled broad-gradient tangential peeling demonstration on the same native Y geometry with fixed material; establish displacement/traction identity and monitor a coherent interface front before any run. No new solve in this task.',wall_seconds=time.perf_counter()-started)
    dump(OUT/'FRAGMENTATION_MECHANISM_DIAGNOSIS.json',report)
    (OUT/'FRAGMENTATION_MECHANISM_DIAGNOSIS.md').write_text(f'''# 原生破坏机制审查

CURRENT_FAILURE_MODE = LOCAL_DEBRIS_SHEDDING_WITH_ATTACHED_CORE。

12810次原生断键中，bulk {report['broken_bulk_bonds']} / {report['initial_bulk_bonds']}，interface {report['broken_interface_bonds']} / {report['initial_interface_bonds']}。按参考键中点到真实初始露出三角面的距离，{report['bulk_break_fraction_within100um']:.2%}的bulk断键位于100µm以内，{report['bulk_break_fraction_within150um']:.2%}位于150µm以内。各深度层同时列出全部初始键数和断键比例，避免把表面积较大误当应力集中。该距离是初始自由表面深度代理，不是变形裂面厚度。

初始有{report['interface_endpoint_nodes']}个clot节点带直接interface键；{report['interface_nodes_losing_all_direct_interface_bonds']}个节点失去全部直接interface键，其中{report['interface_fully_lost_but_still_bulk_attached']}个仍通过bulk连到锚层。节点图与三维断键分布已输出，不能把局部失去直接粘结当成整个组件脱落。未建立宏观脱粘前缘。

末态总附着1702节点（88.323819%），包括最大核心1701节点与一个仍通过interface锚定的单节点。最大核心仍含{report['retained_core_surviving_bulk_bonds']}条bulk键、{report['retained_core_surviving_interface_bonds']}条interface键并连接{report['retained_core_anchor_nodes']}个固定锚节点。三轴中央80%范围内的45个参考截面仍至少有{report['minimum_surviving_cross_section_bonds']}条核心bulk键穿过。该检查只排查指定轴向截面的完全断开，不排除任意斜向局部裂纹，更不能从散点直接重建一个已验证的连续裂面。

逐个原生断键时刻重建连通图：{len(series)}个时刻、{len(records)}个完全相同成员的游离生命周期。最大短暂游离组件13节点（约{report['largest_ever_detached_d_eq_um']:.3f}µm reference-volume equivalent diameter），末态仅一个6节点组件与207个1–2节点质量包；全程没有>=34节点/200µm或>=114节点/300µm的游离大块。{report['episodes_without_saved_COM']}个短生命周期没有任何匹配的保存帧COM，这些值明确置空；保留精确native触发断键位置、成员和时间，不插值制造首次出现位置。

因此不能把本例称作血栓碎裂为几个大块。Gaussian200Pa、当前Liu-LPS critical stretch和曲面锚定组合没有产生能隔离宏观体积的切割/剥离路径；这是当前观察与连通拓扑支持的描述，不是唯一原因证明。增大力并不保证大块，也可能增加节点级粉化。

下一阶段仅列一个候选：在同一原生Y几何和固定材料上，单独标记“宽尺度梯度切向剥离对照”，明确空间载荷、方向与界面开裂指标后再研究连贯脱粘前缘。本轮不执行该对照，不加入预裂缝、normal lift或任意额外力。
''')
    print('MECHANISM',report,flush=True)

if __name__=='__main__':run()
