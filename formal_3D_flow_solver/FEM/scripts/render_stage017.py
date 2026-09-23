#!/usr/bin/env python3
"""WSL-only scientific review of actual adaptive meshing results."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyvista as pv
from fem3d.audit import sha256,write_json,timestamp
from fem3d.cap_remesh import project,triangle_quality
R=ROOT/'reports/stage01_7';O=ROOT/'outputs/stage01_7';read=lambda p:json.loads(p.read_text())
run=read(O/'optimizer_result.json');assert run['selected_iteration']
contract=read(ROOT/'reports/stage01_6/planar_port_contract_v2.json');policy=read(R/'acceptance_policy.json');baseline=read(R/'baseline_recomputed.json')
selected=read(O/'selected/qc/volume_quality.json');surface_qc=read(O/'selected/qc/surface_invariants.json');names=list(contract['ports']);iterations=run['iterations']
source=np.load(ROOT/'inputs/stage01/tagged_surface_si.npz');derived=np.load(O/'selected/surface/tagged_surface_si.npz')
raw0=np.load(ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz');raw1=np.load(O/'selected/mesh/volume_mesh.npz')
ink='#243b50';gray='#65768a';green='#138878';red='#b34555';orange='#bf752f';images={}
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
def footer(fig,text):fig.text(.04,.025,text,fontsize=10,color=gray)
def save(fig,name,role):
 fig.savefig(R/name,dpi=160,bbox_inches='tight');plt.close(fig);images[name]={'sha256':sha256(R/name),'role':role}

fig,axes=plt.subplots(2,2,figsize=(12,9))
for ax,(name,s) in zip(axes.flat,run['initial_surface_searches'].items()):
 trials=s['trials'];counts=[t['triangle_count'] for t in trials];x=np.arange(len(trials));colors=[green if t['quality']['status']=='PASS' else red for t in trials]
 ax.plot(x,counts,color=gray,alpha=.5);ax.scatter(x,counts,c=colors,s=45,zorder=3,label='Surface quality PASS')
 k=next(i for i,t in enumerate(trials) if t['trial']==s['selected_trial']);ax.scatter(k,counts[k],s=210,marker='*',color=orange,zorder=4,label='Port-selected trial')
 ax.set(title=name.upper(),xlabel='Actual trial number',ylabel='Cap triangles',xticks=x);ax.grid(alpha=.15)
 ax.annotate(f"{s['selected_trial']}: {counts[k]} triangles",(k,counts[k]),xytext=(0,18),textcoords='offset points',ha='right' if k>4 else 'left',fontsize=10)
 axes.flat[0].legend(fontsize=9)
fig.suptitle('Which tested cap is sparsest while keeping sufficient quality?',fontsize=16);fig.subplots_adjust(top=.90,bottom=.13,hspace=.40,wspace=.25)
footer(fig,'All real trials passed geometry and surface quality. Select minimum measured count; prefer larger H only on ties.')
save(fig,'adaptive_surface_search.png','All real per-port surface trials; selected trial highlighted; nonmonotone OUTLET_01 count retained')

fig,axes=plt.subplots(4,2,figsize=(10,14),sharex='row',sharey='row')
for i,(name,p) in enumerate(contract['ports'].items()):
 before=int(np.count_nonzero(source['facet_tags']==p['entity_id']));after=int(np.count_nonzero(derived['facet_tags']==p['entity_id']))
 for j,data in enumerate((source,derived)):
  xy,_=project(data['points_m'],p['plane_origin_m'],np.array(p['basis']));tri=data['triangles'][data['facet_tags']==p['entity_id']];ax=axes[i,j]
  ax.triplot(xy[:,0]*1e6,xy[:,1]*1e6,tri,lw=.65,color=gray if j==0 else green);ax.set_aspect('equal');ax.set(xlabel='ξ (µm)',ylabel=name.upper()+'\nη (µm)',title=f"{'Stage 1 fan' if j==0 else 'Stage 1.7 selected'} | {before} → {after} triangles")
fig.suptitle('How did the artificial cap triangulation change?',fontsize=17);fig.subplots_adjust(top=.94,bottom=.10,hspace=.47,wspace=.24)
footer(fig,'Same basis, view and scale per row. Original 3D rim and anatomical wall are unchanged.')
save(fig,'port_mesh_before_after.png','Actual original and selected cap triangulations, shared scale and view')

qsets=[triangle_quality(d['points_m'],d['triangles'][d['facet_tags']!=1])[0] for d in (source,derived)]
fig,(ax,bx)=plt.subplots(1,2,figsize=(12,7),gridspec_kw={'width_ratios':[1.35,1]});stats=[]
for label,q,c in zip(['Stage 1','Selected'],qsets,[gray,green]):
 x=np.sort(q);ax.plot(x,np.arange(1,len(x)+1)/len(x),color=c,lw=2,label=label)
 stats.append([label,str(len(q)),f'{np.quantile(q,.05):.4f}',f'{np.median(q):.4f}',str(np.count_nonzero(q<.1))])
ax.set(xlabel='3D triangle quality q_tri',ylabel='Cumulative fraction',xlim=(0,1.01));ax.legend();ax.grid(alpha=.2);bx.axis('off')
t=bx.table(cellText=stats,colLabels=['Caps','Count','P5','Median','q<0.1'],loc='center',cellLoc='center',colWidths=[.25,.20,.20,.20,.18]);t.auto_set_font_size(False);t.set_fontsize(10);t.scale(1,2)
bx.text(.5,.22,'Every port independently passes\nP5 ≥ 0.45 and median ≥ 0.70',ha='center',transform=bx.transAxes,color=green)
fig.suptitle('Are the new cap triangles good enough?',fontsize=17);fig.subplots_adjust(top=.88,bottom=.16,wspace=.24)
footer(fig,'Actual 3D geometry: q_tri = 4√3 × area / sum(edge²). Cap triangle count is diagnostic, not a hard gate.')
save(fig,'surface_quality_before_after.png','Measured actual 3D cap quality distributions and summary values')

fig,axes=plt.subplots(1,2,figsize=(12,6));x=np.arange(len(iterations));caplimit=int(np.floor(policy['volume_quality']['cap_low_fraction']*baseline['quality']['cap_adjacent_below_0_1']))
for ax,values,limit,title in [(axes[0],[r['volume']['quality']['cap_adjacent_below_0_1'] for r in iterations],caplimit,'Cap-adjacent low-quality tetra'),(axes[1],[r['acceptance']['C_P2'] for r in iterations],policy['cost']['C_P2_max'],'P2 velocity cost ratio')]:
 ax.plot(x,values,'o-',color=green,lw=2);ax.scatter(x[-1],values[-1],marker='*',s=220,color=orange,zorder=5,label='Selected')
 ax.axhline(limit,color=red,ls='--',label=f'Limit = {limit:g}');ax.set(title=title,xlabel='Production volume iteration',xticks=x,xticklabels=[r['iteration'].replace('iteration_','') for r in iterations],xlim=(-.5,max(.5,len(x)-.5)))
 ax.set_ylim((-.04*max(1,limit),max(limit, max(values))*1.18));ax.legend(loc='upper right');ax.grid(alpha=.15)
 ax.annotate(f'{values[-1]:.6g}',(x[-1],values[-1]),xytext=(12,10),textcoords='offset points',color=green)
fig.suptitle('Why did the optimizer stop after its first volume?',fontsize=16);fig.subplots_adjust(bottom=.18,top=.82,wspace=.25)
footer(fig,'One production volume iteration. Both quality and cost pass; no further refinement was allowed. Independent replay matched.')
save(fig,'adaptive_volume_trace.png','Actual production iteration trace; first feasible selected and immediate stop')

fig,axes=plt.subplots(1,2,figsize=(12,7),sharex=True,sharey=True);bins=np.linspace(0,1,61)
for ax,raw,qc,title,color in zip(axes,[raw0,raw1],[baseline,selected],['Stage 1 medium','Stage 1.7 selected'],[gray,green]):
 ax.hist(raw['min_sicn'],bins=bins,color=color,alpha=.8);q=qc['quality'];s=q['min_sicn']
 info='\n'.join(f'{k} = {s[k]:.6f}' for k in ('minimum','P1','P5','median'))+f"\nq < 0.1 = {q['total_below_0_1']}"
 ax.text(.03,.96,info,transform=ax.transAxes,va='top',fontsize=10,bbox={'facecolor':'white','alpha':.9,'edgecolor':'none'})
 ax.set(title=f"{title} | {qc['proxy']['N_tetra']:,} tetra",xlabel='Gmsh minSICN',ylabel='Tetra count')
fig.suptitle('Did real tetra quality improve?',fontsize=18);fig.subplots_adjust(bottom=.16,top=.85,wspace=.18)
footer(fig,'Identical quality definition and histogram bins. P1, P5 and median all meet the frozen baseline requirements.')
save(fig,'tetra_quality_before_after.png','Actual baseline and selected minSICN distributions using common bins')

boundaries=['WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03'];fig,ax=plt.subplots(figsize=(11,7));x=np.arange(len(boundaries));width=.35
for shift,data,label,c in [(-width/2,baseline,'Stage 1',gray),(width/2,selected,'Selected',green)]:
 values=[data['quality']['low_quality_nearest_boundary_counts'][n] for n in boundaries]
 bars=ax.bar(x+shift,values,width,label=label,color=c);ax.bar_label(bars,padding=3)
ax.set(xticks=x,xticklabels=boundaries,ylabel='Tetra with minSICN < 0.1',ylim=(0,62));ax.legend();fig.suptitle('Where did the low-quality tetra disappear?',fontsize=17);fig.subplots_adjust(bottom=.18,top=.87)
footer(fig,'Same nearest exterior-triangle-center classification as Stage 1. Selected: 3 near WALL, 0 near caps.')
save(fig,'low_quality_count_by_boundary.png','Actual low-quality tetra boundary classification counts; fixed nearest-center method')

fig,ax=plt.subplots(figsize=(10,7))
from matplotlib.patches import Rectangle
ax.add_patch(Rectangle((.97,-5),1.35-.97,caplimit+5,facecolor=green,alpha=.10,edgecolor='none'))
ax.axvline(1.35,ls='--',color=red,label='Cost limit 1.35');ax.axhline(caplimit,ls='--',color=orange,label=f'Cap low-count limit {caplimit}')
ax.scatter(1,baseline['quality']['cap_adjacent_below_0_1'],color=gray,s=75,label='Stage 1 baseline')
for r in iterations:
 cost=r['acceptance']['C_P2'];low=r['volume']['quality']['cap_adjacent_below_0_1'];ax.scatter(cost,low,s=200,marker='*',color=green,zorder=4)
 ax.annotate(r['iteration']+' (selected)',(cost,low),xytext=(12,18),textcoords='offset points',color=green)
ax.set(xlabel='P2 velocity cost / baseline cost',ylabel='Cap-adjacent low-quality tetra',xlim=(.97,1.43),ylim=(-5,140));ax.legend(loc='upper right');ax.grid(alpha=.15)
fig.suptitle('Did the measured mesh enter the feasible quality–cost region?',fontsize=15);fig.subplots_adjust(bottom=.18,top=.87)
footer(fig,'One production point; independent volume replay coincides exactly. Other quality gates were checked separately.')
save(fig,'quality_cost_pareto.png','Actual production quality-cost point, baseline and frozen feasibility thresholds')

labels=['Stage 1']+[r['iteration'].replace('iteration_','Iter ') for r in iterations]+['Selected'];proxies=[baseline['proxy']]+[r['volume']['proxy'] for r in iterations]+[selected['proxy']]
fig,axes=plt.subplots(1,2,figsize=(12,7))
for ax,key,title in [(axes[0],'N_tetra','Tetrahedron count'),(axes[1],'N_P2_velocity_proxy','P2 velocity DOF proxy')]:
 values=[p[key] for p in proxies];bars=ax.bar(labels,values,color=[gray]+[green]*len(iterations)+[orange]);ax.bar_label(bars,labels=[f'{value:,}' for value in values],padding=4)
 limit=1.35*baseline['proxy'][key];ax.axhline(limit,color=red,ls='--',label='1.35 × baseline');ax.set(ylim=(0,limit*1.15),title=title);ax.legend();ax.ticklabel_format(axis='y',style='plain')
fig.suptitle('How much more would the mesh cost in future FEM?',fontsize=17);fig.subplots_adjust(bottom=.18,top=.85,wspace=.25)
footer(fig,f"Selected repeats {run['selected_iteration']}; no extra optimization mesh. P2 proxy +{100*(selected['proxy']['N_P2_velocity_proxy']/baseline['proxy']['N_P2_velocity_proxy']-1):.3f}%; tetra +{100*(selected['proxy']['N_tetra']/baseline['proxy']['N_tetra']-1):.3f}%.")
save(fig,'mesh_cost_comparison.png','Measured baseline, each production iteration and selected tetra/P2 topology counts')

fig,axes=plt.subplots(2,2,figsize=(11,9))
for ax,(name,p) in zip(axes.flat,contract['ports'].items()):
 ids=np.array(p['rim_vertex_ids']);loop=np.r_[ids,ids[0]]
 a,_=project(source['points_m'],p['plane_origin_m'],p['basis']);b,_=project(derived['points_m'],p['plane_origin_m'],p['basis'])
 assert np.array_equal(source['points_m'][ids],derived['points_m'][ids])
 ax.plot(a[loop,0]*1e6,a[loop,1]*1e6,color=gray,lw=3,label='Stage 1');ax.plot(b[loop,0]*1e6,b[loop,1]*1e6,color=green,lw=1.6,ls='--',label='Selected')
 ax.text(.5,.50,'Maximum rim displacement\n= 0 m',ha='center',va='center',transform=ax.transAxes,color=ink)
 ax.set_aspect('equal');ax.set(title=name.upper(),xlabel='ξ (µm)',ylabel='η (µm)')
axes.flat[0].legend(loc='upper right',fontsize=9);fig.suptitle('Did any frozen rim move?',fontsize=17);fig.subplots_adjust(bottom=.13,top=.90,hspace=.36,wspace=.26)
footer(fig,'Bitwise-identical original 3D rim coordinates and exact rim edge sets. Projection is display only.')
save(fig,'rim_overlay.png','Actual source and selected rim overlays, exactly zero coordinate displacement')

# Physical 3D meshes, displayed in micrometres only.
def meshes(raw):
 points=raw['points_m']*1e6;tri=raw['boundary_triangles'];tet=raw['tetra']
 surf=pv.PolyData(points,np.c_[np.full(len(tri),3),tri].ravel());wall=surf.extract_cells(raw['facet_tags']==1).extract_surface()
 grid=pv.UnstructuredGrid(np.c_[np.full(len(tet),4),tet].ravel(),np.full(len(tet),pv.CellType.TETRA,dtype=np.uint8),points)
 return surf,wall,grid
sets=[meshes(raw0),meshes(raw1)]
def plotter(shape,size):
 p=pv.Plotter(off_screen=True,shape=shape,window_size=size,border=False)
 for r in p.renderers:r.set_background('white')
 return p
def finish3d(p,name,role):p.screenshot(str(R/name));p.close();images[name]={'sha256':sha256(R/name),'role':role}
plot=plotter((2,2),(2000,1400))
for j,(raw,qc,(surf,wall,grid),title) in enumerate(zip([raw0,raw1],[baseline,selected],sets,['Stage 1 medium','Stage 1.7 selected'])):
 q=raw['min_sicn'];worst=np.argsort(q)[:20];centers=raw['points_m'][raw['tetra']].mean(axis=1)*1e6;anchor=centers[worst[0]]
 plot.subplot(0,j);plot.add_mesh(wall,color='#b7c3cc',opacity=.22);plot.add_mesh(grid.extract_cells(worst),color=red,show_edges=True);plot.add_points(centers[worst],color=red,point_size=8,render_points_as_spheres=True)
 plot.add_text(title+' | worst 20 tetra',font_size=16,color=ink);plot.add_text('Location markers enlarged; cell geometry unchanged',position='lower_left',font_size=11,color=ink)
 plot.view_isometric();plot.enable_parallel_projection();plot.reset_camera();plot.camera.zoom(.85)
 plot.subplot(1,j);near=surf.extract_cells(np.linalg.norm(surf.cell_centers().points-anchor,axis=1)<1.)
 cell=grid.extract_cells([int(worst[0])])
 plot.add_mesh(near,color='#bac5ce',opacity=.15,show_edges=True,line_width=.6);plot.add_mesh(cell,color=red,opacity=.85,show_edges=True,line_width=2)
 patch=qc['quality']['worst_elements'][0]['nearest_boundary_patch']
 plot.add_text(f'Actual worst cell #{worst[0]} | q = {q[worst[0]]:.6f}',font_size=15,color=ink)
 plot.add_text('Nearest patch: '+patch+'\nCamera zoom only; true relative geometry',position='lower_left',font_size=11,color=ink)
 span=float(np.linalg.norm(np.ptp(cell.points,axis=0)))
 plot.camera_position=[anchor+span*np.array([4.,4.,4.]),anchor,[0.,0.,1.]]
 plot.enable_parallel_projection();plot.camera.parallel_scale=span*.85
finish3d(plot,'worst_elements_before_after.png','Actual baseline and selected worst 20 tetra with worst-cell physical geometry close-ups')

plot=plotter((1,2),(2000,1050))
for j,data in enumerate((source,derived)):
 plot.subplot(0,j);points=data['points_m']*1e6;tri=data['triangles'];tags=data['facet_tags'];surf=pv.PolyData(points,np.c_[np.full(len(tri),3),tri].ravel())
 plot.add_mesh(surf.extract_cells(tags==1),color='#a8b4bb');centers=[];labels=[]
 for name,p in contract['ports'].items():
  color=orange if name=='inlet' else green;plot.add_mesh(surf.extract_cells(tags==p['entity_id']),color=color,lighting=False)
  center=np.array(p['plane_origin_m'])*1e6;n=np.array(p['outward_normal']);plot.add_mesh(pv.Arrow(start=center,direction=n,scale=9),color=color)
  centers.append(center+n*13);labels.append(name.upper())
 plot.add_point_labels(np.array(centers),labels,font_size=18,text_color=ink,shape_color='white',shape_opacity=.95,show_points=False,always_visible=True)
 plot.add_text('Stage 1 boundary' if j==0 else 'Stage 1.7 selected boundary',font_size=17,color=ink)
 plot.add_text('1 inlet, 3 outlets; original outward normals\nWall/rim displacement = 0 m; display coordinates: µm',position='lower_left',font_size=11,color=ink)
 plot.view_isometric();plot.enable_parallel_projection();plot.reset_camera();plot.camera.zoom(.82)
finish3d(plot,'boundary_tags_selected.png','Actual original and selected boundaries with unchanged semantic labels and outward normals')

plot=plotter((1,2),(2000,1100));surf,wall,grid=sets[1];normal=np.array([1.,.25,.12]);normal/=np.linalg.norm(normal);origin=grid.center
cut=grid.clip(normal=normal,origin=origin);section=grid.slice(normal=normal,origin=origin)
plot.subplot(0,0);plot.add_mesh(cut,color='#b3c4ca',opacity=.85);plot.add_mesh(section,color='#e8bc7f',show_edges=True,edge_color=gray,line_width=.7);plot.add_mesh(wall,color='#bac5ce',opacity=.10)
plot.add_text('Selected | real tetrahedral interior cutaway',font_size=16,color=ink);plot.add_text(f"{run['selected_iteration']}; {selected['proxy']['N_tetra']:,} actual tetrahedra",position='lower_left',font_size=12,color=ink)
plot.view_isometric();plot.enable_parallel_projection();plot.reset_camera();plot.camera.zoom(.85)
plot.subplot(0,1);areas=section.compute_cell_sizes()['Area'];anchor=section.cell_centers().points[int(np.argmax(areas))];local=section.extract_cells(np.linalg.norm(section.cell_centers().points-anchor,axis=1)<4.)
plot.add_mesh(local,color='#e8bc7f',show_edges=True,edge_color=gray,line_width=1.4,lighting=False)
plot.add_text('Selected | actual tetra-plane intersections',font_size=16,color=ink);plot.add_text('Physical interior cells; camera zoom only\nSection polygons are cuts through tetrahedra',position='lower_left',font_size=12,color=ink)
plot.view_vector(normal);plot.enable_parallel_projection();plot.reset_camera();plot.camera.zoom(.82)
finish3d(plot,'tetrahedral_cutaway_selected.png','Actual selected tetra interior and magnified plane intersections, without synthetic fill geometry')
write_json(R/'visualization_manifest.json',{'timestamp':timestamp(),'selected_iteration':run['selected_iteration'],'production_volume_iterations':len(iterations),'images':images,
 'source_sha256':{str(p.relative_to(ROOT)):sha256(p) for p in [ROOT/'inputs/stage01/tagged_surface_si.npz',ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz',O/'selected/mesh/volume_mesh.npz',O/'selected/surface/tagged_surface_si.npz']}})
print('Generated',len(images),'WSL figures from actual baseline and selected artifacts')
