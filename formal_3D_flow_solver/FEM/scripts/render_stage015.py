#!/usr/bin/env python3
"""WSL review of a rejected surface study: never fabricate a selected volume."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyvista as pv
from fem3d.cap_remesh import project
from fem3d.audit import sha256,timestamp,write_json
R=ROOT/'reports/stage01_5';O=ROOT/'outputs/stage01_5'
def read(p): return json.loads(p.read_text())
selection=read(R/'candidate_selection.json')
assert selection['selected_candidate'] is None,'This renderer documents the geometry-gate failure branch only'
planes=read(O/'baseline/cap_projection.json')['ports'];names=list(planes)
source=np.load(ROOT/'inputs/stage01/tagged_surface_si.npz')
derived={k:np.load(O/k/'surface/tagged_surface_si.npz') for k in ('candidate_A','candidate_B','candidate_C')}
qc={k:read(O/k/'qc/surface_invariants.json') for k in derived}
raw=np.load(ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz')
baseline=read(ROOT/'outputs/stage01/medium/qc/geometry_qc.json')
colors=['#be6e2d','#138878','#456eb1'];ink='#1f3448';red='#ab4054';images={}
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
def note(fig,text): fig.text(.04,.015,'Stage 1 baseline vs Stage 1.5 | No candidate selected; geometry gate failed.\n'+text,fontsize=10,color='#505b68')
def save(fig,name,role):
    fig.savefig(R/name,dpi=160,bbox_inches='tight');plt.close(fig)
    images[name]={'sha256':sha256(R/name),'selected_candidate':None,'role':role}
def unavailable(ax,message='No admissible Stage 1.5 volume'):
    ax.axis('off');ax.text(.5,.65,message,transform=ax.transAxes,ha='center',fontsize=15,color=red,weight='bold')
    ax.text(.5,.38,'All A/B/C cap surfaces failed the fixed area gate.\nVolume meshing and winner selection were blocked.\n\nNo tetra quality or round-trip result exists.\nMissing values are not zero.',transform=ax.transAxes,ha='center',fontsize=11,linespacing=1.8,color=ink)

fig,axes=plt.subplots(4,2,figsize=(10,14),sharex='row',sharey='row')
for i,name in enumerate(names):
    plane=planes[name]
    for j,data in enumerate((source,derived['candidate_A'])):
        xy,_=project(data['points_m'],plane['origin_m'],np.array(plane['basis']));tri=data['triangles'][data['facet_tags']==plane['entity_id']]
        ax=axes[i,j];ax.triplot(xy[:,0]*1e6,xy[:,1]*1e6,tri,lw=.55,color='#4e6475' if j==0 else colors[0]);ax.set_aspect('equal');ax.set_ylabel(f'{name.upper()}\nη (µm)');ax.set_xlabel('ξ (µm)')
        ax.set_title(f"{'Stage 1 original' if j==0 else 'Candidate A — REJECTED'} | {len(tri)} triangles",fontsize=11)
fig.suptitle('Are cap triangles more uniform?',fontsize=18)
fig.subplots_adjust(top=.94,bottom=.10,hspace=.48,wspace=.20)
note(fig,'A is shown only as the first predefined diagnostic example, not a winner. Same plane basis, scale and view per row.')
save(fig,'cap_triangulation_before_after.png','Actual baseline vs rejected A cap triangulations; no selected cap exists')

fig,(ax,bx)=plt.subplots(1,2,figsize=(12,7),gridspec_kw={'width_ratios':[1.6,1]})
arrays=[np.load(O/'baseline/cap_quality.npz')]+[np.load(O/k/'qc/cap_quality.npz') for k in derived]
labels=['Stage 1 original','Candidate A','Candidate B','Candidate C'];stats=[]
for d,label,c in zip(arrays,labels,['#525d6a']+colors):
    values=np.concatenate([d[name+'_q_tri'] for name in names]);v=np.sort(values)
    ax.plot(v,np.arange(1,len(v)+1)/len(v),label=label,color=c,lw=2)
    stats.append([label,f'{np.quantile(v,.05):.4f}',f'{np.median(v):.4f}',str(np.count_nonzero(v<.1))])
ax.set(xlabel='Normalized triangle quality q_tri',ylabel='Cumulative fraction',title='Full distributions across four caps',xlim=(0,1.01));ax.legend(loc='upper left');ax.grid(alpha=.2)
bx.axis('off');tbl=bx.table(cellText=stats,colLabels=['Cap set','P5','Median','q<0.1'],loc='center',cellLoc='center',colWidths=[.43,.20,.22,.20]);tbl.auto_set_font_size(False);tbl.set_fontsize(10);tbl.scale(1,2)
bx.text(.5,.21,'Selected: NONE\nBetter triangles do not override\nthe failed area invariant.',ha='center',transform=bx.transAxes,color=red,fontsize=12)
fig.suptitle('Cap triangle quality improved; geometric acceptance did not',fontsize=16);fig.subplots_adjust(bottom=.18,top=.88,wspace=.2)
note(fig,'q_tri = 4√3 × area / sum(edge length²). Baseline and A/B/C use the same 3D definition.')
save(fig,'cap_triangle_quality.png','Measured cap quality distributions for baseline and all rejected candidates')

fig,(ax,bx)=plt.subplots(1,2,figsize=(12,7));q=raw['min_sicn'];s=baseline['quality']['gmsh_min_sicn']
ax.hist(q,bins=np.linspace(0,1,61),color='#138878',alpha=.85)
for key,c in [('minimum',red),('P1','#be6e2d'),('P5','#456eb1'),('median',ink)]:
    ax.axvline(s[key],color=c,lw=1.4,label=f'{key} = {s[key]:.6f}')
ax.set(xlabel='Gmsh minSICN',ylabel='Tetra count',title='Stage 1 medium | 147,569 tetra');ax.legend(fontsize=9)
ax.text(.05,.64,'q < 0.1: 153 cells',transform=ax.transAxes,va='top',color=red)
unavailable(bx);fig.suptitle('Did tetra quality improve? Not evaluated after geometry rejection',fontsize=15);fig.subplots_adjust(bottom=.18,top=.87,wspace=.28)
note(fig,'Common histogram bins would be used for a valid winner. Only the frozen baseline exists.')
save(fig,'tetra_quality_before_after.png','Baseline histogram with explicitly unavailable Stage 1.5 comparison')

fig,(ax,bx)=plt.subplots(1,2,figsize=(12,7));boundaries=['WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03'];counts=baseline['quality']['low_quality_nearest_boundary_counts']
y=np.arange(5);values=[counts[n] for n in boundaries];ax.barh(y,values,color=['#64748b']+[colors[0]]*4)
ax.set_yticks(y,boundaries);ax.invert_yaxis();ax.set(xlabel='Tetra with minSICN < 0.1',title='Stage 1 medium: 153 total, 129 cap-adjacent',xlim=(0,62))
for i,v in enumerate(values):ax.text(v+1,i,str(v),va='center')
unavailable(bx);fig.suptitle('Did low-quality tetra counts decrease?',fontsize=17);fig.subplots_adjust(bottom=.18,top=.87,wspace=.3)
note(fig,'Classification remains nearest boundary triangle center, exactly as in Stage 1. No new counts are claimed.')
save(fig,'low_quality_count_by_boundary.png','Frozen boundary counts and explicit missing candidate values')

fig,axes=plt.subplots(2,2,figsize=(11,9))
for ax,name in zip(axes.flat,names):
    plane=planes[name];ids=np.array(plane['ccw_rim_ids']);loop=np.r_[ids,ids[0]]
    xy,_=project(source['points_m'],plane['origin_m'],np.array(plane['basis']))
    ax.plot(xy[loop,0]*1e6,xy[loop,1]*1e6,color=ink,lw=3,label='Stage 1 rim')
    for (candidate,data),c,style in zip(derived.items(),colors,['--',':','-.']):
        after,_=project(data['points_m'],plane['origin_m'],np.array(plane['basis']))
        assert np.array_equal(data['points_m'][ids],source['points_m'][ids])
        ax.plot(after[loop,0]*1e6,after[loop,1]*1e6,color=c,lw=1.3,ls=style,label=candidate.replace('candidate_','Candidate '))
    ax.set_aspect('equal');ax.set(title=name.upper(),xlabel='ξ (µm)',ylabel='η (µm)')
    ax.text(.5,.48,'Maximum rim displacement\n= 0 m (A/B/C)',ha='center',va='center',transform=ax.transAxes,color=ink,fontsize=10)
axes.flat[0].legend(fontsize=8,loc='upper right');fig.suptitle('Did the port outlines move? Exact coordinate differences are zero',fontsize=15);fig.subplots_adjust(bottom=.14,top=.91,hspace=.36,wspace=.3)
note(fig,'Original vertex IDs and complete rim edge sets are identical. No projected rim coordinates replaced the original 3D points.')
save(fig,'rim_overlay.png','Actual exact rim overlays for all three candidates; none selected')

# Actual tetra shapes and surface contexts, in display-only micrometres.
points=raw['points_m']*1e6;tri=raw['boundary_triangles'];tet=raw['tetra']
surf=pv.PolyData(points,np.c_[np.full(len(tri),3),tri].ravel());wall=surf.extract_cells(raw['facet_tags']==1).extract_surface()
grid=pv.UnstructuredGrid(np.c_[np.full(len(tet),4),tet].ravel(),np.full(len(tet),pv.CellType.TETRA,dtype=np.uint8),points)
centers=points[tet].mean(axis=1);worst=np.argsort(q)[:20];low=centers[worst[0]]
def make(shape,size):
    plot=pv.Plotter(off_screen=True,shape=shape,window_size=size,border=False)
    for r in plot.renderers:r.set_background('white')
    return plot
def pvsave(plot,name,role):
    plot.screenshot(str(R/name));plot.close();images[name]={'sha256':sha256(R/name),'selected_candidate':None,'role':role}
def blank(plot,title):
    plot.add_text(title,position='upper_left',font_size=15,color=ink)
    plot.add_text('NOT GENERATED\n\nAll A/B/C surfaces failed\nthe fixed cap area gate.\n\nNo selected volume exists.\nNo comparison is claimed.',position=(50,100),font_size=17,color=red)
plot=make((2,2),(2000,1400))
plot.subplot(0,0);plot.add_mesh(wall,color='#b4bcc5',opacity=.25);plot.add_mesh(grid.extract_cells(worst),color='#d6513d',show_edges=True)
plot.add_points(centers[worst],color='#d6513d',point_size=8,render_points_as_spheres=True)
plot.add_text('Stage 1 medium | worst 20 tetra',font_size=16,color=ink)
plot.add_text('Location markers enlarged; tetra shapes unchanged',position='lower_left',font_size=11,color=ink)
plot.view_isometric();plot.enable_parallel_projection();plot.reset_camera();plot.camera.zoom(.85)
plot.subplot(0,1);blank(plot,'Stage 1.5 | worst 20: unavailable')
plot.subplot(1,0)
local=surf.extract_cells(np.linalg.norm(surf.cell_centers().points-low,axis=1)<1.0)
plot.add_mesh(local,color='#b9c5ce',opacity=.28,show_edges=True,line_width=.6)
plot.add_mesh(grid.extract_cells([int(worst[0])]),color='#d6513d',opacity=.8,show_edges=True,line_width=2)
plot.add_text(f'Stage 1 actual worst cell #{worst[0]} | q = {q[worst[0]]:.6f}',font_size=15,color=ink)
plot.add_text('Actual relative geometry; camera zoom only\nNearest patch: INLET; no tetra enlargement',position='lower_left',font_size=11,color=ink)
plot.view_isometric();plot.enable_parallel_projection();plot.reset_camera();plot.camera.zoom(.88)
plot.subplot(1,1);blank(plot,'Stage 1.5 | worst-cell shape: unavailable')
pvsave(plot,'worst_elements_before_after.png','Actual baseline worst 20 and worst-cell shape; candidate volume panels explicitly unavailable')

plot=make((1,2),(2000,1050))
for j,data in enumerate((source,derived['candidate_A'])):
    plot.subplot(0,j);p=data['points_m']*1e6;t=data['triangles'];tags=data['facet_tags']
    surface=pv.PolyData(p,np.c_[np.full(len(t),3),t].ravel())
    plot.add_mesh(surface.extract_cells(tags==1),color='#a4adb5')
    pos=[];labels=[]
    for name,plane in planes.items():
        color=colors[0] if name=='inlet' else colors[1]
        plot.add_mesh(surface.extract_cells(tags==plane['entity_id']),color=color,lighting=False)
        center=np.array(plane['origin_m'])*1e6;n=np.array(plane['basis'])[2]
        plot.add_mesh(pv.Arrow(start=center,direction=n,scale=9),color=color)
        pos.append(center+n*13);labels.append(name.upper())
    plot.add_point_labels(np.array(pos),labels,font_size=18,text_color=ink,shape_color='white',shape_opacity=.95,show_points=False,always_visible=True)
    plot.add_text('Stage 1 original boundary' if j==0 else 'Stage 1.5 candidate A | REJECTED, not selected',font_size=15,color=ink if j==0 else red)
    plot.add_text('1 inlet, 3 outlets, same frozen wall\nOriginal contract normals; display coordinates: µm',position='lower_left',font_size=11,color=ink)
    plot.view_isometric();plot.enable_parallel_projection();plot.reset_camera();plot.camera.zoom(.82)
pvsave(plot,'boundary_tags_selected.png','Source vs rejected A boundary tags; filename retained for requested review list, no selected surface claimed')

plot=make((1,2),(2000,1050));plot.subplot(0,0)
origin=grid.center;normal=np.array([1.,.25,.12]);normal/=np.linalg.norm(normal)
clipped=grid.clip(normal=normal,origin=origin);sliced=grid.slice(normal=normal,origin=origin)
plot.add_mesh(clipped,color='#acc1c7',show_edges=False,opacity=.8)
plot.add_mesh(sliced,color='#e8bc7f',show_edges=True,edge_color='#6f777c',line_width=.55)
plot.add_mesh(wall,color='#b6c0c9',opacity=.10)
plot.add_text('Stage 1 medium | actual interior cutaway',font_size=16,color=ink)
plot.add_text('Only cells intersecting the cut show edges\nExisting baseline, not a Stage 1.5 mesh',position='lower_left',font_size=12,color=ink)
plot.view_isometric();plot.enable_parallel_projection();plot.reset_camera();plot.camera.zoom(.85)
plot.subplot(0,1);blank(plot,'Stage 1.5 selected cutaway: unavailable')
pvsave(plot,'tetrahedral_cutaway_selected.png','Existing baseline interior cutaway; no candidate volume or fill result exists')
write_json(R/'visualization_manifest.json',{'timestamp':timestamp(),'study_status':'FAIL','selected_candidate':None,'preview_candidate':'candidate_A, first predefined candidate, for diagnostic illustration only','limitation':'Four volume-dependent selected comparisons cannot be produced. Their requested filenames contain real baseline evidence and clearly labelled unavailable panels, not invented new meshes.','images':images,'source_sha256':{str(p.relative_to(ROOT)):sha256(p) for p in [ROOT/'inputs/stage01/tagged_surface_si.npz',ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz']+[O/k/'surface/tagged_surface_si.npz' for k in derived]}})
print('Saved 8 diagnostic figures; selected-volume comparisons explicitly unavailable')
