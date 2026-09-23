#!/usr/bin/env python3
"""WSL-only review figures from fetched FEM fields and quadrature diagnostics."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyvista as pv
from fem3d.audit import sha256,timestamp,write_json

OUT=ROOT/'reports/stage02'
CASES=ROOT/'outputs/stage02/cases'
profiles=['coarse','medium','fine']
colors=['#de883b','#4282b1','#177d6d']
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})
def read(case,kind='analytic_comparison'):
    return json.loads((CASES/case/'qc'/f'{kind}.json').read_text())
reference=[read(f'pipe_{p}_reference') for p in profiles]
natural=[read(f'pipe_{p}_natural') for p in profiles]
images={}
def save(fig,name,description,sources):
    fig.savefig(OUT/name,dpi=160,bbox_inches='tight'); plt.close(fig)
    images[name]={'sha256':sha256(OUT/name),'description':description,'source_cases':sources}
def note(fig,text):
    fig.text(.05,.01,'SYNTHETIC NUMERICAL BENCHMARK · NOT EXPERIMENTAL CONDITION\n'+text,fontsize=10,color='#4b5563')

# Actual quadratic mesh, display coordinates only converted to micrometres.
saved=np.load(CASES/'pipe_medium_natural/solution/visualization_mesh.npz')
points=saved['points_m'][:,[2,0,1]]*1e6
mesh=pv.UnstructuredGrid(saved['topology'],saved['cell_types'],points)
mesh.point_data['speed_mm_s']=np.linalg.norm(saved['velocity_m_s'],axis=1)*1e3
surface=mesh.extract_surface(nonlinear_subdivision=1).triangulate()
tri=surface.faces.reshape(-1,4)[:,1:]
centers=surface.points[tri].mean(axis=1)
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
fig=plt.figure(figsize=(12,7.5)); ax=fig.add_subplot(111,projection='3d')
facecolors=np.full(len(tri),'#b6c5d3',dtype=object)
facecolors[np.abs(centers[:,0])<1e-8]='#e08b40'
facecolors[np.abs(centers[:,0]-50)<1e-8]='#239b89'
ax.add_collection3d(Poly3DCollection(surface.points[tri],facecolors=facecolors,edgecolors='none',alpha=.86))
ax.quiver(0,0,0,-8,0,0,color='#b45309',arrow_length_ratio=.25,linewidth=2)
ax.quiver(50,0,0,8,0,0,color='#00796b',arrow_length_ratio=.25,linewidth=2)
ax.quiver(4,0,8,10,0,0,color='#1f2937',arrow_length_ratio=.2,linewidth=2)
ax.text(5,0,9,'flow +axis',fontsize=10)
ax.set(xlim=(-9,59),ylim=(-6,6),zlim=(-6,10),xlabel='Axial s (µm)',ylabel='Local y (µm)',zlabel='Local z (µm)')
ax.set_yticks([-5,0,5]); ax.set_zticks([-5,0,5,10])
ax.xaxis.labelpad=16; ax.yaxis.labelpad=10; ax.zaxis.labelpad=10
ax.tick_params(labelsize=9)
ax.set_box_aspect((5,1,1.15)); ax.view_init(elev=21,azim=-68)
ax.set_title('Formal boundary conditions on the actual pipe mesh',pad=12)
fig.text(.06,.18,'INLET · tag 4\nOutward n points upstream\n∫ u·n dS = −Q;  σn = −λn',color='#9a4e14',fontsize=11)
fig.text(.40,.18,'WALL · tag 1\nNo slip: u = 0\nR = 5 µm; L = 50 µm',color='#374151',fontsize=11)
fig.text(.72,.18,'OUTLET · tag 2\nOutward n points downstream\nσn = 0; no pressure pin',color='#006b5c',fontsize=11)
fig.subplots_adjust(bottom=.26,top=.93,left=.02,right=.98)
note(fig,'Q > 0 is the inflow magnitude; the signed inlet flux is negative. No inlet velocity profile is prescribed.')
save(fig,'pipe_geometry_and_bc.png','Actual medium mesh with formal natural-traction BCs',['pipe_medium_natural'])

fig,axes=plt.subplots(2,3,figsize=(13,8),sharex=True,gridspec_kw={'height_ratios':[2.2,1]})
for j in range(3):
    r=np.array(reference[-1]['velocity_profiles'][j]['r_over_R'])
    axes[0,j].plot(r,2*(1-r*r),color='#202b3a',lw=2.5,label='Poiseuille exact')
    for p,c,data in zip(profiles,colors,reference):
        profile=data['velocity_profiles'][j]; x=np.array(profile['r_over_R'])
        y=np.array(profile['axial_velocity_m_s'])/data['mean_velocity_m_s']
        exact=np.array(profile['analytic_axial_velocity_m_s'])/data['mean_velocity_m_s']
        axes[0,j].plot(x,y,color=c,lw=1.1,ls='--',label=p)
        axes[1,j].plot(x,(y-exact)*1e3,color=c,lw=1.2)
    axes[0,j].set_title(f"s/L = {reference[-1]['velocity_profiles'][j]['s_over_L']:.2f}")
    axes[1,j].set_xlabel('Signed diameter coordinate r/R')
    axes[1,j].axhline(0,color='gray',lw=.7)
    for ax in axes[:,j]: ax.grid(alpha=.18)
axes[0,0].set_ylabel('Axial velocity / Umean'); axes[1,0].set_ylabel('Error / Umean × 1000')
axes[0,0].legend(fontsize=9)
fig.suptitle('Traction-consistent analytic reference: three internal velocity sections')
fig.subplots_adjust(bottom=.15,top=.9,hspace=.2,wspace=.28)
note(fig,'The reference adds the exact end-face tangential traction. Velocity is solved, never imposed at either end.')
save(fig,'velocity_profile.png','Velocity profiles and resolved signed errors in the analytic reference',[f'pipe_{p}_reference' for p in profiles])

fig,(ax,bx)=plt.subplots(2,1,figsize=(10,8),sharex=True,gridspec_kw={'height_ratios':[2,1]})
s=np.linspace(0,1,100); dp=reference[-1]['analytic_delta_p_pa']
ax.plot(s,dp*(1-s),color='#202b3a',label='Poiseuille exact',lw=2)
for p,c,data in zip(profiles,colors,reference):
    sections=data['pressure_sections']; x=[v['s_over_L'] for v in sections]; y=np.array([v['mean_pressure_pa'] for v in sections])
    ax.plot(x,y,'o',ms=4,mfc='none',color=c,label=f'{p} reference')
    bx.plot(x,(y-dp*(1-np.array(x)))*1e3,'o-',ms=3,color=c,label=p)
n=natural[-1]['pressure_sections']
ax.plot([v['s_over_L'] for v in n],[v['mean_pressure_pa'] for v in n],'--',color='#ab5865',label='Formal BC, fine (end effect)')
ax.set_ylabel('Section mean pressure (Pa)'); ax.legend(fontsize=10)
bx.set(xlabel='Axial position s/L',ylabel='Reference error (mPa)'); bx.axhline(0,color='gray',lw=.7)
for a in (ax,bx): a.grid(alpha=.18)
fig.suptitle('Pressure gauge from traction: means over nine internal disks')
fig.subplots_adjust(bottom=.15,top=.92,hspace=.14)
note(fig,'Each disk uses 12 radial Gauss × 64 azimuthal samples. No pressure DOF is pinned.')
save(fig,'pressure_axial.png','Area-weighted pressure section means; formal end effect shown separately',[f'pipe_{p}_reference' for p in profiles]+['pipe_fine_natural'])

fig,ax=plt.subplots(figsize=(10,7))
q=np.array([.5,1,2]); vals=[read(name)['lambda_pa'] for name in ('q_half','pipe_medium_natural','q_double')]
x=np.linspace(0,2.1,100)
ax.plot(x,dp*x,'--',color='#64748b',label='Poiseuille Δp (different end traction)')
ax.plot(x,natural[1]['lambda_pa']*x,color='#177d6d',label='Formal BC: line through Q0 result')
ax.plot(q,vals,'o',color='#177d6d',ms=8,label='Three independent FEM solves')
ax.set(xlabel='Q / Q0   (Q0 = 10⁻¹⁴ m³/s)',ylabel='Inlet normal traction multiplier λ (Pa)',title='Formal boundary conditions: λ scales linearly with total flow')
ax.grid(alpha=.18); ax.legend(loc='upper left')
ax.text(.97,.12,'Finite pipe end effect at Q0:\nλ is 1.361% below Poiseuille Δp',transform=ax.transAxes,ha='right',fontsize=11)
fig.subplots_adjust(bottom=.18)
note(fig,'Medium mesh fixed; μ = 0.003 Pa·s. The Poiseuille line is a comparison, not the exact formal-BC solution.')
save(fig,'lambda_vs_Q.png','Actual Q scaling results for formal BC',['q_half','pipe_medium_natural','q_double'])

fig,ax=plt.subplots(figsize=(10,7)); cells=[6014,19541,55505]
for key,label,c,marker in [('velocity_L2_relative_error','Velocity section L2','#177d6d','o'),('pressure_profile_error','Pressure section L2','#4282b1','s'),('lambda_relative_error','Multiplier λ','#de883b','^')]:
    ax.loglog(cells,[r[key] for r in reference],marker+'-',color=c,label=label,lw=1.8,ms=7)
ax.axhline(1e-4,color='#7c5267',ls='--',label='Fine-mesh acceptance: 10⁻⁴')
ax.set_xticks(cells,[f'{p}\n{n:,} tetra' for p,n in zip(profiles,cells)])
ax.set(xlabel='Three independent quadratic-geometry pipe meshes',ylabel='Relative analytic error',title='Traction-consistent reference: measured mesh convergence')
ax.grid(which='both',alpha=.17); ax.legend(fontsize=10)
fig.subplots_adjust(bottom=.20)
note(fig,'Only canonical pipe convergence is assessed. This makes no convergence claim about Stage 1 vascular meshes.')
save(fig,'mesh_convergence.png','Three-mesh measured analytic error decrease',[f'pipe_{p}_reference' for p in profiles])

fig,ax=plt.subplots(figsize=(11,7)); names=[f'pipe_{p}_{m}' for m in ('natural','reference') for p in profiles]
for offset,key,label,c in [(-.10,'relative_inlet_constraint_error','Inlet constraint','#177d6d'),(.10,'relative_mass_closure','Global mass closure','#4282b1')]:
    values=np.array([read(name,'flux')[key] for name in names])
    # Exact floating-point zero is explicitly marked at the display floor.
    ax.semilogy(np.arange(6)+offset,np.maximum(values,1e-17),'o',color=c,label=label,ms=8)
    for j,v in enumerate(values): ax.annotate(f'{v:.1e}',(j+offset,max(v,1e-17)),xytext=(0,12 if offset<0 else -19),textcoords='offset points',ha='center',fontsize=9,color=c)
ax.axhline(1e-10,color='#ad5162',ls='--',label='Unchanged hard gate: 10⁻¹⁰')
ax.set_xticks(range(6),[f'{m}\n{p}' for m in ('formal BC','reference') for p in profiles]); ax.set_ylim(3e-18,4e-10)
ax.set(ylabel='Relative error normalized by Q0',title='Fluxes integrated directly from the FEM velocity field')
ax.grid(axis='y',alpha=.18); ax.legend(fontsize=10,loc='upper right')
fig.subplots_adjust(bottom=.20)
note(fig,'Zero is displayed at 10⁻¹⁷ and labelled 0.0e+00. Q0 = 10⁻¹⁴ m³/s; no target values are substituted for integrals.')
save(fig,'flux_balance.png','Direct FEM inlet and mass closure errors',names)

# Rendering of fetched FEM values; no analytic field enters this panel.
plot=pv.Plotter(off_screen=True,window_size=(1500,950))
plot.set_background('white')
plot.add_mesh(surface,color='#c2ccd5',opacity=.12)
cut=mesh.slice(normal=(0,1,0),origin=(25,0,0))
plot.add_mesh(cut,scalars='speed_mm_s',cmap='viridis',lighting=False,clim=(0,float(mesh['speed_mm_s'].max())),
    scalar_bar_args={'title':'Speed (mm/s)','vertical':False,'position_x':.32,'position_y':.12,'width':.48,'height':.07,'title_font_size':17,'label_font_size':15,'color':'#1f2937'})
plot.add_mesh(surface.outline(),color='#9ca3af',line_width=1)
plot.add_text('Actual FEM velocity | formal boundary conditions | medium mesh',position='upper_left',font_size=19,color='#1f2937')
plot.add_text('SYNTHETIC NUMERICAL BENCHMARK · NOT EXPERIMENTAL CONDITION\nAxial slice through the pipe; coordinates displayed in µm. No velocity profile imposed.',position='lower_left',font_size=12,color='#4b5563')
plot.show_bounds(grid=False,location='outer',xtitle='s (µm)',ytitle='y (µm)',ztitle='z (µm)',font_size=12,color='#4b5563',all_edges=False)
plot.camera_position=[(68,70,37),(25,0,0),(0,0,1)]
plot.enable_parallel_projection(); plot.camera.zoom(1.2)
plot.screenshot(str(OUT/'velocity_slice_3d.png')); plot.close()
images['velocity_slice_3d.png']={'sha256':sha256(OUT/'velocity_slice_3d.png'),'description':'Actual medium formal-BC FEM velocity on an axial 3D slice','source_cases':['pipe_medium_natural']}
write_json(OUT/'visualization_manifest.json',{'timestamp':timestamp(),'renderer':'WSL matplotlib + PyVista, actual fetched FEM fields/diagnostics','display_units_only':'coordinates m to µm; speed m/s to mm/s; solver remains SI','images':images,'source_sha256':{str(p.relative_to(ROOT)):sha256(p) for p in CASES.rglob('*') if p.is_file() and (p.name=='analytic_comparison.json' or p.name=='visualization_mesh.npz' or p.name=='flux.json')}})
print('Rendered all seven Stage 2 figures in WSL')
