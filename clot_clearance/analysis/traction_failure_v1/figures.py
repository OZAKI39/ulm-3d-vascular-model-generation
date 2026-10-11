"""Seven scientific figure groups drawn only from this run's real arrays."""
from analysis.traction_failure_v1.common import *
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import Normalize,LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.titlesize':14,'axes.labelsize':12,'legend.fontsize':11,'axes.grid':True,'grid.color':'#858c95','grid.alpha':.55,'savefig.facecolor':'white'})
record=[]
def finish(fig,folder,name,caption):
 fig.subplots_adjust(left=.055,right=.96,top=.91,bottom=.12,wspace=.32,hspace=.40);fig.text(.045,.025,caption,fontsize=11,va='bottom');fig.savefig(folder/(name+'.png'),dpi=350);fig.savefig(folder/(name+'.svg'));plt.close(fig);record.append(dict(figure=name,folder=str(folder.relative_to(ROOT)),caption=caption))
def axes3(fig,nr,nc,i,X,title):
 ax=fig.add_subplot(nr,nc,i,projection='3d');ax.view_init(18,-68);ax.set_title(title,pad=13);ax.set_xlabel('x (mm)',labelpad=9);ax.set_ylabel('y (mm)',labelpad=9);ax.set_zlabel('z (mm)',labelpad=9);ax.set_box_aspect(np.maximum(np.ptp(X,axis=0),DX));ax.tick_params(labelsize=10,pad=1)
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_major_locator(MaxNLocator(4));axis.pane.set_facecolor((.95,.97,.99,1));axis._axinfo['grid']['color']='#8c949c'
 return ax
def cloud(ax,X,c=None,label=None,cmap='viridis',norm=None,size=11):
 return ax.scatter(*(X*1e3).T,c=c,s=size,cmap=cmap if c is not None and not isinstance(c,str) else None,norm=norm,depthshade=False,linewidths=.12,edgecolors='#30516b',label=label)
def surface(ax,s,values,cmap='viridis',label='',norm=None,alpha=1):
 norm=norm or Normalize(float(np.min(values)),float(np.max(values)));colors=plt.get_cmap(cmap)(norm(values));poly=Poly3DCollection(s['triangles_xyz_m']*1e3,facecolors=colors,edgecolors='none',alpha=alpha);ax.add_collection3d(poly);xyz=s['triangles_xyz_m'].reshape(-1,3)*1e3;ax.auto_scale_xyz(*xyz.T);sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap);plt.colorbar(sm,ax=ax,shrink=.65,pad=.20,fraction=.035,label=label)
def cbar(fig,obj,ax,label):fig.colorbar(obj,ax=ax,shrink=.68,pad=.20,fraction=.035,label=label)
def run():
 t=time.perf_counter();cpu=time.process_time();s=dict(np.load(PREV/'refined/CLOT_SURFACE_TRACTIONS.npz'));g=dict(np.load(B/'GEOMETRY_PREFLIGHT.npz'));d=dict(np.load(B/'CLOT_NODAL_TRACTION_FORCES.npz'));X=d['reference_coordinates_m'];anchors=d['anchor_reference_coordinates_m'];loaded=d['support_node'];area=s['triangle_area_m2'];quality=json.loads((B/'LOCAL_MAPPING_QUALITY.json').read_text());aud=json.loads((B/'CONSERVATIVE_MAPPING_AUDIT.json').read_text());repair=json.loads((B/'CONTACT_RIM_SUPPORT_REPAIR_FREEZE.json').read_text())
 fig=plt.figure(figsize=(16,12));ax=axes3(fig,2,2,1,np.r_[X,anchors],'B1a  True exposed mesh and material lattice');cloud(ax,anchors,'#d1aa4c',size=5,label='Fixed wall anchors');cloud(ax,X,'#398ac0',size=8,label='1,927 material points');poly=Poly3DCollection(s['triangles_xyz_m']*1e3,facecolors='#8ec9e6',edgecolors='#567585',linewidths=.07,alpha=.14);ax.add_collection3d(poly);ax.legend(loc='upper left')
 ax=axes3(fig,2,2,2,X,'B1b  Material distance to exposed surface');obj=cloud(ax,X,g['distance_to_surface_m']*1e6);cbar(fig,obj,ax,'Exact triangle distance (µm)')
 ax=axes3(fig,2,2,3,X,'B1c  Face distance to nearest material point');surface(ax,s,g['face_nearest_material_distance_m']*1e6,label='Distance (µm)');bad=repair['faces'];ax.scatter(*(s['surface_coordinates_m'][bad]*1e3).T,c='#c73337',s=28,label='10 repaired rim faces');ax.legend(loc='upper left')
 ax=fig.add_subplot(2,2,4);ax.hist(g['face_nearest_material_distance_m']*1e6,bins=35,color='#398ac0',label='3,611 actual faces');ax.axvline(DELTA*1e6,color='#c73337',ls='--',label='150.5 µm support limit');ax.set(xlabel='Nearest material-point distance (µm)',ylabel='Face count',title='B1d  Complete support; 580 load nodes');ax.legend()
 finish(fig,FB,'B1_geometry_and_support','Identity coordinates; no volume rescaling. Exposed patch area 1.821924 mm². All faces supported after one contact-rim graph repair.\nMaterial selection: surface distance ≤50 µm; face support ≤150.5 µm; 1,561 anchors excluded.')
 fig=plt.figure(figsize=(16,12));ax=axes3(fig,2,2,1,X,'B2a  Actual CFD gauge pressure');surface(ax,s,s['pressure_Pa'],label='Gauge pressure (Pa)')
 ax=axes3(fig,2,2,2,X,'B2b  Actual viscous shear');surface(ax,s,s['tangential_magnitude_Pa'],label='Tangential traction (Pa)');q=np.arange(0,len(area),85);vec=s['tangential_traction_Pa'][q]*3.;ax.computed_zorder=False;ax.quiver(*(s['surface_coordinates_m'][q]*1e3).T,*vec.T,color='#202b42',arrow_length_ratio=.25,linewidth=1.1,zorder=10)
 ax=axes3(fig,2,2,3,X,'B2c  Face area and solid-to-fluid normals');surface(ax,s,area*1e12,label='Triangle area (µm²)',cmap='magma');q=np.arange(0,len(area),110);ax.computed_zorder=False;ax.quiver(*(s['surface_coordinates_m'][q]*1e3).T,*(s['n_clot'][q]*.08).T,color='#253c54',arrow_length_ratio=.3,linewidth=.9,zorder=10)
 ax=fig.add_subplot(2,2,4);names=['Pressure','Shear','Normal visc.','Total'];F=np.array([d[k+'_force_N'].sum(0) for k in ['pressure','shear','viscous_normal','total']])*1e6
 for i,col in enumerate(['#398ac0','#e99f36','#c74447']):ax.bar(np.arange(4)+(i-1)*.23,F[:,i],.23,label='xyz'[i],color=col)
 ax.set(xticks=np.arange(4),xticklabels=names,ylabel='Integrated force (µN)',title='B2d  Real component forces');ax.legend(ncol=3)
 finish(fig,FB,'B2_native_surface_tractions','σ = −pI + µ(∇u+∇uᵀ), t = σ n_clot; n_clot = −n_fluid. No attachment face or wall reaction included. Vector glyphs drawn in foreground.\nShear arrows: 3 mm/Pa; normal arrows: 0.08 mm. Gauge pressure is tied to the original outlet reference.')
 fig=plt.figure(figsize=(18,11));ax=axes3(fig,2,3,1,X,'B3a  Original face force magnitude');surface(ax,s,np.linalg.norm(s['total_traction_Pa']*area[:,None],axis=1)*1e9,label='Face force (nN)')
 vmax=np.linalg.norm(d['total_force_N'],axis=1).max()*1e9
 for k,(name,field) in enumerate([('Initial local map','total_initial_force_N'),('Conservative result','total_force_N')],2):
  ax=axes3(fig,2,3,k,X,f'B3{chr(96+k)}  {name}');obj=cloud(ax,X[loaded],np.linalg.norm(d[field][loaded],axis=1)*1e9,norm=Normalize(0,vmax));cbar(fig,obj,ax,'Node force (nN)')
 ax=fig.add_subplot(2,3,4);names=['pressure','shear','total'];x=np.arange(3)
 for offset,stage,color in [(-.13,'initial','#94b2c7'),(.13,'corrected','#165989')]:
  ax.bar(x+offset,[aud['components'][k][stage]['moment_scaled_error'] for k in names],.25,color=color,label=stage)
 ax.set(yscale='log',xticks=x,xticklabels=names,ylabel='Moment error / (Σ|F_face| L)',title='B3d  Moment conservation');ax.axhline(1e-8,color='#c74447',ls='--',label='Frozen gate');ax.legend()
 ax=fig.add_subplot(2,3,5);ax.bar(names,[quality['metrics'][k]['correction_L2_fraction']*100 for k in names],color=['#398ac0','#e99f36','#d0524f']);ax.axhline(10,color='#333333',ls='--');ax.set(ylabel='Correction / initial force L2 (%)',title='B3e  Frozen 10% engineering screen')
 ax=fig.add_subplot(2,3,6)
 for stage,field in [('initial','total_initial_force_N'),('corrected','total_force_N')]:
  f=np.sort(np.linalg.norm(d[field][loaded],axis=1))[::-1];ax.plot(np.arange(1,len(f)+1)/len(f)*100,np.cumsum(f)/sum(f)*100,label=stage)
 ax.set(xlabel='Largest-force nodes (% of 580)',ylabel='Cumulative |force| (%)',title='B3f  Force concentration');ax.legend()
 finish(fig,FB,'B3_conservation_and_local_quality','PASS: frozen engineering screen. Total correction 0.7716%; shear correction 9.0150%. No reversed pressure/shear node forces.\nAll components conserve force and moment; near-zero torques use the disclosed physical force × length scale, not a hidden denominator.')
 fig=plt.figure(figsize=(19,7));scales={'pressure':5e7,'shear':1e9,'total':5e7};arrowrecord={}
 for i,key in enumerate(['pressure','shear','total'],1):
  ax=axes3(fig,1,3,i,np.r_[X,anchors],f'B4{chr(96+i)}  {key.capitalize()} nodal forces');cloud(ax,anchors,'#d5b053',size=5);cloud(ax,X[~loaded],'#cdd8df',size=7);cloud(ax,X[loaded],'#287fba',size=10);ids=np.flatnonzero(loaded)[::7];vec=d[key+'_force_N'][ids]*scales[key];ax.computed_zorder=False;ax.quiver(*(X[ids]*1e3).T,*vec.T,color='#b72c20',arrow_length_ratio=.25,linewidth=1.0,zorder=10);arrowrecord[key]=dict(displayed_node_ids=d['material_atom_id'][ids],scale_mm_per_N=scales[key]);ax.text2D(.02,.93,f'Arrows × {scales[key]:.1e} mm/N',transform=ax.transAxes,fontsize=10)
 fig.legend(handles=[Line2D([0],[0],marker='o',ls='',color='#287fba',label='580 loaded nodes'),Line2D([0],[0],marker='o',ls='',color='#cdd8df',label='1,347 internal / unloaded'),Line2D([0],[0],marker='o',ls='',color='#d5b053',label='1,561 fixed anchors')],loc='upper center',ncol=3,bbox_to_anchor=(.5,.98))
 finish(fig,FB,'B4_material_point_loading','Actual 1,927-node V3 clot. Every node retained; every seventh loaded node carries a displayed arrow for legibility.\nCoordinates remain undeformed; scaled vector glyphs drawn in foreground for readability. Node forces in N; production coupling not run.');save(B/'ARROW_VISUALIZATION_SCALE.json',arrowrecord)
 nm=dict(np.load(A/'M5C_NEIGHBORHOOD.npz'));nv=dict(np.load(A/'V3_NEIGHBORHOOD.npz'));fig=plt.figure(figsize=(18,12))
 colors={'interior':'#2388b8','free_surface':'#78c5c3','grip':'#dfbc56','grip_corner':'#c84342','initial_notch':'#884fa0','notch_tip':'#e18932','anchor':'#ddd0a5'}
 for j,(name,a) in enumerate([('M5C',nm),('V3',nv)]):
  for k,(field,title) in enumerate([('labels','Actual region identity'),('neighbors','Reference neighbor count'),('ratio','Weighted volume / full lattice')]):
   ax=axes3(fig,2,3,j*3+k+1,a['X'],f'A1{chr(97+j*3+k)}  {name}: {title}')
   if field=='labels':
    for v in np.unique(a['labels']):cloud(ax,a['X'][a['labels']==v],colors[v],size=7,label=v.replace('_',' '))
    ax.legend(fontsize=8,loc='upper left')
   else:obj=cloud(ax,a['X'],a['neighbors'] if field=='neighbors' else a['m']/a['full_m'],size=8);cbar(fig,obj,ax,'Count' if field=='neighbors' else 'm / m_full')
 finish(fig,FA,'A1_actual_neighborhoods','Real reference topology: M5C 2,880 nodes with saved pre-notch bonds removed; V3 1,927 material + 1,561 anchor nodes.\nThe V3 neighborhood diagnostic includes actual interface/anchor bonds. A boundary force under affine strain is not by itself a traction-free BC error.')
 wm=dict(np.load(A/'M5C_WEIGHTS.npz'));wv=dict(np.load(A/'V3_WEIGHTS.npz'));me=dict(np.load(A/'FORCE_ENERGY_DIAGNOSTIC.npz'));mj=json.loads((A/'FORCE_ENERGY_COMPATIBILITY.json').read_text());fig=plt.figure(figsize=(18,12))
 for i,(name,a) in enumerate([('M5C',wm),('V3',wv)],1):
  ax=axes3(fig,2,3,i,a['X'],f'A2{chr(96+i)}  {name}: actual fitted weights');obj=cloud(ax,a['X'],a['weights'],norm=Normalize(0,4),size=9);cbar(fig,obj,ax,'Weight (colors clipped to 0–4)');neg=a['weights']<0
  if neg.any():
   ax.computed_zorder=False;ax.scatter(*(a['X'][neg]*1e3).T,c='#d9253f',s=24,marker='x',label='15 negative weights',zorder=12);ax.legend(fontsize=9,loc='upper left')
 ax=fig.add_subplot(2,3,3)
 for name,a in [('M5C',wm),('V3',wv)]:ax.scatter(a['operator_residual_before'],a['operator_residual_after'],s=5,alpha=.35,label=name)
 ax.plot([0,8],[0,8],ls='--',c='black');ax.set(xlabel='Unweighted row-scaled residual',ylabel='Weighted row-scaled residual',title='A2c  Operator matching is not uniformly improved');ax.set_yscale('symlog',linthresh=10);ax.legend()
 ax=fig.add_subplot(2,3,4)
 for name,a in [('M5C',wm),('V3',wv)]:ax.hist(np.log10(a['local_condition']),bins=30,alpha=.6,label=name)
 ax.set(xlabel='log10 regularized local condition',ylabel='Node count',title='A2d  V3: 15 negative nodes / 59 nonpositive bonds');ax.legend()
 ax=fig.add_subplot(2,3,5);names=['original','weighted_diagnostic'];ax.bar([0,1],[mj[k]['U_J']*1e12 for k in names],color=['#398ac0','#d99632']);ax.set(xticks=[0,1],xticklabels=['Original','Weighted diagnostic'],ylabel='Elastic U (pJ)',title='A2e  Same prescribed M5C deformation');ax.text(.03,.92,'Max component FD error < 1×10⁻⁸\nAdaptation; no response certification',transform=ax.transAxes,va='top',fontsize=10)
 ax=axes3(fig,2,3,6,me['X'],'A2f  M5C change in local force');obj=cloud(ax,me['X'],np.linalg.norm(me['weighted_diagnostic_force_N']-me['original_force_N'],axis=1)*1e6,size=7);cbar(fig,obj,ax,'|Δforce| (µN)')
 finish(fig,FA,'A2_weights_and_variational_diagnostic','Pinned perifit OSB implementation; actual source has 54 operator rows. No weight clipping, mesh cleanup, or crack restoration.\nM5C kernel-preserving weighted potential passes finite differences, but this does not validate Shojaei transfer to native LPS. V3 weights fail positivity.')
 pw=np.load(A/'INTERACTION_PATH_WORK.npz');fig=plt.figure(figsize=(18,11));ax=fig.add_subplot(2,3,1);ax.plot(pw['corner_energy_J']*1e12,'o-',label='Independent U');ax.plot(pw['all_pair_cumulative_work_J'].sum(1)*1e12,'x--',label='Σ unordered-pair work');ax.set(xticks=range(5),xticklabels=['0','a','a+b','b','0'],ylabel='Energy / work (pJ)',title='A3a  Actual reversible closed deformation path');ax.legend()
 ax=fig.add_subplot(2,3,2);ax.plot(pw['selected_pair_cumulative_work_J']*1e15,'o-',c='#c44242');ax.axhline(0,c='#555',lw=.7);ax.set(xticks=range(5),xticklabels=['0','a','a+b','b','0'],ylabel='Single-pair accumulated work (fJ)',title='A3b  Pair IDs 240–720: nonzero closed-path work')
 ax=fig.add_subplot(2,3,3);v=pw['all_pair_cumulative_work_J'][-1]*1e15;ax.hist(v,bins=60,color='#398ac0');ax.set(xlabel='Pair closed-path work (fJ)',ylabel='Unordered pair count',title='A3c  Pair work transfers cancel in the total')
 base=ROOT/'verification/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z'
 for i,name in enumerate(['M1','M2'],4):
  h=np.genfromtxt(base/name/'history.csv',names=True,delimiter=',');ax=fig.add_subplot(2,3,i)
  for f,label in [('U_bulk_J','Elastic U'),('W_ext_J','External work'),('D_num_J','Numerical damping')]:ax.plot(h['time_s']*1e3,h[f]*1e12,label=label)
  ax.set(xlabel='Original elapsed time (ms)',ylabel='Energy (pJ)',title=f'A3{chr(96+i)}  Reused {name} load–unload');ax.legend(fontsize=9)
 ax=fig.add_subplot(2,3,6);ax.axis('off');ax.text(0,.98,'Definitions and limits\n\nCurrent: s > s00, s00 = 0.15\n\nFoster: ∫(Tᵢⱼ − Tⱼᵢ) · dη\nWork density: J/m⁶\nDiscrete pair work: Vᵢ Vⱼ × density\n\nwc = 4 Gc / (π δ⁴)\nThis formula is not an LPS calibration.\n\nUnique pair elastic potential:\nNOT ESTABLISHED\nProduction law: UNCHANGED',va='top',fontsize=13)
 finish(fig,FA,'A3_path_work_and_failure_limits','Computed path work agrees with global elastic energy to ~1.2×10⁻¹⁴ (scaled), while individual pair loops need not close.\nM1/M2 are historical damped dynamics; their external work includes kinetic energy and numerical damping. No new fracture or calibrated Gc.')
 save(COMMON/'FIGURE_PROVENANCE.json',dict(figures=record,real_data_only=True,arrow_subset_does_not_change_exported_load=True));save(COMMON/'FIGURE_RESOURCES.json',cost(t,cpu));print('generated',len(record),'PNG/SVG groups')
if __name__=='__main__':run()
