"""Reduce accepted native trajectories; unknown events are never zero thresholds."""
from common import *
from prescribed_traction import Load,MODES,END,RISE
import csv,html
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection,Poly3DCollection
from matplotlib.ticker import MaxNLocator

def read(p):return json.loads(Path(p).read_text())
def writecsv(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def collect():
 cases=[];rows=[]
 for p in sorted((O/'cases').glob('*/SUMMARY.json')):
  z=read(p);z['path']=p.parent;z['tag']=p.parent.name;cases.append(z)
 cases.sort(key=lambda z:(MODES.index(z['params']['mode']),z['params']['peak_pa']))
 for z in cases:
  q=z['last_row'];ev=z['events'];cs=z['final_free_components'];largest=max(cs,key=lambda c:c['nodes']) if cs else None
  rows.append(dict(specimen_id='precracked_1792_clot_336_bottom_anchors',load_mode=z['params']['mode'],peak_traction_pa=z['params']['peak_pa'],peak_total_force_n=z['loaded_area_m2']*z['params']['peak_pa'],observed_peak_total_force_n=q['force_norm_n'],application_area_m2=z['loaded_area_m2'],duration_s=END,accepted_duration_s=z['time_s'],impulse_ns=q['impulse_ns'],external_work_j=q['W_J'],first_softening_s=ev.get('first_softening',{}).get('time_s'),first_failure_s=ev.get('first_failure',{}).get('time_s'),first_multinode_s=ev.get('first_multinode',{}).get('time_s'),max_detached_nodes=z['history_max_free_nodes'],final_detached_nodes=q['largest_detached_nodes'],final_total_detached_nodes=q['detached_nodes'],unresolved_debris_mass=q['debris_1_2_mass_kg'],internal_bond_survival=largest['internal_alive_fraction'] if largest else None,energy_residual=q['R_J'],max_relative_energy_residual=z['max_relative_energy_residual'],classification=z['classification'],numerical_status=z['status'],physical_status='MECHANICAL_REQUIREMENT_ONLY',evidence_path=str(z['path'].relative_to(R))))
 for z,row in zip(cases,rows):
  hist=list(csv.DictReader((z['path']/'history.csv').open()))
  for key in ['N_ge14','N_ge34','N_ge114','N_ge200um','N_ge300um']:
   hits=[float(q['time_s']) for q in hist if int(q[key])>0];row['first_'+key+'_s']=hits[0] if hits else None;row['history_max_'+key]=max(int(q[key]) for q in hist);row['final_'+key]=int(hist[-1][key])
  dh=read(z['path']/'DETACHED_COMPONENT_HISTORY.json');peaks=[(q,c) for q in dh for c in q['components'] if c['nodes']==z['history_max_free_nodes'] and c['nodes']>0]
  row['history_largest_internal_alive_fraction']=peaks[0][1]['internal_alive_fraction'] if peaks else None
  row['history_max_speed_m_s']=max(float(q['max_speed_m_s']) for q in hist);row['history_max_damage_increment']=max(float(q['max_dd']) for q in hist);row['accepted_sub_base_dt_steps']=sum(float(q['dt_s'])<BASE*(1-1e-10) for q in hist)
  row['history_largest_first_time_s']=peaks[0][0]['time_s'] if peaks else None
  row['first_softening_instantaneous_pa']=z['params']['peak_pa']*float(__import__('prescribed_traction').envelope(row['first_softening_s'])) if row['first_softening_s'] is not None else None
  row['first_failure_instantaneous_pa']=z['params']['peak_pa']*float(__import__('prescribed_traction').envelope(row['first_failure_s'])) if row['first_failure_s'] is not None else None
 writecsv(O/'MECHANICAL_LOAD_REQUIREMENT_TABLE.csv',rows)
 brackets=[]
 for mode in MODES:
  cc=[z for z in cases if z['params']['mode']==mode];valid=[z for z in cc if z['status']=='COMPLETED']
  result={'load_mode':mode,'no_damage_observed_pa':json.dumps([z['params']['peak_pa'] for z in valid if not z['events']]),'softening_observed_pa':json.dumps([z['params']['peak_pa'] for z in cc if 'first_softening' in z['events']]),'failure_observed_pa':json.dumps([z['params']['peak_pa'] for z in cc if 'first_failure' in z['events']]),'multinode_observed_pa':json.dumps([z['params']['peak_pa'] for z in cc if 'first_multinode' in z['events']]),'secondary_observed_pa':json.dumps([z['params']['peak_pa'] for z in cc if z['secondary_fracture_events']]),'unresolved_pa':json.dumps([z['params']['peak_pa'] for z in cc if z['status']!='COMPLETED'])}
  for label in ['no_damage_observed_pa','softening_observed_pa','failure_observed_pa','multinode_observed_pa','secondary_observed_pa']:result[label]=json.dumps(sorted(set(json.loads(result[label]))))
  result['failed_attempts_pa']=result['unresolved_pa'];valid_peaks={z['params']['peak_pa'] for z in valid};result['unresolved_pa']=json.dumps(sorted({z['params']['peak_pa'] for z in cc if z['status']!='COMPLETED' and z['params']['peak_pa'] not in valid_peaks}))
  for label,key in [('softening','first_softening'),('failure','first_failure'),('multinode','first_multinode')]:
   observed=[z['params']['peak_pa'] for z in cc if key in z['events']];upper=min(observed) if observed else None;neg=[z['params']['peak_pa'] for z in valid if key not in z['events'] and (upper is None or z['params']['peak_pa']<upper)]
   result[label+'_negative_lower_pa']=max(neg) if neg else None;result[label+'_positive_upper_pa']=upper
  result['interpretation']='Observed finite-window cases; no-event side requires completed window; positive accepted-prefix event survives later numerical failure. Not universal critical stress.';brackets.append(result)
 writecsv(O/'LOAD_BRACKETING_SUMMARY.csv',brackets)
 dump(O/'FRACTURE_EVENTS.json',{z['tag']:z['events'] for z in cases})
 dump(O/'DETACHED_COMPONENT_HISTORY.json',{z['tag']:dict(evidence=str((z['path']/'DETACHED_COMPONENT_HISTORY.json').relative_to(R)),sha256=sha(z['path']/'DETACHED_COMPONENT_HISTORY.json'),history_max_free_nodes=z['history_max_free_nodes'],final_free_components=z['final_free_components'],secondary_fracture_events=z['secondary_fracture_events']) for z in cases})
 dump(O/'ENERGY_LEDGER.json',dict(definition='W = delta U + delta K + D + delta Uc - Jc + R; trapezoidal external work on each actual displacement; R is not dissipation',cases={z['tag']:dict(final_accepted={k:z['last_row'][k] for k in ['W_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J','relative_R','impulse_ns','max_speed_m_s']},max_relative_residual=z['max_relative_energy_residual'],history=str((z['path']/'history.csv').relative_to(R))) for z in cases}))
 dump(O/'RUN_RESOURCES.json',dict(native_wall_s=sum(z['wall_s'] for z in cases),native_cpu_s=sum(z['cpu_s'] for z in cases),peak_RSS_bytes=max(z['peak_RSS_bytes'] for z in cases),accepted_steps=sum(z['accepted_steps'] for z in cases),event_substeps=sum(z['predictive_substeps'] for z in cases),predictive_rejections=sum(z['predictive_rejections']-(read(Path(z['params']['replay_source'])/'SUMMARY.json')['predictive_rejections'] if z['params'].get('replay_source') else 0) for z in cases),native_replayed_steps=sum(z.get('prefix_replay_steps',0) for z in cases),normal_modes_max_concurrent_cases=2,replay_may_overlap=True,peak_RSS_scope='Maximum per native process; independent normal modes and necessary replay may overlap',ranks=1,threads=1,new_full_V3_runs=0,new_CFD_runs=0,new_dependencies=[],cases={z['tag']:{k:z[k] for k in ['wall_s','cpu_s','peak_RSS_bytes','accepted_steps','predictive_substeps','predictive_rejections','status']} for z in cases}))
 resource=read(O/'RUN_RESOURCES.json');preflights=[read(p) for p in (O/'cases').glob('*/FAILED_PREFLIGHT.json')];resource['failed_preflights']=preflights;resource['additional_failed_preflight_native_wall_s']=sum(q['native_total_wall_s'] for q in preflights);resource['native_wall_s_including_failed_preflights']=resource['native_wall_s']+resource['additional_failed_preflight_native_wall_s'];resource['completed_unique_physical_conditions']=len({(z['params']['mode'],z['params']['peak_pa']) for z in cases if z['status']=='COMPLETED'});resource['unique_completed_accepted_steps']=sum(z['accepted_steps'] for z in cases if z['status']=='COMPLETED');resource['cpu_total_scope']='Measured case CPU sum; failed timer-preflight CPU is unavailable and not imputed';resource['total_native_replayed_steps_including_failed_preflight']=resource['native_replayed_steps']+sum(q['replayed_native_steps'] for q in preflights);dump(O/'RUN_RESOURCES.json',resource)
 deformation={}
 for z in cases:
  st=load_npz(z['path']/'keys'/('final.npz' if z['status']=='COMPLETED' else 'last_accepted.npz'));ref=load_npz(O/'provenance/reference.npz');load=Load(ref,z['params']);q=load.quads[load.mask];xx=st['x'][q]
  av=.5*(np.cross(xx[:,1]-xx[:,0],xx[:,2]-xx[:,0])+np.cross(xx[:,2]-xx[:,0],xx[:,3]-xx[:,0]));area=np.linalg.norm(av,axis=1);normal=av/np.maximum(area[:,None],1e-30);direction=load.audit(ref)['direction'];normal_fraction=np.abs(normal@direction)
  deformation[z['tag']]=dict(scope='Diagnostic original quad geometry only; not a reconstructed fluid surface after fracture',area_ratio_min=float(np.min(area/load.face_area[load.mask])),area_ratio_max=float(np.max(area/load.face_area[load.mask])),normal_rotation_max_deg=float(np.degrees(np.arccos(np.clip(normal[:,2],-1,1))).max()),applied_direction_current_normal_fraction_min=float(normal_fraction.min()),applied_direction_current_normal_fraction_max=float(normal_fraction.max()),interpretation='Frozen reference nominal directions; large rotation changes the current shear/normal decomposition. Never claimed as follower pressure or current-surface pure shear.')
 load_partition={};g=module('load_partition_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
 for z in cases:
  st=load_npz(z['path']/'keys'/('final.npz' if z['status']=='COMPLETED' else 'last_accepted.npz'));ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);lab,sz,att=g.topology(ref,st['d']);free=~att[lab];force=st['external_force'];total=float(np.linalg.norm(force.sum(0)));largest=max(z['final_free_components'],key=lambda c:c['nodes']) if z['final_free_components'] else None;ids=np.array(largest['atom_ids'],int)-1 if largest else np.array([],int)
  load_partition[z['tag']]=dict(accepted_time_s=float(st['time']),total_force_n=total,force_on_detached_n=float(np.linalg.norm(force[:n][free].sum(0))),force_fraction_on_detached=float(np.linalg.norm(force[:n][free].sum(0))/max(total,1e-30)),largest_component_force_n=float(np.linalg.norm(force[ids].sum(0))),kinetic_energy_detached_j=float(.5*RHO*V*np.sum(st['velocity'][:n][free]**2)),total_kinetic_energy_j=float(st['K']),interpretation='Actual final accepted native forces/velocities; nominal material load persists on detached original nodes')
 dump(O/'LOAD_PARTITION_AUDIT.json',load_partition)
 secondary_audit={}
 for z in cases:
  births={}
  for q in read(z['path']/'DETACHED_COMPONENT_HISTORY.json'):
   for c in q['components']:births.setdefault(tuple(c['atom_ids']),q['time_s'])
  secondary_audit[z['tag']]=[dict(time_s=e['time_s'],parent_nodes=len(e['parent_atom_ids']),children_nodes=[len(c) for c in e['children_atom_ids']],parent_first_recorded_s=births.get(tuple(e['parent_atom_ids'])),lifetime_s=e['time_s']-births[tuple(e['parent_atom_ids'])] if tuple(e['parent_atom_ids']) in births else None,original_parent_atom_ids=e['parent_atom_ids']) for e in z['secondary_fracture_events']]
 chunk_audit={}
 for z in cases:
  records=read(z['path']/'DETACHED_COMPONENT_HISTORY.json');peak=z['history_max_free_nodes'];first=next(((q,c) for q in records for c in q['components'] if c['nodes']==peak and peak>0),None)
  if first is None:chunk_audit[z['tag']]=dict(observed=False);continue
  q,c=first;ids=set(c['atom_ids']);same=next((f for f in z['final_free_components'] if set(f['atom_ids'])==ids),None)
  chunk_audit[z['tag']]=dict(observed=True,history_largest_nodes=peak,first_time_s=q['time_s'],atom_ids=c['atom_ids'],initial_internal_bonds=c['internal_initial_alive_bonds'],alive_at_formation=c['internal_alive_bonds'],survival_at_formation=c['internal_alive_fraction'],same_original_members_connected_at_final=same is not None,final_internal_alive_fraction=same['internal_alive_fraction'] if same else None,observed_persistence_s=z['time_s']-q['time_s'] if same else None,final_descendant_sizes=[len(ids & set(f['atom_ids'])) for f in z['final_free_components'] if ids & set(f['atom_ids'])],largest_secondary_parent_nodes=max((len(e['parent_atom_ids']) for e in z['secondary_fracture_events']),default=0))
 dump(O/'COHESIVE_CHUNK_AUDIT.json',chunk_audit)
 dump(O/'SECONDARY_FRAGMENTATION_AUDIT.json',secondary_audit)
 dump(O/'LOAD_DEFORMATION_AUDIT.json',deformation)
 return cases,rows,brackets

def figures(cases,geometry_only=False):
 # Superseded failed attempts remain in tables/resources, not the physical response curve.
 valid_keys={(z['params']['mode'],z['params']['peak_pa']) for z in cases if z['status']=='COMPLETED'}
 cases=[z for z in cases if z['status']=='COMPLETED' or (z['params']['mode'],z['params']['peak_pa']) not in valid_keys]
 F.mkdir(parents=True,exist_ok=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':15,'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':11,'savefig.dpi':220,'svg.fonttype':'none','grid.color':'#838b94','grid.alpha':.55})
 gallery=[];ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);b=ref['bonds'];bt=ref['bond_types'];cut=ref['initial_cut'];X=ref['X']*1e6
 def save(fig,name,caption):
  fig.text(.99,.003,'MECHANICAL_REQUIREMENT_ONLY',ha='right',va='bottom',fontsize=9,color='#5a646c')
  for ext in ['png','svg']:fig.savefig(F/(name+'.'+ext),bbox_inches='tight')
  plt.close(fig);gallery.append(dict(name=name,caption=caption))
 def axis3(ax,limits):
  for i,axis in enumerate([ax.xaxis,ax.yaxis,ax.zaxis]):axis.set_major_locator(MaxNLocator(3));axis.labelpad=10;axis.pane.set_facecolor((.96,.97,.98,1));axis._axinfo['grid']['color']=(.42,.46,.50,.65)
  ax.set(xlabel='X (μm)',ylabel='Y (μm)',zlabel='Z (μm)',xlim=limits[:,0],ylim=limits[:,1],zlim=limits[:,2]);ax.set_box_aspect(np.ptp(limits,axis=0));ax.view_init(16,-63)
 lim=np.array([X.min(0)-80,X.max(0)+150]);fig=plt.figure(figsize=(21,7))
 for k,mode in enumerate(MODES):
  ax=fig.add_subplot(1,3,k+1,projection='3d');ax.scatter(*X[:n].T,s=7,c='#9b3150',alpha=.65,label='Material (original IDs)');ax.scatter(*X[n:].T,s=6,c='#626e7a',alpha=.65,label='Fixed anchors')
  load=Load(ref,dict(mode=mode,peak_pa=100));ax.add_collection3d(Poly3DCollection(X[load.quads[load.mask]],facecolor='#39a6b7',edgecolor='#327581',linewidth=.25,alpha=.5))
  mids=np.abs(X[b].mean(1)[:,1])<30;idx=np.flatnonzero(cut&mids)[::7];ax.add_collection3d(Line3DCollection(X[b[idx]],colors='#de8c28',linewidth=1.1,alpha=.8))
  pos=load.positions[load.mask][::7]*1e6;dire=np.array([1,0,0]) if k==0 else np.array([0,0,-1 if k==1 else 1]);vec=np.tile(dire*120,(len(pos),1));ax.quiver(*pos.T,*vec.T,color='#174870',arrow_length_ratio=.28,linewidth=1.2)
  ax.quiver(325,175,375,0,0,90,color='#222',arrow_length_ratio=.3);ax.text(325,175,480,'n₀',fontsize=12);axis3(ax,lim);ax.set_title(mode.replace('_',' ')+'\n100 Pa reference arrows')
  if k==1:ax.legend(loc='upper left',bbox_to_anchor=(.02,1.02))
 fig.subplots_adjust(bottom=.16,top=.84,wspace=.13);fig.text(.5,.015,'Arrow scale: 100 Pa → 120 μm   |   Orange: original pre-cut bonds   |   Cyan: loaded reference surface',ha='center',fontsize=13)
 fig.suptitle('Prescribed nominal traction · real 1792 + 336 nodes · MECHANICAL REQUIREMENT ONLY',y=.97);save(fig,'Fig01_load_geometry','Reference-area dead traction. 91 top quads for shear; 42 left top quads for local inward/outward loading. No extra cuts. Directions stay frozen; not follower pressure.')
 if geometry_only:return gallery
 fig,axes=plt.subplots(2,2,figsize=(13,9));colors=['#246b93','#be7040','#963c65'];fields=['max_tensile_kappa','new_softened','failed','history_max_free_nodes'];titles=['Maximum tensile κ','Softened bonds (including failed)','New fully failed bonds','History maximum free component (nodes)']
 for ax,field,title in zip(axes.flat,fields,titles):
  for mode,color in zip(MODES,colors):
   cc=[z for z in cases if z['params']['mode']==mode];xx=[z['params']['peak_pa'] for z in cc];yy=[z['history_max_free_nodes'] if field=='history_max_free_nodes' else z['last_row']['bulk_failed']+z['last_row']['interface_failed'] if field=='failed' else z['last_row'][field] for z in cc]
   ax.plot(xx,yy,'o-',color=color,label=mode.replace('_',' '),lw=1.5)
   for z,x,y in zip(cc,xx,yy):
    if z['status']!='COMPLETED':ax.scatter(x,y,s=150,marker='x',c='black',zorder=9)
  ax.set_xscale('log');ax.set(xlabel='Peak nominal traction (Pa)',ylabel=title);ax.grid()
 axes[0,0].set_yscale('log')
 for ax in [axes[0,1],axes[1,0],axes[1,1]]:ax.set_yscale('symlog',linthresh=1);ax.set_ylim(bottom=-.05)
 axes[0,0].axhline(.15,c='#777',ls=':',label='Bulk s₀');axes[0,0].legend(fontsize=9);fig.suptitle('Accepted native response · × marks an incomplete numerical window');fig.tight_layout();save(fig,'Fig02_load_response','Pre-cut bonds excluded. Maximum κ includes already failed bonds and is not an intact-material strain measure. Count axes use symlog to retain actual zero responses. Lines join sampled cases, not a fitted threshold. Incomplete windows are explicit; a zero damage prefix cannot establish a full-window negative result.')
 completed=[z for z in cases if z['status']=='COMPLETED'];chosen=max(completed or cases,key=lambda z:(z['history_max_free_nodes'],z['last_row']['bulk_failed']+z['last_row']['interface_failed']))
 names=['initial','first_softening','first_failure','first_multinode','max_free','final' if chosen['status']=='COMPLETED' else 'last_accepted'];states=[load_npz(chosen['path']/'keys'/(key+'.npz')) if (chosen['path']/'keys'/(key+'.npz')).exists() else None for key in names];coords=np.concatenate([st['x']*1e6 for st in states if st is not None]);graph=module('fig_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py');core=[]
 for st in states:
  if st is None:continue
  ll,ss,aa=graph.topology(ref,st['d']);keep=np.r_[ss[ll]>=n/4,np.ones(len(ref['X'])-n,bool)];core.append(st['x'][keep]*1e6)
 core=np.concatenate(core);limits=np.array([core.min(0)-75,core.max(0)+75]);fig=plt.figure(figsize=(20,14));fig.subplots_adjust(hspace=.35,wspace=.18,top=.9,bottom=.07)
 for i,(key,s) in enumerate(zip(names,states)):
  ax=fig.add_subplot(2,3,i+1,projection='3d');axis3(ax,limits);ax.set_yticks([-100,100]);ax.tick_params(axis='y',pad=6)
  if s is None:ax.text2D(.08,.5,'NOT OBSERVED\nwithin the window',transform=ax.transAxes);ax.set_title(key.replace('_',' '));continue
  xx=s['x']*1e6;lab,sz,att=graph.topology(ref,s['d']);cc=np.full(n,'#a03c50');softnodes=np.unique(b[(s['d']>0)&(s['d']<1)&~cut]);softnodes=softnodes[softnodes<n];cc[softnodes]='#d48c24';cc[~att[lab]]='#237dad';inside=np.all((xx>=limits[0])&(xx<=limits[1]),axis=1);ax.scatter(*xx[:n][inside[:n]].T,s=10,c=cc[inside[:n]],alpha=.7);ax.scatter(*xx[n:].T,s=6,c='#68737e',alpha=.5)
  if key in chosen['events'] and 'original_bond_ids' in chosen['events'][key]:
   eb=b[np.array(chosen['events'][key]['original_bond_ids'])-1];ax.add_collection3d(Line3DCollection(xx[eb],colors='#192b3b',linewidth=2.2,linestyle='--' if key=='first_failure' else '-'));ei=np.unique(eb);ax.scatter(*xx[ei].T,s=32,c='#192b3b',depthshade=False)
  middle=np.abs(ref['X'][:,1])<=DX*.51;draw=(bt==1)&middle[b[:,0]]&middle[b[:,1]]&(s['d']<1)&~cut&inside[b[:,0]]&inside[b[:,1]];idx=np.flatnonzero(draw)[::8];ax.add_collection3d(Line3DCollection(xx[b[idx]],colors='#547278',linewidth=.4,alpha=.45))
  ax.set_title(key.replace('_',' ')+'\n'+format(float(s['time'])*1e3,'.4f')+' ms\n'+str(int(np.sum(~inside[:n])))+' material nodes outside fixed ROI')
 fig.suptitle(chosen['params']['mode'].replace('_',' ').title()+f" · peak {chosen['params']['peak_pa']:.3g} Pa · actual displacement ×1 · fixed physical ROI / camera",y=.985);fig.text(.5,.02,'Red: attached   |   Gold: softening endpoints   |   Blue: detached   |   Gray: anchors   |   Black: first-event bonds',ha='center',fontsize=13);save(fig,'Fig03_actual_fracture_sequence','Actual accepted states at displacement ×1 and a common physical ROI/camera. ROI is the joint envelope of components with at least 448 nodes plus anchors, across all shown frames. Outside-ROI material counts are explicit; the full final extent appears in Fig03b. Gold marks nodes incident to softening bonds; blue marks detached nodes. Sparse actual live bonds and first-event bonds only; no invented continuous crack surface.')
 st=states[-1];xx=st['x']*1e3;ll,ss,aa=graph.topology(ref,st['d']);fig,ax=plt.subplots(figsize=(18,4));ax.scatter(xx[:n,0],xx[:n,2],c=np.where(aa[ll],'#a03c50','#237dad'),s=12,alpha=.7,label='All final material nodes');ax.scatter(xx[n:,0],xx[n:,2],c='#68737e',s=10,label='Anchors')
 from matplotlib.patches import Rectangle
 ax.add_patch(Rectangle((limits[0,0]/1000,limits[0,2]/1000),(limits[1,0]-limits[0,0])/1000,(limits[1,2]-limits[0,2])/1000,fill=False,ec='#38836c',ls='--',lw=1.5,label='Fixed 3D ROI in Fig03'));ax.set(xlabel='X (mm)',ylabel='Z (mm)',title='All final particles · X–Z projection · no displacement magnification');ax.grid();ax.set_aspect('equal',adjustable='datalim');ax.legend(ncol=3,loc='upper center',bbox_to_anchor=(.5,-.23));save(fig,'Fig03b_full_extent_projection','All actual final material and anchor coordinates in X–Z projection. Fast detached debris is not hidden by the fixed main-body ROI used for the six 3D stages. The rectangular outline is a camera ROI, not a crack or material surface.')
 fig,axes=plt.subplots(1,2,figsize=(17,6));labels=[z['tag'].replace('surface_shear','S').replace('local_impingement','I').replace('ideal_local_pull','P') for z in cases];xs=np.arange(len(cases));bottom=np.zeros(len(cases));bins=[(1,1),(2,2),(3,13),(14,33),(34,113),(114,1792)]
 for (lo,hi),color in zip(bins,['#a7b5c0','#7d99ac','#4f8099','#2d7685','#de9550','#a43f54']):
  vals=np.array([sum(lo<=c['nodes']<=hi for c in z['final_free_components']) for z in cases]);axes[0].bar(xs,vals,bottom=bottom,label=str(lo) if lo==hi else f'{lo}–{hi}',color=color);bottom+=vals
 axes[0].set(ylabel='Final free components',title='Actual connected sizes; contacts excluded');axes[0].legend(title='Nodes',ncol=3,fontsize=9);axes[1].bar(xs-.18,[z['history_max_free_nodes'] for z in cases],.36,label='History largest');axes[1].bar(xs+.18,[z['last_row']['largest_detached_nodes'] for z in cases],.36,label='Final largest');axes[1].set(ylabel='Nodes',title='History versus final largest');axes[1].legend()
 for ax in axes:ax.set_xticks(xs,labels,rotation=60,ha='right',fontsize=9);ax.grid(axis='y')
 fig.tight_layout();save(fig,'Fig04_fragment_distribution','1–2 nodes are unresolved debris even when a two-node bond survives. 34 and 114 nodes exceed 200 and 300 μm volume-equivalent diameter. These are discretization measures, not clinical fragment certification.')
 fig,axes=plt.subplots(2,3,figsize=(21,10))
 for j,(ax,mode) in enumerate(zip(axes[0],MODES)):
  cc=[z for z in cases if z['params']['mode']==mode]
  if not cc:ax.text(.5,.5,'No completed case yet',ha='center',transform=ax.transAxes);ax.set_title(mode);continue
  z=max(cc,key=lambda z:z['params']['peak_pa']);h=np.genfromtxt(z['path']/'history.csv',delimiter=',',names=True,dtype=None,encoding='utf8')
  for field,label in [('W_J','External W'),('U_J','Stored U'),('K_J','Kinetic K'),('D_J','Damage D'),('Uc_J','Contact Uc'),('R_J','Residual R')]:ax.plot(h['time_s']*1e3,h[field]*1e9,label=label,lw=1.5)
  ax.set(xlabel='Time (ms)',ylabel='Energy (nJ)',title=mode.replace('_',' ')+'\n'+str(z['params']['peak_pa'])+' Pa');ax.grid()
  low=axes[1,j]
  for field,label,color in [('U_J','Stored U','C1'),('D_J','Damage D','red'),('Uc_J','Contact Uc','C4'),('R_J','Residual R','C5')]:low.plot(h['time_s']*1e3,h[field]*1e9,label=label,color=color,lw=1.5)
  low.set_yscale('symlog',linthresh=.001);low.axhline(0,c='#777',lw=.5);low.set(xlabel='Time (ms)',ylabel='U, D, Uc, R (nJ; symlog)',title='Smaller ledger terms; signed R retained');low.grid()
 axes[0,1].legend(ncol=2,fontsize=10);axes[1,1].legend(ncol=2,fontsize=10);fig.tight_layout();save(fig,'Fig05_energy_ledger','Highest tested amplitude per mode. Upper panels show the full linear energy ledger; lower panels resolve U, D, Uc and signed R on symlog axes (linear threshold 0.001 nJ). Full CSV ledgers for every amplitude accompany the figure. Contact activation offset Jc is recorded separately. R is not material dissipation.')
 categories=['NO_DAMAGE','SOFTENING_ONLY','BOND_FAILURE_NO_DETACHMENT','DETACHED_SINGLETON_DEBRIS','COHESIVE_MULTINODE_DETACHMENT','SECONDARY_FRAGMENTATION','NUMERICALLY_UNRESOLVED'];palette=['#619b83','#dab04a','#d77c4b','#b95d71','#397fa1','#82386e','#777'];fig,ax=plt.subplots(figsize=(13,6))
 for z in cases:
  j=categories.index(z['classification']);ax.scatter(z['params']['peak_pa'],MODES.index(z['params']['mode']),c=palette[j],s=160,marker='X' if j==6 else 'o',edgecolor='#222');ax.annotate(str(z['history_max_free_nodes'])+' nodes',(z['params']['peak_pa'],MODES.index(z['params']['mode'])),xytext=(0,12),textcoords='offset points',ha='center',fontsize=10)
 ca=read(O/'COHESIVE_CHUNK_AUDIT.json')
 for z in cases:
  q=ca[z['tag']]
  if q.get('same_original_members_connected_at_final') and q.get('history_largest_nodes',0)>=34:ax.scatter(z['params']['peak_pa'],MODES.index(z['params']['mode']),s=330,facecolors='none',edgecolors='#397fa1',lw=2.3)
 for name,color in zip(categories,palette):ax.scatter([],[],c=color,s=80,label=name.replace('_',' '))
 ax.scatter([],[],s=120,facecolors='none',edgecolors='#397fa1',lw=2,label='≥34-node maximum retained')
 ax.set_xscale('log');ax.set(xlabel='Peak nominal traction (Pa)',yticks=range(3),yticklabels=[m.replace('_',' ') for m in MODES],ylim=(-.5,2.7),title='Observed mechanical regimes · no acoustic feasibility implied');ax.grid();ax.legend(loc='upper center',bbox_to_anchor=(.5,-.18),ncol=2,fontsize=10);save(fig,'Fig06_mechanical_regime_map','Classification follows accepted history; unresolved numerical cases are not powdering. Labels give the history maximum free size. A blue ring means that the same original members of a ≥34-node maximum remain connected at the end. Secondary fragmentation can concern smaller debris while that large component survives.')
 dump(F/'FIGURES.json',gallery);return gallery

def main():
 cases,rows,brackets=collect();gallery=figures(cases)
 cards=''.join('<section><h2>'+g['name']+'</h2><a href="'+g['name']+'.svg"><img src="'+g['name']+'.png"></a><p>'+html.escape(g['caption'])+'</p></section>' for g in gallery)
 page='<!doctype html><meta charset="utf-8"><title>Mechanical load requirement</title><style>body{font:17px sans-serif;max-width:1600px;margin:30px auto;background:#f6f7f8;color:#243447}section{background:white;padding:24px;margin:24px 0}img{width:100%}p{line-height:1.6}</style><h1>Clot mechanical load requirement bracketing</h1><p>PRESCRIBED_TRACTION_MECHANICAL_REQUIREMENT · finite 2.5 ms observation · nominal reference-area dead loads · not validated acoustic streaming.</p>'+cards
 base='../../../verification/'+TASK+'/'+RUN+'/'
 links='<p><a href="'+base+'MECHANICAL_LOAD_REQUIREMENT_TABLE.csv">Mechanical requirement CSV</a> · <a href="'+base+'LOAD_BRACKETING_SUMMARY.csv">Brackets CSV</a> · <a href="'+base+'REVIEW_CLOT_MECHANICAL_LOAD_REQUIREMENT_BRACKETING_V1.md">中文科研报告</a> · <a href="'+base+'FINAL_STATUS.json">Final status</a></p>'
 page=page.replace('<h1>',links+'<h1>',1)
 (F/'OPEN_RESULTS.html').write_text(page);(O/'OPEN_RESULTS.html').write_text('<meta http-equiv="refresh" content="0;url=../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html">')
 print(json.dumps(dict(cases=len(cases),brackets=brackets),default=cv))
if __name__=='__main__':main()
