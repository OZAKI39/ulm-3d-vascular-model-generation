"""Diagnostic figures from measured CSVs/raw facets only; no field smoothing."""
import os
os.environ.setdefault('LP_NUM_THREADS','2');os.environ.setdefault('VTK_SMP_MAX_THREADS','2')
from pathlib import Path
import csv,json,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyvista as pv
V=Path(__file__).resolve().parents[1];FIG=V/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
COLORS=['#2966a2','#d27d21','#27834b']
def rows(p):return list(csv.DictReader(Path(p).open()))
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(FIG/(name+'.'+ext),dpi=240,bbox_inches='tight')
 plt.close(fig)
def pipe():
 finest='pipe_nr16' if (V/'stage2/pipe_nr16/reports/pipe_validation.csv').exists() else 'pipe_nr16_cpu_mpi8';cases=['pipe_nr4','pipe_nr8',finest];fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
 for case,col in zip(cases,COLORS):
  pp=V/'stage2'/case/'reports/velocity_profile.csv'
  if not pp.exists():continue
  rr=rows(pp);xx=np.array([float(r['x_um']) for r in rr]);yy=np.array([float(r['velocity_z_m_s'])*1e3 for r in rr]);ok=np.array([int(r['valid']) for r in rr],bool);ax.plot(xx[ok],yy[ok],color=col,label=case)
 xx=np.linspace(-4,4,401);ax.plot(xx,4*(1-(xx/4)**2),'k--',lw=1.5,label='Continuous Poiseuille')
 ax.set(xlabel='x at y=0, z=20 [um]',ylabel='Axial velocity [mm/s]',title='Actual CFD velocity profiles | common physical section');ax.legend();ax.grid(alpha=.2);save(fig,'Figure_01_pipe_velocity')
 data=rows(V/'data/pipe_validation.csv');fig,axs=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
 for kind,mark,ls in [('CFD','o','-'),('analytic_nodal','s','--')]:
  rr=[next(r for r in data if r['case']==case and r['field_type']==kind) for case in cases if any(r['case']==case and r['field_type']==kind for r in data)];xx=[float(r['radial_step_um']) for r in rr]
  for ax,key in zip(axs,['wss_mean_signed_error_pct','velocity_volume_L2_error_pct','central_delta_p_signed_error_pct']):ax.plot(xx,[abs(float(r[key])) for r in rr],marker=mark,ls=ls,label=kind)
 for ax,title,limit in zip(axs,['Central mean WSS error','Central velocity L2 error','Central pressure-drop error'],[5,1,1]):
  ax.axhline(limit,color='gray',ls=':',label=f'Work target {limit}%');ax.set(xlabel='Radial interval [um]',ylabel='Absolute relative error [%]',title=title);ax.grid(alpha=.2);ax.invert_xaxis();ax.legend(fontsize=8)
 save(fig,'Figure_02_pipe_refinement')
def vessel():
 rr=rows(V/'data/vessel_region_wss.csv');cases=['vessel_baseline','vessel_medium'];reg=['real_ROI','J1_r5um','J2_r5um','bend_r5um','worst_tet_r5um','O1_extension','O2_extension','O3_extension'];fig,axs=plt.subplots(2,3,figsize=(16,9),layout='constrained')
 for ax,key in zip(axs.flat,['mean_Pa','p05_Pa','p95_Pa','below1_area_pct','below2_area_pct','above30_area_pct']):
  for k,(case,col) in enumerate(zip(cases,COLORS)):
   if not any(r['case']==case for r in rr):continue
   vals=[float(next(r[key] for r in rr if r['case']==case and r['region']==rgn)) for rgn in reg];ax.plot(np.arange(len(reg)),vals,'o-',color=col,label=case)
  ax.set(xticks=np.arange(len(reg)),xticklabels=[r.replace('_extension',' ext').replace('_r5um','') for r in reg],ylabel='Area fraction [%]' if 'area_pct' in key else 'WSS [Pa]',title=key);ax.tick_params(axis='x',rotation=35);ax.grid(alpha=.2)
 axs[0,0].legend(fontsize=8);save(fig,'Figure_03_vessel_region_metrics')
 meshes={};shown_max=0.
 for case in cases:
  mesh=pv.read(V/'stage3'/case/'wss/data/wall_wss_si.vtp');mesh.points*=1e6;meshes[case]=mesh;cent=mesh.cell_centers().points
  sel=np.logical_or.reduce([np.linalg.norm(cent-np.array(ctr),axis=1)<8.5 for ctr in [[92,49,111],[130.04,82.04,87.18]]])
  shown_max=max(shown_max,float(np.asarray(mesh.cell_data['WSS_raw_Pa'])[sel].max()))
 vmax=max(55.,float(5*np.ceil(shown_max/5)))
 (FIG/'Figure_04_metadata.json').write_text(json.dumps(dict(cases=cases,field='WSS_raw_Pa',association='cell',unit='Pa',shared_linear_limits=[0,vmax],maximum_in_selected_local_meshes_Pa=shown_max,lighting=False,smoothing=False,normalization=False,camera_direction=[1,-2.2,1.1],camera_up=[0,0,1],parallel_scale_um=6.,selection_radius_um=8.5,statistics_radius_um=5.),indent=2)+'\n')
 def render(mesh,ctr):
  center=np.array(ctr);cent=mesh.cell_centers().points;local=mesh.extract_cells(np.linalg.norm(cent-center,axis=1)<8.5).extract_surface(algorithm='dataset_surface');plot=pv.Plotter(off_screen=True,window_size=(1000,850));plot.set_background('white');plot.add_mesh(local,scalars='WSS_raw_Pa',preference='cell',cmap='viridis',clim=(0,vmax),lighting=False,show_edges=True,edge_color='#303c45',line_width=.35,show_scalar_bar=False);direction=np.array([1,-2.2,1.1]);direction/=np.linalg.norm(direction);plot.camera.position=center+250*direction;plot.camera.focal_point=center;plot.camera.up=(0,0,1);plot.enable_parallel_projection();plot.camera.parallel_scale=6.;im=plot.screenshot();plot.close();return im
 fig,axs=plt.subplots(2,len(cases),figsize=(11,10));fig.subplots_adjust(wspace=.02,hspace=.03,right=.92)
 for j,case in enumerate(cases):
  p=V/'stage3'/case/'wss/data/wall_wss_si.vtp'
  if not p.exists():continue
  mesh=meshes[case]
  for i,(name,ctr) in enumerate([('J1',[92,49,111]),('J2',[130.04,82.04,87.18])]):axs[i,j].imshow(render(mesh,ctr));axs[i,j].axis('off');axs[i,j].set_title(name+' | '+case)
 cax=fig.add_axes([.94,.18,.015,.64]);fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,vmax),cmap='viridis'),cax=cax,label='Raw facet WSS [Pa] | shared linear scale');save(fig,'Figure_04_J1_J2_raw_mesh')
 profile=rows(V/'data/outlet2_profile.csv');adj=rows(V/'data/adjacent_jump_summary.csv');fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained')
 for case,col in zip(cases,COLORS):
  aa=[r for r in profile if r['case']==case];axs[0].plot([float(r['s_mid_um']) for r in aa],[float(r['mean_Pa']) for r in aa],'o-',ms=3,color=col,label=case)
  bb=[r for r in adj if r['case']==case and r['region']=='known_jump_r1um']
  if bb:axs[1].plot(['Median','P95','Maximum'],[float(bb[0][k]) for k in ['jump_p50_Pa','jump_p95_Pa','jump_max_Pa']],'o-',color=col,label=case)
 path=json.loads((V.parent/'wss_audit/data/outlet2_path.json').read_text());axs[0].axvline(path['real_cut_arc_um'],c='gray',ls='--',label='Artificial extension begins');axs[0].set(xlabel='Fixed path distance from J1 [um]',ylabel='Area-mean WSS [Pa]',title='O2 fixed path: 1 um bins');axs[1].set(ylabel='Absolute neighbor-facet jump [Pa]',title='Within 1 um of original J1 jump location')
 for ax in axs:ax.legend(fontsize=8);ax.grid(alpha=.2)
 save(fig,'Figure_05_vessel_local_response')
def boundary():
 flows=rows(V/'data/boundary_flow_responses.csv');resp=rows(V/'data/boundary_region_responses.csv');fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained');roles=['OUTLET_01','OUTLET_02','OUTLET_03'];reg=['J1_r5um','J2_r5um','O1_extension','O2_extension','O3_extension']
 for col,sign,label in [('#2966a2','minus','O2 pressure -1%'),('#d27d21','plus','O2 pressure +1%')]:
  xx=np.arange(3);vals=[float(next(r[sign+'_relative_change_pct'] for r in flows if r['boundary']==role and r['metric']=='Q_outward_m3_s')) for role in roles];axs[0].plot(xx,vals,'o-',color=col,label=label)
  for metric,ls in [('mean_Pa','-'),('p05_Pa','--'),('p95_Pa',':')]:vals=[float(next(r[sign+'_relative_change_pct'] for r in resp if r['region']==region and r['metric']==metric)) for region in reg];axs[1].plot(np.arange(len(reg)),vals,'o',ls=ls,color=col,label=label+' '+metric)
 axs[0].set(xticks=np.arange(3),xticklabels=['O1','O2','O3'],ylabel='Flow change from fixed-grid baseline [%]',title='Fixed mesh: only O2 pressure changes');axs[1].set(xticks=np.arange(len(reg)),xticklabels=['J1','J2','O1 ext','O2 ext','O3 ext'],ylabel='WSS metric change from fixed-grid baseline [%]',title='Raw area-weighted WSS response')
 for ax in axs:ax.axhline(0,color='gray',lw=.7);ax.grid(alpha=.2);ax.legend(fontsize=7)
 save(fig,'Figure_06_boundary_sensitivity')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--part',choices=['pipe','vessel','boundary','all'],default='all');a=p.parse_args()
 for name,fn in [('pipe',pipe),('vessel',vessel),('boundary',boundary)]:
  if a.part in [name,'all']:
   if name=='boundary' and (V/'stage4/user_cancellation.json').exists():
    print('Boundary sensitivity figure SKIPPED_BY_USER; no pressure CFD exists');continue
   fn()
