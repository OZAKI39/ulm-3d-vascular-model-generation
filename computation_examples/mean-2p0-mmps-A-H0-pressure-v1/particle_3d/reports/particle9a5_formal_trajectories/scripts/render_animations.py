"""Saved-trajectory replay only; deterministic smallest-ID outlet representatives."""
from pathlib import Path
import sys,json,subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
from matplotlib.lines import Line2D
sys.path.insert(0,str(Path(__file__).resolve().parent))
import render_formal as v
from particle_3d.formal_cohort_p9a5 import canonical,digest
R=v.R;D=v.D;OUT=R/'animations';OUT.mkdir(exist_ok=True)


def sphere_faces(center,radius):
 u=np.linspace(0,2*np.pi,13);w=np.linspace(0,np.pi,9)
 p=np.stack([np.outer(np.cos(u),np.sin(w)),np.outer(np.sin(u),np.sin(w)),np.outer(np.ones_like(u),np.cos(w))],axis=-1)*radius+center
 return [np.array([p[i,j],p[i+1,j],p[i+1,j+1],p[i,j+1]]) for i in range(12) for j in range(8)]


def position(a,t):
 t=np.clip(t,a[0,0],a[-1,0]);j=min(max(np.searchsorted(a[:,0],t,side='right')-1,0),len(a)-2)
 f=(t-a[j,0])/(a[j+1,0]-a[j,0]) if a[j+1,0]>a[j,0] else 0.
 return (a[j,1:4]+f*(a[j+1,1:4]-a[j,1:4]))*1e6


def main():
 s,co,rows=v.load();wall,_,surfaces=v.geometry();simple=wall.decimate_pro(.85,preserve_topology=True)
 triangles=simple.points[simple.faces.reshape(-1,4)[:,1:]]*1e6
 arrays={r['particle_id']:np.load(R/r['track_relative_path']/'trajectory.npz')['samples'] for r in rows}
 reps={o:min([r['particle_id'] for r in rows if r['status']=='COMPLETED' and r['outlet']==o],default=None) for o in ['O1','O2','O3']}
 extra=[min((r['particle_id'] for r in rows if r['status']==status),default=None) for status in ['SUPPORTED_STATIONARY','LONG_RESIDENCE_CENSORED']]
 overview=[p for p in list(reps.values())+extra if p is not None]
 jobs=[('formal_overview',overview,f'正式 {len(rows):,} 条微泡轨迹总览',True)]
 jobs.extend((o+'_representative',[pid],f'{o} 自然流经代表 · particle_id {pid}',False) for o,pid in reps.items() if pid is not None)
 records=[];frames=144;fps=18
 for name,ids,title,is_overview in jobs:
  fig=plt.figure(figsize=(12.8,7.2),facecolor='#0a1019');ax=fig.add_subplot(111,projection='3d',facecolor='#0a1019');fig.subplots_adjust(left=.01,right=.86,bottom=.09,top=.9)
  ax.add_collection3d(Poly3DCollection(triangles,facecolor='#7291a4',edgecolor='none',alpha=.075))
  v.bounds(ax,wall.points*1e6);ax.set_axis_off()
  chosen_rows=[r for r in rows if is_overview or r['particle_id'] in ids]
  for category in v.COLORS:
   paths=[]
   for r in chosen_rows:
    if v.category(r)==category:
     xyz=arrays[r['particle_id']][:,1:4]*1e6;idx=np.unique(np.linspace(0,len(xyz)-1,min(90,len(xyz)),dtype=int));paths.append(xyz[idx])
   if paths:ax.add_collection3d(Line3DCollection(paths,colors=v.COLORS[category],linewidths=.55 if is_overview else 1.7,alpha=.07 if is_overview else .55))
  artists=[];lookup={r['particle_id']:r for r in rows}
  for pid in ids:
   r=lookup[pid];a=arrays[pid];color=v.COLORS[v.category(r)]
   ax.plot(*(a[:,1:4]*1e6).T,color=color,lw=1.9,alpha=.95)
   ball=Poly3DCollection(sphere_faces(a[0,1:4]*1e6,r['radius_m']*1e6),facecolor='#f9f6db',edgecolor='none',alpha=1);ax.add_collection3d(ball);artists.append((pid,ball))
  for name_role,surf in surfaces.items():
   if name_role.startswith('OUTLET_'):
    center=surf.points.mean(axis=0)*1e6;o='O'+str(int(name_role[-2:]));ax.text(*center,o,color=v.COLORS[o],fontsize=13,fontweight='bold')
  fig.suptitle(title,color='#f3f6fb',fontsize=20,y=.97)
  handles=[Line2D([0],[0],color=v.COLORS[v.category(lookup[pid])],lw=2,label=f'{v.LABELS[v.category(lookup[pid])]}\nID {pid}') for pid in ids]
  legend=fig.legend(handles=handles,loc='center right',bbox_to_anchor=(.995,.52),frameon=False,fontsize=10)
  for text in legend.get_texts():text.set_color('#edf3f8')
  footer=fig.text(.5,.035,'',color='#d5e2ee',ha='center',fontsize=11)
  video=OUT/(name+'.mp4');writer=FFMpegWriter(fps=fps,codec='libx264',bitrate=4500,extra_args=['-pix_fmt','yuv420p','-movflags','+faststart'])
  with writer.saving(fig,str(video),dpi=100):
   for frame in range(frames):
    phase=frame/(frames-1)
    for pid,artist in artists:
     r=lookup[pid];a=arrays[pid];artist.set_verts(sphere_faces(position(a,phase*a[-1,0]),r['radius_m']*1e6))
    ax.view_init(elev=22,azim=-65+(18*np.sin(2*np.pi*phase) if is_overview else 0))
    if is_overview:footer.set_text(f'独立轨迹按归一化进度重放 {phase*100:5.1f}% · 非同时注入 · 球体使用真实半径')
    else:
     pid=ids[0];footer.set_text(f'轨迹年龄 {phase*arrays[pid][-1,0]:.4f} / {arrays[pid][-1,0]:.4f} s · 完成者中最小 ID · 不外推出口后的运动')
    writer.grab_frame(facecolor=fig.get_facecolor())
  plt.close(fig)
  subprocess.run(['ffmpeg','-v','error','-i',str(video),'-f','null','-'],check=True)
  probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-select_streams','v:0','-show_entries','stream=width,height,nb_frames,r_frame_rate,duration','-of','json',str(video)]))
  records.append(dict(file=str(video.relative_to(R)),sha256=digest(video),particle_ids=ids,all_paths_drawn=len(rows) if is_overview else len(ids),frames=frames,probe=probe,full_decode_pass=True))
  print('ANIMATION_READY',video.name,flush=True)
 manifest=dict(formal_cohort_sha256=digest(D/'FINAL_FORMAL_COHORT.json'),representatives=reps,selection='MINIMUM_COMPLETED_PARTICLE_ID_PER_NATURALLY_OBSERVED_OUTLET',
  all_outlets_observed=s['all_outlets_naturally_observed'],videos=records,actual_trajectory_data_only=True,
  overview_semantics='Independent trajectories aligned by normalized age, not simultaneous infusion',geometry_render_only_decimation=.85,original_scientific_geometry_unchanged=True)
 (D/'animation_manifest.json').write_bytes(canonical(manifest))
if __name__=='__main__':main()
