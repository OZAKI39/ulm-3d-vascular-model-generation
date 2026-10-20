"""Two 1080p videos from saved FP64 accepted states only; no trajectory interpolation."""
from common import *
import pyvista as pv
import matplotlib
import matplotlib.colors as mcolors
import pandas as pd
import subprocess,time,resource,csv
from PIL import Image,ImageDraw,ImageFont
from vtkmodules.vtkCommonCore import vtkStringArray
from figures import damage,DAMAGE
from shear_field_library import ShearField

graph=module('video_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
FPS=24;WIDTH=1920;HEIGHT=1080;HOLD_FRAMES=60
BLUE='#37779b';ANCHOR='#6d7884'
DISPLAY_DRIVE=ShearField('uniform')
STATUS={m:json.loads((O/'cases'/m/'SUMMARY.json').read_text())['status'] for m in CASES} if all((O/'cases'/m/'SUMMARY.json').exists() for m in CASES) else {}

def load_states(mode):
 case=O/'cases'/mode;paths=[case/'keys/initial.npz']+sorted((case/'keys').glob('frame_*.npz'))+[case/'keys'/('final.npz' if (case/'keys/final.npz').exists() else 'last_accepted.npz')]
 # In addition to regular states, exact onset/failure keys are eligible for nearest alignment.
 
 for key in ['first_softening','first_failure','first_secondary','first_multinode','max_multifragment','first_instant_target']:
  p=case/'keys'/(key+'.npz')
  if p.exists():paths.append(p)
 states=[dict(path=p,state=load_npz(p)) for p in paths];return sorted(states,key=lambda z:float(z['state']['time']))
def topology(ref,s,registry):
 lab,sz,att=graph.topology(ref,s['d']);free=~att[lab];rgb=np.tile(np.array([174,187,197],np.uint8),(len(lab),1));comp_ids=[]
 for label in np.flatnonzero(~att):
  ids=np.flatnonzero(lab==label);cid,color_seed=registry[tuple((ids+1).tolist())];comp_ids.append(cid);rgb[ids]=(255*mcolors.hsv_to_rgb([(color_seed*.61803398875)%1,.67,.76])).astype(np.uint8)
 return dict(labels=lab,sizes=sz,attached=att,free=free,rgb=rgb,component_ids=comp_ids)

def color_registry(mode):
 case=O/'cases'/mode;original=json.loads((case/'COMPONENT_REGISTRY.json').read_text());summary=json.loads((case/'SUMMARY.json').read_text())
 seeds={int(cid):int(cid) for cid in original};inheritance=[]
 for event in summary['secondary_events']:
  parent=int(event['parent_component_id']);children=[int(cid) for cid in event['children_component_ids']]
  dominant=max(children,key=lambda cid:original[str(cid)]['nodes'])
  if original[str(dominant)]['nodes']>original[str(parent)]['nodes']/2:
   seeds[dominant]=seeds[parent];inheritance.append(dict(time_s=event['time_s'],parent_component_id=parent,child_component_id=dominant,color_seed=seeds[parent]))
 # Color inheritance changes no member set or connectivity. At most one child
 # inherits a parent color; every other child receives its own distinct seed.
 registry={tuple(q['atom_ids']):(int(cid),seeds[int(cid)]) for cid,q in original.items()}
 return registry,dict(color_seed_by_exact_component_id=seeds,dominant_child_inheritance=inheritance)
def lines(points,bonds):
 mesh=pv.PolyData(points)
 # PolyData(points) also creates one vertex cell per input point. A line actor
 # must not redisplay all unselected material nodes through those extra cells.
 mesh.verts=np.empty(0,dtype=np.int64)
 mesh.lines=np.c_[np.full(len(bonds),2),bonds].ravel()
 assert mesh.n_verts==0 and mesh.n_lines==len(bonds)
 return mesh

def add_free_nodes(pl,xyz,top,prefix,display_mask=None):
 # All free nodes remain visible. Fine-debris glyphs use a muted class color;
 # components of >=3 nodes retain original-member lineage colors.
 sizes=top['sizes'][top['labels']]
 for name,mask,point_size in [('free_fines',sizes<=2,3.),('free_small',(sizes>=3)&(sizes<20),4.),('free_target',(sizes>=20)&(sizes<=200),7.),('free_large',sizes>200,5.)]:
  selected=top['free']&mask
  if display_mask is not None:selected&=display_mask
  if not selected.any():continue
  args=dict(color='#b29a72') if name=='free_fines' else dict(scalars=top['rgb'][selected],rgb=True)
  pl.add_points(xyz[:len(sizes)][selected],**args,point_size=point_size,render_points_as_spheres=True,show_scalar_bar=False,name=prefix+name)

def fitted_parallel_scale(bounds,camera):
 # Fit the projected physical bounding box to the actual 960 × 648 viewport.
 # Both panels retain this one fixed scale for the full real trajectory.
 corners=np.array([(x,y,z) for x in bounds[:2] for y in bounds[2:4] for z in bounds[4:6]])
 forward=np.asarray(camera[1])-np.asarray(camera[0]);forward=forward/np.linalg.norm(forward)
 right=np.cross(forward,np.asarray(camera[2]));right=right/np.linalg.norm(right);up=np.cross(right,forward)
 aspect=(WIDTH/2)/(.60*HEIGHT)
 return float(1.36*max(np.ptp(corners@up)/2,np.ptp(corners@right)/(2*aspect)))

def base_view(pl,index,ref,bounds,camera):
 pl.subplot(0,index);pl.renderer.SetViewport(index*.5,.22,(index+1)*.5,.82);pl.set_background('#fcfcfc');pl.camera_position=camera;pl.enable_parallel_projection()
 pl.camera.parallel_scale=fitted_parallel_scale(bounds,camera)
 span=max(bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]);pl.camera.clipping_range=(.01*span,20*span)
 axes=pl.show_bounds(bounds=bounds,grid='all',location='outer',ticks='outside',xtitle='X (um)',ytitle='Y (um)',ztitle='Z (um)',font_size=14,color='#425565',n_xlabels=4,n_ylabels=4,n_zlabels=4,all_edges=True,fmt='%.0f',use_3d_text=False)
 # The lowest Y label collides with the highest X label at their shared corner.
 # Keep the physical bounds/ticks unchanged and omit only that redundant label.
 ylabels=vtkStringArray()
 for k,value in enumerate(np.linspace(bounds[2],bounds[3],4)):ylabels.InsertNextValue('' if k==0 else f'{value:.0f}')
 axes.SetAxisLabels(1,ylabels)
 pl.add_mesh(pv.Box(bounds=bounds).outline(),color='#9aa6af',line_width=.7,name=f'panel{index}_fixed_outline')
 pl.add_axes(interactive=False,line_width=2)
def add_nodes(pl,ref,s,top,kind,prefix):
 n=int(ref['clot_count']);xyz=s['x']*1e6
 pl.add_points(xyz[n:],color=ANCHOR,point_size=4,render_points_as_spheres=True,opacity=.45,name=prefix+'anchors')
 if kind=='damage':
  pl.add_points(xyz[:n],scalars=damage(ref,s['d'])[:n],cmap=DAMAGE,clim=(0,1),point_size=3,render_points_as_spheres=True,show_scalar_bar=False,scalar_bar_args=dict(title='Mean incident damage d',vertical=False,position_x=.12,position_y=.045,width=.72,height=.07,title_font_size=16,label_font_size=14,n_labels=3),name=prefix+'material')
  add_free_nodes(pl,xyz,top,prefix)
 elif kind=='fragment':
  pl.add_points(xyz[:n][~top['free']],color='#b4c1ca',point_size=4,render_points_as_spheres=True,opacity=.30,name=prefix+'body')
  add_free_nodes(pl,xyz,top,prefix)
 else:
  pl.add_points(xyz[:n],color='#8a9eaa',point_size=3,render_points_as_spheres=True,opacity=.30,name=prefix+'material')
  add_free_nodes(pl,xyz,top,prefix)
def overlay(image,video,index,timeline,states,tops,ref,final_hold=False):
 im=Image.fromarray(image).convert('RGB');d=ImageDraw.Draw(im);font='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf';regular=ImageFont.truetype(font,25);small=ImageFont.truetype(font,21);large=ImageFont.truetype(font,32)
 d.rectangle((0,0,WIDTH,194),fill='#fcfcfc');d.rectangle((0,842,WIDTH,HEIGHT),fill='#fcfcfc')
 titles={1:'Uniform vs biphasic material · native fracture',2:'Biphasic material · true damage and cohesive fragments',3:'True free-component evolution · no visual regrouping'}
 d.text((40,14),titles[video],font=large,fill='#243d4e')
 d.text((40,59),'Accepted coordinates ×1 | uniform +X shear, 1000 Pa nominal peak | seed 20261018 | exploratory, uncalibrated',font=small,fill='#536675')
 labels=['U · Uniform reference','H · Biphasic material'] if video in [1,3] else ['H · Damage and surviving bonds','H · Frozen phases and true free components']
 n=int(ref['clot_count'])
 for j,m in enumerate(['uniform','biphasic'] if video in [1,3] else ['biphasic','biphasic']):
  s=states[m];top=tops[m];free=top['sizes'][~top['attached']];target=int(np.sum((free>=20)&(free<=200)));fine=int(free[free<=2].sum());total=int(free.sum());largest=int(free.max(initial=0));ratio=f'{100*fine/total:.1f}%' if total else 'undefined'
  nominal=1000*float(DISPLAY_DRIVE.envelope(float(s['time'])))
  panel_label=labels[j]
  if video==2 and j==1:panel_label='H · True free components >=20 nodes' if np.any(free>=20) else 'H · Attached material phases'
  x=j*960+35;d.text((x,115),panel_label+(' · INCOMPLETE' if STATUS.get(m)!='COMPLETED' else ''),font=regular,fill='#293f4d');d.text((x,151),f'Actual t = {float(s["time"])*1e3:.5f} ms | nominal +X: {nominal:.1f} Pa',font=regular,fill='#293f4d')
  d.rectangle((j*960+15,906,(j+1)*960-15,1058),fill='#fcfcfc')
  if video==2 and j==1:
   detail=(f'Displayed free >=20: {np.sum(free>=20)} | largest {largest} nodes' if np.any(free>=20) else 'No free body >=20; attached phases shown')
   text=[f'TARGET_FRAGMENT_COUNT_20_200 = {target}',detail,('TARGET_FRAGMENT_NOT_OBSERVED' if target==0 else 'All internal live bonds shown in displayed components')]
  else:text=[f'Free components {len(free)} | largest {largest} nodes | detached {100*total/n:.2f}%',f'TARGET_FRAGMENT_COUNT_20_200 = {target} | fines {ratio}',('TARGET_FRAGMENT_NOT_OBSERVED' if target==0 else 'True target-band connected material components')]
  for k,line in enumerate(text):d.text((x,917+k*31),line,font=small,fill='#9f475b' if k==2 and target==0 else '#435a6a')
 
 # Manual scalar bars avoid renderer/viewport scaling and remain readable at 1080p.
 for j in range(2):
  cmap=DAMAGE if video==1 else None
  if cmap is not None:
   left=j*960+250;right=left+460;d.text((left,842),'Attached-node mean damage d' if video==1 else 'Frozen reference q',font=small,fill='#435a6a')
   for px in range(460):d.line((left+px,868,left+px,880),fill=tuple(int(c*255) for c in cmap(px/459)[:3]))
   for k,label in enumerate(['0','0.5','1'] if video==1 else ['−1','0','1']):d.text((left+k*230-8,882),label,font=small,fill='#435a6a')
 if video==2:
  for j,(label,color) in enumerate([('Intact','#347e91'),('Softened','#d48a25'),('Failed trace','#b23b55')]):
   d.line((40,211+j*28,69,211+j*28),fill=color,width=4);d.text((78,198+j*28),label,font=small,fill=color)
 if video==2:
  has_cohesive=np.any(tops['biphasic']['sizes'][~tops['biphasic']['attached']]>=20)
  d.text((1030,854),'Whole free components >=20; lineage colors' if has_cohesive else 'Attached phases: weak (blue) / strong (red)',font=small,fill='#435a6a')
  d.text((1030,883),'Full scene and fines are retained on the left.' if has_cohesive else 'No free component >=20 is present at this time.',font=small,fill='#435a6a')
 d.text((40,1029),'Anchors: gray. Fines: small tan points. A >50% child inherits its parent color. All original IDs remain recorded; glyphs are not shells.',font=small,fill='#536675')
 d.text((40,1058),'Supported-facet loading is a proxy; graph cohesion alone does not certify hydrodynamic stability.',font=ImageFont.truetype(font,18),fill='#7c5361')
 return np.array(im)
def render(which,selected,ref,field,bounds,camera,registries,only_frame=False,smoke=False):
 starts=time.perf_counter();cpu=time.process_time();name={1:'VIDEO_01_UNIFORM_VS_BIPHASIC',2:'VIDEO_02_TRUE_COHESIVE_FRAGMENTS',3:'VIDEO_03_TRUE_FRAGMENT_EVOLUTION'}[which];out=F/'videos';out.mkdir(exist_ok=True,parents=True)
 if smoke:out=out/'render_smoke';out.mkdir(exist_ok=True)
 cmd=['ffmpeg','-hide_banner','-loglevel','error','-y','-f','rawvideo','-vcodec','rawvideo','-pix_fmt','rgb24','-s',f'{WIDTH}x{HEIGHT}','-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','medium','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart',str(out/(name+'.mp4'))]
 log=(out/(name+'.ffmpeg.log')).open('w');enc=None if only_frame else subprocess.Popen(cmd,stdin=subprocess.PIPE,stderr=log)
 records=[];last_image=None;pl=pv.Plotter(shape=(1,2),off_screen=True,window_size=(WIDTH,HEIGHT),border=False);pl.ren_win.SetMultiSamples(4)
 try:
  seq=[len(selected)-1] if only_frame else range(len(selected))
  for index in seq:
   row=selected[index];s={m:row[m]['state'] for m in CASES};tops={m:topology(ref,s[m],registries[m]) for m in CASES};panels=CASES if which in [1,3] else ['biphasic','biphasic']
   for j,m in enumerate(panels):
    pl.subplot(0,j);state=s[m];top=tops[m];xyz=state['x']*1e6;prefix=f'panel{j}_'
    # Keep the camera, bounds and orientation widgets alive across frames.
    # Clearing and recreating those VTK widgets crashed the offscreen renderer.
    for actor in ['anchors','material','free','free_fines','free_small','free_target','free_large','body','q','live_bonds','soft_bonds','failed_bonds','internal_bonds']:pl.renderer.remove_actor(prefix+actor,render=False)
    if which==2 and j==1:
     n=int(ref['clot_count']);cohesive=top['free']&(top['sizes'][top['labels']]>=20);shown=cohesive if cohesive.any() else ~top['free']
     if shown.any():pl.add_points(xyz[:n][shown],scalars=field['material_phase'][:n][shown],cmap=['#7eabba','#b83e62'],clim=(0,1),point_size=3,render_points_as_spheres=True,show_scalar_bar=False,name=prefix+'q')
     pl.add_points(xyz[int(ref['clot_count']):],color=ANCHOR,point_size=4,opacity=.4,name=prefix+'anchors')
    else:
     add_nodes(pl,ref,state,top,'damage' if which==1 else 'fragment' if which==3 else 'bonds',prefix)
     b=ref['bonds'];near=(ref['r0']<=1.01*DX)&~ref['initial_cut']
     if which==2:
      for name_bond,sel,color,opacity,width in [('live_bonds',near&(state['d']==0),'#347e91',.28,1.),('soft_bonds',near&(state['d']>0)&(state['d']<1),'#d48a25',.8,2.),('failed_bonds',near&(state['d']==1),'#b23b55',.12,1.)]:
       if sel.any():pl.add_mesh(lines(xyz,b[sel]),color=color,opacity=opacity,line_width=width,lighting=False,name=prefix+name_bond)
     if which in [2,3]:
      inside=np.r_[top['free'],np.zeros(len(ref['X'])-len(top['free']),bool)];sel=inside[b].all(1)&(state['d']<1)&(ref['bond_types']==1)
      if sel.any():pl.add_mesh(lines(xyz,b[sel]),color='#254d63',line_width=.8,opacity=.16,lighting=False,name=prefix+'internal_bonds')
    if which==2 and j==1:
     n=int(ref['clot_count']);free=top['free']
     add_free_nodes(pl,xyz,top,prefix,display_mask=cohesive)
     b=ref['bonds'];inside=np.r_[cohesive,np.zeros(len(ref['X'])-n,bool)];sel=inside[b].all(1)&(state['d']<1)&(ref['bond_types']==1)
     if sel.any():pl.add_mesh(lines(xyz,b[sel]),color='#254d63',line_width=.8,opacity=.16,lighting=False,name=prefix+'internal_bonds')
    if index==seq[0]:base_view(pl,j,ref,bounds,camera)

   for j in range(2):
    key='body' if which==3 else 'q' if which==2 and j==1 else 'material'
    if which==2 and j==1 and not shown.any():continue # A physically empty display selection is explicit, never populated with fake nodes.
    assert f'panel{j}_{key}' in pl.renderers[j].actors,('missing material actor',which,index,j)
   img=pl.screenshot(return_img=True);img=overlay(img,which,index,row['requested_time_s'],s,tops,ref);last_image=img
   if index in [0,len(selected)//2,len(selected)-1]:Image.fromarray(img).save(out/(name+f'_frame_{index:03d}.png'))
   if index==len(selected)-1:Image.fromarray(img).save(out/(name+'_COVER.png'))
   if enc:enc.stdin.write(img.tobytes())
   records.append(dict(frame=index,playback_s=index/FPS,requested_time_s=row['requested_time_s'],uniform_actual_time_s=float(s['uniform']['time']),biphasic_actual_time_s=float(s['biphasic']['time']),uniform_source=str(row['uniform']['path'].relative_to(O)),biphasic_source=str(row['biphasic']['path'].relative_to(O)),right_selection=('FREE_COMPONENTS_GE20' if cohesive.any() else 'ATTACHED_PHASES') if which==2 else 'ALL_MATERIAL_NODES',right_selected_material_nodes=int(shown.sum()) if which==2 else int(ref['clot_count']),state_interpolation=False,hold=False))
   if index%24==0:print('VIDEO_FRAME',which,index,flush=True)
  if enc:
   for k in range(HOLD_FRAMES):
    enc.stdin.write(last_image.tobytes());q=records[-1].copy();q.update(frame=len(selected)+k,playback_s=(len(selected)+k)/FPS,hold=True);records.append(q)
   enc.stdin.close();code=enc.wait();assert code==0,code
 finally:
  pl.close();log.close()
  if enc and enc.poll() is None:enc.kill();enc.wait()
 if only_frame:return
 with (out/(name+'_FRAME_TIMES.csv')).open('w') as f:w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
 (out/(name+'_FFMPEG_COMMAND.txt')).write_text(' '.join(cmd)+'\nInput is raw RGB frames piped by render_videos.py, not interpolated trajectories.\n')
 probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-show_format','-of','json',str(out/(name+'.mp4'))],text=True));v=probe['streams'][0];assert v['codec_name']=='h264' and v['width']==WIDTH and v['height']==HEIGHT and int(v['nb_frames'])==len(records)
 return dict(file=name+'.mp4',sha256=sha(out/(name+'.mp4')),frames=len(records),native_timeline_frames=len(selected),final_hold_seconds=HOLD_FRAMES/FPS,fps=FPS,resolution=[WIDTH,HEIGHT],codec='H.264',material_actor_presence_checked_every_frame=True,wall_s=time.perf_counter()-starts,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,ffprobe=probe)
def main():
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--preview',action='store_true');p.add_argument('--smoke',action='store_true');args=p.parse_args()
 ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'BIPHASIC_FIELD.npz');states={m:load_states(m) for m in CASES};events=[float(z['state']['time']) for arr in states.values() for z in arr if z['path'].stem in ['first_softening','first_failure','first_secondary','first_multinode','max_multifragment','first_instant_target']];last_time=max(float(arr[-1]['state']['time']) for arr in states.values());requested=np.unique(np.r_[np.linspace(0,last_time,145),events]);selected=[dict(requested_time_s=t,**{m:min(states[m],key=lambda z:abs(float(z['state']['time'])-t)) for m in CASES}) for t in requested]
 points=np.concatenate([z['state']['x'] for arr in states.values() for z in arr])*1e6;lo=points.min(0);hi=points.max(0);pad=np.maximum((hi-lo)*.06,25);lo-=pad;hi+=pad;bounds=np.array([lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]]);center=(lo+hi)/2;span=np.max(hi-lo);camera=[center+np.array([1.4,-2.9,1.0])*span,center,[0,0,1]]
 registries={};color_evidence={}
 for m in CASES:registries[m],color_evidence[m]=color_registry(m)
 dump(O/'VIDEO_COMPONENT_COLOR_LINEAGE.json',dict(rule='Exact surviving-graph member sets retain their registry IDs. For components of at least 3 nodes, only a child containing strictly more than half the parent nodes inherits its color; all other children use new color seeds. Fine debris of 1–2 nodes uses a muted tan class color while all member IDs remain recorded. No proximity regrouping, no change to component definitions, no cross-case identity correspondence.',cases=color_evidence))
 results=[]
 if args.smoke:selected=[selected[i] for i in np.linspace(0,len(selected)-1,5).astype(int)]
 for which in [1,2]:results.append(render(which,selected,ref,field,bounds,camera,registries,args.preview,args.smoke))
 if not args.preview and not args.smoke:
  dump(O/'VIDEO_GENERATION.json',dict(status='PASS',videos=results,fixed_global_bounds_um=bounds,camera=camera,fixed_parallel_scale_um=fitted_parallel_scale(bounds,camera),camera_fit='Orthographic projected full bounding box with 36% viewport margin; same fixed scale in both panels for all frames',axis_label_policy='Physical bounds and ticks unchanged; minimum Y text omitted to prevent shared-corner collision with maximum X text.',video2_right_selection='Whole true free components >=20 nodes, or attached phase view before such a body exists. Left panorama retains every material node and all fines. Same fixed camera and physical scale in both panels.',glyph_point_sizes_px={'fines_1_2':3,'small_3_19':4,'target_20_200':7,'large_gt200':5},actual_displacement_scale=1,actual_accepted_states_only=True,alignment='Nearest saved native accepted state to common requested time; each panel and CSV gives its own actual timestamp. No trajectory interpolation.',max_time_alignment_difference_s=max(abs(float(row['uniform']['state']['time'])-float(row['biphasic']['state']['time'])) for row in selected),display_bond_threshold_m=1.01*DX,topology='ALL original bulk/interface live bonds; anchors determine attachment',fragment_colors='Components of at least 3 nodes use original-member lineage colors; a strictly >50% dominant child inherits its parent color. Fine debris of 1–2 nodes uses small tan glyphs. All exact member sets retain their registry IDs. Mapping: VIDEO_COMPONENT_COLOR_LINEAGE.json. Colors across U/H do not imply matched material.',failed_bond_render='Translucent red current-endpoint traces, not load-carrying solid bonds',glyphs='Display point glyphs, no reconstructed shells; 1–2 node debris unresolvable as 3D continuum fragments'))
  print('VIDEOS_DONE',flush=True)
if __name__=='__main__':main()
