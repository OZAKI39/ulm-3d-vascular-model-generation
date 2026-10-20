"""Native accepted states, fixed Z-up cameras, actual volume markers and bonds."""
from common import *
from analyze_result import topology,state_paths,choose,CASE,OLD
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont
import subprocess,argparse,hashlib
WIDTH=1920;HEIGHT=1080;FPS=24
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
PALETTE=np.array([[146,24,48],[202,55,67],[111,12,38],[224,86,81],[172,37,66],[190,70,88],[128,33,58]],np.uint8)

def state(path):
 with np.load(path) as q:return {k:q[k].copy() for k in ['x','d','atom_ids','time','step']}
def line_mesh(points,edges):
 ids,inverse=np.unique(edges,return_inverse=True);mesh=pv.PolyData(points[ids]);mesh.verts=np.empty(0,np.int64);mesh.lines=np.c_[np.full(len(edges),2),inverse.reshape(-1,2)].ravel();return mesh
def grid_points(bounds):
 lo=np.array(bounds)[::2];hi=np.array(bounds)[1::2];out=[]
 for axis in range(3):
  other=[j for j in range(3) if j!=axis]
  for side in [lo[axis],hi[axis]]:
   for along,across in [other,other[::-1]]:
    for value in np.linspace(lo[across],hi[across],4):
     a=lo.copy();b=hi.copy();a[axis]=b[axis]=side;a[across]=b[across]=value;out.extend([a,b])
 return np.array(out)
def setup(bounds,panels):
 pl=pv.Plotter(off_screen=True,window_size=(WIDTH,HEIGHT),shape=(1,panels));lo=np.array(bounds)[::2];hi=np.array(bounds)[1::2];center=(lo+hi)/2;span=max(hi-lo);direction=np.array([1.2,-3.,-.6]);corners=np.array([(x,y,z) for x in bounds[:2] for y in bounds[2:4] for z in bounds[4:]])
 for j in range(panels):
  pl.subplot(0,j);pl.set_background('white');pl.renderer.SetViewport(j/panels+.02,.14,(j+1)/panels-.02,.9);pl.camera_position=[center+span*direction,center,(0,0,1)];pl.enable_parallel_projection();forward=-direction/np.linalg.norm(direction);right=np.cross(forward,[0,0,1]);right/=np.linalg.norm(right);up=np.cross(right,forward);aspect=(WIDTH*(1/panels-.04))/(HEIGHT*.76);pl.camera.parallel_scale=float(1.32*max(np.ptp(corners@up)/2,np.ptp(corners@right)/(2*aspect)));pl.camera.clipping_range=(span*.01,span*20)
  g=grid_points(bounds);pl.add_mesh(line_mesh(g,np.arange(len(g)).reshape(-1,2)),color='#76848f',opacity=.5,line_width=.8,lighting=False,name=f'grid{j}')
 return pl
def scene(pl,j,ref,s,draw_edges,opacity=1.):
 pl.subplot(0,j);n=int(ref['clot_count']);N=len(ref['X']);lab,attached=topology(ref,s['d']);rgb=np.tile([171,27,49],(N,1)).astype(np.uint8);rgb[n:]=[124,145,160];counts=np.bincount(lab);free=~attached[lab];order=np.argsort(lab,kind='stable');starts=np.r_[0,np.cumsum(counts)]
 for label in np.flatnonzero(~attached):
  ids=np.sort(order[starts[label]:starts[label+1]]);code=int.from_bytes(hashlib.sha256((ids+1).astype('<i8').tobytes()).digest()[:4],'little');rgb[ids]=PALETTE[code%len(PALETTE)]
 points=pv.PolyData(s['x']*1e3);points['radius']=np.cbrt(3*ref['V']/(4*np.pi))*1e3;points['rgb']=rgb
 glyph=points.glyph(geom=pv.Sphere(radius=1.,theta_resolution=10,phi_resolution=10),scale='radius',orient=False)
 pl.add_mesh(glyph,scalars='rgb',rgb=True,name=f'nodes{j}',show_scalar_bar=False,smooth_shading=True,ambient=.25,diffuse=.75,specular=.12,opacity=opacity)
 b=ref['bonds'];target=free&(counts[lab]>=2);alive=(s['d']<1)&target[b].all(1)
 if draw_edges and alive.any():pl.add_mesh(line_mesh(s['x']*1e3,b[alive]),name=f'edges{j}',color='#510e25',line_width=.8,opacity=.65,lighting=False)
 else:pl.remove_actor(f'edges{j}',reset_camera=False)
 return dict(time_s=float(s['time']),step=int(s['step']),free_nodes=int(free[:n].sum()),free_components=int((~attached).sum()),N20_200=int(np.sum(~attached&(counts>=20)&(counts<=200))),shown_live_internal_bonds=int(alive.sum()) if draw_edges else 0)
def annotate(pl,raw,bounds,stats,names,panels):
 im=Image.new('RGB',(WIDTH,HEIGHT),'white');source=Image.fromarray(raw).convert('RGB');draw=ImageDraw.Draw(im);font=lambda size:ImageFont.truetype(FONT,size)
 for j in range(panels):
  box=(round(WIDTH*(j/panels+.02)),round(HEIGHT*.10),round(WIDTH*((j+1)/panels-.02)),round(HEIGHT*.86));im.paste(source.crop(box),box);pl.subplot(0,j)
  draw.text((WIDTH*(j+.5)/panels,28),names[j],font=font(30),fill='#202832',anchor='mt');draw.text((WIDTH*(j+.5)/panels,70),f't = {stats[j]["time_s"]*1e3:.5f} ms',font=font(26),fill='#202832',anchor='mt')
  def project(p):
   pl.renderer.SetWorldPoint(*p,1);pl.renderer.WorldToDisplay();x,y,_=pl.renderer.GetDisplayPoint();return np.array([x,HEIGHT-y])
  lo=np.array(bounds)[::2];hi=np.array(bounds)[1::2]
  probe=(lo+hi)/2;probe2=probe.copy();probe2[2]+=1.
  z_error=float(abs(project(probe2)[0]-project(probe)[0]))
  assert z_error<1e-7,('World Z must stay vertical',z_error)
  stats[j]['world_Z_horizontal_error_px']=z_error
  for axis in range(3):
   other=[k for k in range(3) if k!=axis];candidates=[]
   for a in [lo[other[0]],hi[other[0]]]:
    for b in [lo[other[1]],hi[other[1]]]:
     point=(lo+hi)/2;point[other]=[a,b];candidates.append(point)
   # Bottom silhouettes for horizontal axes, right silhouette for vertical Z.
   anchor=max(candidates,key=lambda p:project(p)[0 if axis==2 else 1]).copy()
   for value in np.linspace(lo[axis],hi[axis],4):
    if axis==1 and value in (lo[axis],hi[axis]):continue # avoid shared X/Y and Y/Z corner labels
    p=anchor.copy();p[axis]=value;xy=project(p)+np.array([0,14] if axis<2 else [18,0]);draw.text(tuple(xy),f'{value:.2f}',font=font(23),fill='#384651',anchor='mt' if axis<2 else 'lm')
   p=anchor.copy();p[axis]=(lo[axis]+hi[axis])/2;xy=project(p)+np.array([0,50] if axis==0 else [-28,72] if axis==1 else [65,0]);draw.text(tuple(xy),'XYZ'[axis]+' (mm)',font=font(26),fill='#263540',anchor='mt' if axis<2 else 'lm')
  draw.text((WIDTH*(j+.5)/panels,940),f'20–200 nodes: {stats[j]["N20_200"]} free components',font=font(25),fill='#642136',anchor='mt')
 draw.text((WIDTH/2,1000),'Red: original clot nodes / surviving components    Blue-grey: fixed anchors',font=font(24),fill='#303b44',anchor='mt')
 return im
def movie(pl,bounds,refs,series,names,path,panels):
 encoder=subprocess.Popen(['ffmpeg','-hide_banner','-loglevel','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{WIDTH}x{HEIGHT}','-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','medium','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart',str(path)],stdin=subprocess.PIPE)
 record=[]
 for frame,items in enumerate(series):
  stats=[scene(pl,j,refs[j],state(p),panels==1) for j,p in enumerate(items)];pl.render();im=annotate(pl,pl.screenshot(return_img=True),bounds,stats,names,panels);repeat=24 if frame in [0,len(series)-1] else 2
  for _ in range(repeat):encoder.stdin.write(np.asarray(im).tobytes())
  record.append(dict(frame=frame,repeat=repeat,sources=[str(p.relative_to(R)) for p in items],stats=stats))
  if frame in [0,len(series)//2,len(series)-1]:im.save(F/'visual_review'/f'{path.stem}_{frame:04d}.png')
 encoder.stdin.close();assert encoder.wait()==0;pl.close();return record

def overview(ref,wall):
 pl=pv.Plotter(off_screen=True,window_size=(WIDTH,HEIGHT));pl.set_background('white')
 mesh=wall.wall.copy();mesh.points*=1e3
 pl.add_mesh(mesh,color='#adbdc8',opacity=.14,smooth_shading=True)
 n=int(ref['clot_count']);points=pv.PolyData(ref['X'][:n]*1e3);points['r']=np.cbrt(3*ref['V'][:n]/(4*np.pi))*1e3
 pl.add_mesh(points.glyph(geom=pv.Sphere(radius=1,theta_resolution=10,phi_resolution=10),scale='r',orient=False),color='#a92140',smooth_shading=True)
 center=np.array(mesh.center);span=max(np.ptp(mesh.points,axis=0));pl.camera_position=[center+span*np.array([.2,-1.5,.75]),center,(0,0,1)];pl.enable_parallel_projection();pl.camera.parallel_scale=span*.34
 im=Image.fromarray(pl.screenshot(return_img=True)).convert('RGB');draw=ImageDraw.Draw(im)
 draw.text((WIDTH/2,26),'Original engineered Y vessel + V4 clot',anchor='mt',font=ImageFont.truetype(FONT,32),fill='#283641')
 pixels=5*HEIGHT/(2*pl.camera.parallel_scale);left=100;bottom=HEIGHT-105
 draw.line([(left,bottom),(left+pixels,bottom)],fill='#334956',width=4)
 draw.text((left+pixels/2,bottom+12),'5 mm',anchor='mt',font=ImageFont.truetype(FONT,24),fill='#334956')
 draw.text((WIDTH-100,HEIGHT-48),'Actual geometry ×1; no visual enlargement of the clot',anchor='rt',font=ImageFont.truetype(FONT,24),fill='#465560')
 im.save(F/'figures/FIG_01_Y_VESSEL_LOCATION.png');pl.close()

def force_picture(ref,s,bounds,path):
 from current_surface import CurrentSurface
 surface=CurrentSurface(ref,'uniform');z=surface.evaluate(s['x'],s['velocity'],s['d'],float(s['time']))
 err=float(np.max(abs(z['force']+z['wall_force']-s['external_force'])))
 assert err<1e-15,('saved_force_replay',err)
 pl=setup(bounds,1);st=scene(pl,0,ref,s,False,opacity=.3)
 # All actual current loaded patches, no invented smoothing/continuous shell.
 ids=np.flatnonzero(z['active']&(np.linalg.norm(z['face_force'],axis=1)>0))
 mag=np.linalg.norm(z['traction'][ids],axis=1)
 if len(ids):
  pl.add_mesh(pv.PolyData(z['centroid'][ids]*1e3),scalars=mag,cmap='YlOrRd',clim=(0,float(mag.max())),point_size=8,render_points_as_spheres=True,show_scalar_bar=False,lighting=False,name='traction')
  # Deterministic subset for legibility only; full field saved separately.
  chosen=ids[np.unique(np.linspace(0,len(ids)-1,min(55,len(ids))).astype(int))]
  f=z['face_force'][chosen];scale=.28/max(np.linalg.norm(z['face_force'][ids],axis=1))
  pl.add_arrows(z['centroid'][chosen]*1e3,f,mag=scale,color='#273c4d',show_scalar_bar=False)
 pl.render();im=annotate(pl,pl.screenshot(return_img=True),bounds,[st],['Actual current-surface load'],1);draw=ImageDraw.Draw(im)
 draw.rectangle((0,920,WIDTH,1080),fill='white')
 font=lambda x:ImageFont.truetype(FONT,x)
 import matplotlib
 gradient=(matplotlib.colormaps['YlOrRd'](np.linspace(0,1,500))[:,:3]*255).astype(np.uint8)
 im.paste(Image.fromarray(np.repeat(gradient[None,:,:],18,axis=0)),(320,948))
 draw.text((300,945),'0',anchor='rt',font=font(24),fill='#283b48')
 draw.text((840,945),f'{mag.max(initial=0):.1f} Pa',anchor='lt',font=font(24),fill='#283b48')
 draw.text((570,978),'Current surface traction',anchor='mt',font=font(24),fill='#283b48')
 maxforce=float(np.linalg.norm(z['face_force'],axis=1).max(initial=0))
 draw.text((1370,945),'Arrows: actual face forces',anchor='mt',font=font(25),fill='#283b48')
 draw.text((1370,980),f'0.28 mm arrow = {maxforce*1e6:.3g} µN',anchor='mt',font=font(23),fill='#485963')
 draw.text((WIDTH/2,1030),'Red: clot nodes   Blue-grey: fixed anchors   Translucent markers reveal the load field',anchor='mt',font=font(23),fill='#485963')
 im.save(path);pl.close()
 np.savez_compressed(O/'ACTUAL_LOAD_FIELD.npz',time=s['time'],x=s['x'],velocity=s['velocity'],d=s['d'],centroid=z['centroid'],normal=z['normal'],area=z['area'],weight=z['weight'],active=z['active'],traction=z['traction'],face_force=z['face_force'],nodal_force=z['force'],wall_force=z['wall_force'],shown_arrow_ids=chosen if len(ids) else np.empty(0,int))
 dump(O/'LOAD_VISUALIZATION_AUDIT.json',dict(status='PASS',time_s=float(s['time']),force_replay_max_N=err,active_faces=int(len(ids)),shown_force_arrows=min(55,len(ids)),all_faces_retained_in_npz=True,force_times_area_error_N=float(np.max(abs(z['face_force']-z['traction']*z['area'][:,None]))),source_time_is_actual=True,arrow_scale_mm_per_N=float(scale) if len(ids) else None,maximum_actual_traction_Pa=float(mag.max(initial=0))))
 return surface

def main(preview=False):
 from current_surface import CurrentSurface
 (F/'videos').mkdir(parents=True,exist_ok=True);(F/'visual_review').mkdir(exist_ok=True);(F/'figures').mkdir(exist_ok=True)
 ref=load_npz(O/'provenance/reference.npz');paths=state_paths(CASE)
 if preview:paths=paths[:1]
 else:
  result=json.loads((O/'NATIVE_V4_RESULT.json').read_text());assert result['status'] in ['COMPLETED','NUMERICAL_FAILURE'];paths=[q for q in paths if q[0]<=result['accepted_time_s']+1e-14]
 bounds=np.array([[np.inf,-np.inf]]*3)
 for _,p in paths:
  x=state(p)['x']*1e3;bounds[:,0]=np.minimum(bounds[:,0],x.min(0));bounds[:,1]=np.maximum(bounds[:,1],x.max(0))
 span=bounds[:,1]-bounds[:,0];bounds[:,0]-=.06*span+.03;bounds[:,1]+=.06*span+.03;bounds=bounds.ravel()
 pl=setup(bounds,1)
 for tag,(_,p) in [('INITIAL',paths[0]),('FINAL',paths[-1])]:
  st=scene(pl,0,ref,state(p),True);pl.render();im=annotate(pl,pl.screenshot(return_img=True),bounds,[st],['Y clot — V4 irregular LPS'],1);im.save(F/'figures'/f'FIG_02_{tag}_CLOT.png')
 pl.close()
 if preview:
  overview(ref,CurrentSurface(ref,'uniform').wall);print('RESULT_RENDER_PREVIEW_COMPLETE',flush=True);return
 from analyze_result import table
 history=table(CASE/'history.csv')
 large=max(history,key=lambda q:int(q['largest_nodes']))
 by_time={round(float(q['time_s']),12):q for q in history}
 fragment_state=max(paths,key=lambda q:int(by_time.get(round(q[0],12),{}).get('largest_nodes',0)))
 pl=setup(bounds,1);st=scene(pl,0,ref,state(fragment_state[1]),True);pl.render()
 annotate(pl,pl.screenshot(return_img=True),bounds,[st],['Y clot — saved large-fragment state'],1).save(F/'figures/FIG_02_LARGE_FRAGMENT_CLOT.png');pl.close()
 # Use saved state nearest to the maximum observed actual traction; preserve its time.
 peak=max(history,key=lambda q:float(q['actual_peak_slip_shear_Pa']));selected=min(paths,key=lambda q:abs(q[0]-float(peak['time_s'])))
 s=load_npz(selected[1]);surface=force_picture(ref,s,bounds,F/'figures/FIG_03_ACTUAL_SURFACE_LOAD.png');overview(ref,surface.wall)
 logs=movie(setup(bounds,1),bounds,[ref],[[p] for _,p in paths],['Y clot — V4 irregular LPS'],F/'videos/VIDEO_01_Y_CLOT_V4_FRAGMENTATION.mp4',1)
 dump(O/'DYNAMICS_VIDEO_PROVENANCE.json',dict(status='COMPLETE',width=WIDTH,height=HEIGHT,fps=FPS,world_Z_vertical=True,displacement_scale=1,actual_saved_states_without_interpolation=True,frames=logs,geometry='Original engineered wall-attached Y clot; no patient-specific claim',load_picture_source=str(selected[1].relative_to(R)),history_peak_time_s=float(peak['time_s']),load_picture_actual_time_s=selected[0],large_fragment_picture_source=str(fragment_state[1].relative_to(R)),large_fragment_picture_actual_time_s=fragment_state[0],history_largest_free_component_time_s=float(large['time_s']),bounds_mm=bounds,render_source_sha256=sha(A/'render_results.py')))
 print('RESULT_RENDER_COMPLETE',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--preview',action='store_true');main(p.parse_args().preview)
