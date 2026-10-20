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
 pl=pv.Plotter(off_screen=True,window_size=(WIDTH,HEIGHT),shape=(1,panels));lo=np.array(bounds)[::2];hi=np.array(bounds)[1::2];center=(lo+hi)/2;span=max(hi-lo);direction=np.array([1.2,-3.,.7]);corners=np.array([(x,y,z) for x in bounds[:2] for y in bounds[2:4] for z in bounds[4:]])
 for j in range(panels):
  pl.subplot(0,j);pl.set_background('white');pl.renderer.SetViewport(j/panels+.02,.14,(j+1)/panels-.02,.9);pl.camera_position=[center+span*direction,center,(0,0,1)];pl.enable_parallel_projection();forward=-direction/np.linalg.norm(direction);right=np.cross(forward,[0,0,1]);right/=np.linalg.norm(right);up=np.cross(right,forward);aspect=(WIDTH*(1/panels-.04))/(HEIGHT*.76);pl.camera.parallel_scale=float(1.32*max(np.ptp(corners@up)/2,np.ptp(corners@right)/(2*aspect)));pl.camera.clipping_range=(span*.01,span*20)
  g=grid_points(bounds);pl.add_mesh(line_mesh(g,np.arange(len(g)).reshape(-1,2)),color='#76848f',opacity=.5,line_width=.8,lighting=False,name=f'grid{j}')
 return pl
def scene(pl,j,ref,s,draw_edges):
 pl.subplot(0,j);n=int(ref['clot_count']);N=len(ref['X']);lab,attached=topology(ref,s['d']);rgb=np.tile([171,27,49],(N,1)).astype(np.uint8);rgb[n:]=[124,145,160];counts=np.bincount(lab);free=~attached[lab];order=np.argsort(lab,kind='stable');starts=np.r_[0,np.cumsum(counts)]
 for label in np.flatnonzero(~attached):
  ids=np.sort(order[starts[label]:starts[label+1]]);code=int.from_bytes(hashlib.sha256((ids+1).astype('<i8').tobytes()).digest()[:4],'little');rgb[ids]=PALETTE[code%len(PALETTE)]
 points=pv.PolyData(s['x']*1e3);points['radius']=np.cbrt(3*ref['V']/(4*np.pi))*1e3;points['rgb']=rgb
 glyph=points.glyph(geom=pv.Sphere(radius=1.,theta_resolution=10,phi_resolution=10),scale='radius',orient=False)
 pl.add_mesh(glyph,scalars='rgb',rgb=True,name=f'nodes{j}',show_scalar_bar=False,smooth_shading=True,ambient=.25,diffuse=.75,specular=.12)
 b=ref['bonds'];target=free& (counts[lab]>=20)&(counts[lab]<=200);alive=(s['d']<1)&target[b].all(1)
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
  for axis in range(3):
   anchor=lo.copy()
   if axis==1:anchor[0]=hi[0]
   if axis==2:anchor[:2]=hi[:2]
   for value in np.linspace(lo[axis],hi[axis],4):
    if axis==1 and value in (lo[axis],hi[axis]):continue # avoid shared X/Y and Y/Z corner labels
    p=anchor.copy();p[axis]=value;xy=project(p)+np.array([0,14] if axis==0 else [30,-12] if axis==1 else [18,0]);draw.text(tuple(xy),f'{value:.2f}',font=font(20),fill='#384651',anchor='mt' if axis==0 else 'lt' if axis==1 else 'lm')
   p=anchor.copy();p[axis]=(lo[axis]+hi[axis])/2;xy=project(p)+np.array([0,47] if axis==0 else [47,45] if axis==1 else [57,0]);draw.text(tuple(xy),'XYZ'[axis]+' (mm)',font=font(22),fill='#263540',anchor='mt' if axis==0 else 'lt' if axis==1 else 'lm')
  draw.text((WIDTH*(j+.5)/panels,940),f'20–200 nodes: {stats[j]["N20_200"]} free components',font=font(25),fill='#642136',anchor='mt')
 draw.text((WIDTH/2,1000),'Red: original clot nodes / surviving components    Blue-grey: fixed anchors',font=font(24),fill='#303b44',anchor='mt')
 draw.text((WIDTH/2,1040),'Volume-equivalent nodal markers; not a continuous surface. Displacement ×1.',font=font(20),fill='#59636c',anchor='mt')
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
def main(preview=False):
 (F/'videos').mkdir(parents=True,exist_ok=True);(F/'visual_review').mkdir(exist_ok=True);refs=[load_npz(V3_O/'provenance/reference.npz'),load_npz(O/'provenance/reference.npz')];p0=state_paths(OLD);p1=state_paths(CASE)
 if preview:
  pairs=[[p0[0][1],p1[0][1]]];times=[0.]
 else:
  result=json.loads((O/'NATIVE_V4_RESULT.json').read_text());assert result['status'] in ['COMPLETED','NUMERICAL_FAILURE'];end=result['accepted_time_s'];times=[t for t,_ in p1 if t<=end+1e-14];pairs=[[choose(p0,t)[1],choose(p1,t)[1]] for t in times]
 bounds=np.array([[np.inf,-np.inf]]*3)
 for j in range(2):
  for p in dict.fromkeys(q[j] for q in pairs):
   x=state(p)['x']*1e3;bounds[:,0]=np.minimum(bounds[:,0],x.min(0));bounds[:,1]=np.maximum(bounds[:,1],x.max(0))
 span=bounds[:,1]-bounds[:,0];bounds[:,0]-=.07*span+.04;bounds[:,1]+=.07*span+.04;bounds=bounds.ravel();pl=setup(bounds,2)
 if preview:
  stats=[scene(pl,j,refs[j],state(pairs[0][j]),False) for j in range(2)];pl.render();im=annotate(pl,pl.screenshot(return_img=True),bounds,stats,['V3 regular','V4 irregular'],2);im.save(F/'visual_review/REFERENCE_VIDEO_PREVIEW.png');pl.close();print('VIDEO_PREVIEW_READY',flush=True);return
 logs={};logs['VIDEO_01_REGULAR_VS_IRREGULAR_FRACTURE']=movie(pl,bounds,refs,pairs,['V3 regular','V4 irregular'],F/'videos/VIDEO_01_REGULAR_VS_IRREGULAR_FRACTURE.mp4',2)
 logs['VIDEO_02_IRREGULAR_FRAGMENT_EVOLUTION']=movie(setup(bounds,1),bounds,[refs[1]],[[p[1]] for p in pairs],['V4 irregular — surviving original material graph'],F/'videos/VIDEO_02_IRREGULAR_FRAGMENT_EVOLUTION.mp4',1)
 metadata=dict(status='COMPLETE',width=WIDTH,height=HEIGHT,fps=FPS,bounds_mm=bounds,camera='Fixed orthographic, world Z vertical; same world bounds',displacement_scale=1.,time_matching='Use actual saved V4 times and latest V3 accepted state at or before each; show both actual times. Repeat frames, no interpolation.',color='Deterministic SHA256 of exact original component member IDs selects a red palette; attached nodes red, anchors blue-grey',volume_marker='Sphere volume equals Vi; visualization quadrature marker, not a true spherical fiber or continuous clot surface',shown_bonds='Video 2 shows ALL internal surviving original bonds of 20–200-node free components; no reference failed bonds',frames=logs)
 metadata['render_source_sha256']=sha(A/'render_videos.py');metadata['analysis_source_sha256']=sha(A/'analyze_result.py');dump(O/'VIDEO_PROVENANCE.json',metadata);print('V4_VIDEOS_COMPLETE',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--preview',action='store_true');main(p.parse_args().preview)
