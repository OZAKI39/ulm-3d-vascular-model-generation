"""English physical-time animation; inherited BraVa style and exact 1.30 zoom."""
from pathlib import Path
import argparse,json,csv,os
os.environ.setdefault('VTK_DEFAULT_OPENGL_WINDOW','vtkEGLRenderWindow')
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from style_reference import projected_bounds,enlarge_animation
ap=argparse.ArgumentParser();ap.add_argument('--run-name',required=True);ap.add_argument('--pose',default='0');ap.add_argument('--preview',action='store_true');args=ap.parse_args()
root=Path(__file__).resolve().parents[1];base=root.parent;analysis=root/'analysis'/args.run_name
out=root/('candidate_'+args.pose);out.mkdir(exist_ok=True)
for folder in ['figures','animations','data']:(out/folder).mkdir(exist_ok=True)
rows=list(csv.DictReader((analysis/'trajectory_catalog.csv').open()));summary=json.loads((root/'runs'/args.run_name/'SUMMARY.json').read_text())
data=np.load(analysis/'visualization_data.npz');geometry=json.loads((base/'inputs/geometry.json').read_text())
T=np.array(geometry['pose_transforms'][args.pose]['source_mm_to_candidate_mm']);origin=np.array(geometry['origin_source_mm'])
def pose(x):return (np.asarray(x)*1e3+origin)@T[:3,:3].T+T[:3,3]
config=json.loads((root/'data/environment.json').read_text())
# Authoritative physical surface, no geometry substitution.
wallpath=base/'mesh/SV_MESH/mesh-surfaces/WALL.vtp'
surface=pv.read(config['flow']).extract_surface(algorithm='dataset_surface');surface.points=pose(surface.points)
center=np.array(surface.center);radius=np.linalg.norm(np.ptp(surface.points,axis=0))/2
direction=np.array([-.72,.69,.11]);direction/=np.linalg.norm(direction)
BACKGROUND='#101824';COLORS={'OUTLET_01':'#b994ff','OUTLET_02':'#39d8e6','OUTLET_03':'#ffb153','SUPPORTED_STATIONARY':'#bac2cf','OBSERVATION_CENSORED':'#f672a7','SOLVER_FAILURE':'#ff4b4b'}
labels={'OUTLET_01':'O1','OUTLET_02':'O2','OUTLET_03':'O3','SUPPORTED_STATIONARY':'Supported stop','OBSERVATION_CENSORED':'Censored','SOLVER_FAILURE':'Solver stop'}
keys=[r['outlet'] or r['status'] for r in rows];groups={k:np.array([i for i,q in enumerate(keys) if q==k]) for k in COLORS if k in keys}
font=root/'assets/NotoSansCJKsc-Regular.otf';titlefont=ImageFont.truetype(str(font),30);smallfont=ImageFont.truetype(str(font),23)
plot=pv.Plotter(off_screen=True,window_size=(1920,1080));plot.set_background(BACKGROUND)
plot.add_mesh(surface,color='#b2becf',opacity=.16,smooth_shading=True)
plot.camera.position=center+direction*radius*4;plot.camera.focal_point=center;plot.camera.up=(0,0,1)
plot.camera.parallel_projection=True;plot.camera.parallel_scale=radius*.78
plot.show(auto_close=False);bounds=projected_bounds(plot,surface.points);span=bounds[1]-bounds[0]
plot.camera.parallel_scale*=max(span[0]/1750,span[1]/840)*1.30
layout=enlarge_animation(plot,surface)
caps=plot.ren_win.ReportCapabilities();gpu_strings=[s for s in caps.splitlines() if any(k in s for k in ['vendor string','renderer string','version string'])]
assert any('NVIDIA' in s for s in gpu_strings),gpu_strings
birth=np.array([float(r['birth_time_s']) for r in rows]);end=np.array([float(r['physical_time_s']) for r in rows])
ports={}
import vtk
for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
 p=pv.read(base/'mesh/SV_MESH/mesh-surfaces'/(role+'.vtp'))
 xyz=pose(np.array(p.center));c=vtk.vtkCoordinate();c.SetCoordinateSystemToWorld();c.SetValue(*xyz)
 xy=c.GetComputedDoubleDisplayValue(plot.renderer);ports[role]=[float(xy[0]),float(1080-xy[1])]
def overlay(pixels,t=None,static=False):
 img=Image.fromarray(pixels);d=ImageDraw.Draw(img)
 if static:title='Continuous microbubble perfusion | Complete recorded trajectories'
 else:
  active=int(((birth<=t)&(end>t)).sum());mode='ON' if t<summary['infusion_stop_time_s'] else 'OFF'
  title=f'Continuous microbubble perfusion | Physical time: {t:.3f} s | Infusion: {mode} | Active: {active}'
 assert titlefont.getlength(title)<=1840;d.text((40,22),title,font=titlefont,fill='white')
 for j,(key,ids) in enumerate(groups.items()):
  text=f'{labels[key]}: {len(ids)}';d.text((40+j*min(350,1840/len(groups)),69),text,font=smallfont,fill=COLORS[key])
 footer='Color = final outcome | dt = 0.5 ms | C = 4e9 m^-3 | 18 mL/min | Markers enlarged for display; physical sizes are unchanged.'
 assert smallfont.getlength(footer)<=1840;d.text((40,1018),footer,font=smallfont,fill='#cad4e2')
 # Same neutral side labels, rigidly transformed official cap identity.
 for role,(x,y) in (ports.items() if static else []):
  name='INLET' if role=='INLET' else 'O'+str(int(role[-2:]));left=x<960;tx=55 if left else 1780;ty=float(np.clip(y,135,965))
  d.text((tx,ty-14),name,font=smallfont,fill='white')
  anchor=(tx+85 if left else tx-8,ty+2);d.line([anchor,(x,y)],fill='#dce2ed',width=2)
  vec=np.array(anchor)-[x,y];vec/=max(np.linalg.norm(vec),1e-30);side=np.array([-vec[1],vec[0]])
  d.polygon([tuple([x,y]),tuple(np.array([x,y])+12*vec+4*side),tuple(np.array([x,y])+12*vec-4*side)],fill='#dce2ed')
 return img
# All tracks, with explicit display decimation. Raw trajectories are unaffected.
points=pose(data['track_points_m']);offset=data['track_offsets'];actors=[]
for key,ids in groups.items():
 paths=[points[offset[i]:offset[i+1]] for i in ids];xyz=np.concatenate(paths);lines=[];start=0
 for path in paths:lines.extend([len(path),*range(start,start+len(path))]);start+=len(path)
 poly=pv.PolyData(xyz);poly.lines=np.array(lines);actors.append(plot.add_mesh(poly,color=COLORS[key],line_width=1.35,opacity=.68,lighting=False))
plot.render();img=overlay(plot.screenshot(return_img=True),static=True)
img.save(out/'figures/01_trajectories.png',dpi=(300,300));img.convert('RGB').save(out/'figures/01_trajectories.pdf',resolution=300)
for actor in actors:plot.remove_actor(actor)
clouds={};actors={}
for key,ids in groups.items():
 cloud=pv.PolyData(np.zeros((len(ids),3)));actor=plot.add_mesh(cloud,color=COLORS[key],point_size=9.1,render_points_as_spheres=True,lighting=True)
 clouds[key]=cloud;actors[key]=actor
ft=data['frame_time_s'];fx=pose(data['frame_x_m']);snapshot_index=int(np.argmin(abs(ft-summary['infusion_stop_time_s'])));indices=[snapshot_index] if args.preview else range(len(ft))
if not args.preview:
 import imageio.v2 as imageio
 os.environ['IMAGEIO_FFMPEG_EXE']='/usr/bin/ffmpeg'
 writer=imageio.get_writer(out/'animations/microbubble_physical_time.mp4',fps=24,codec='libx264',quality=8,macro_block_size=1,ffmpeg_params=['-threads','2'])
for j in indices:
 for key,ids in groups.items():clouds[key].points=fx[j,ids]
 plot.render();img=overlay(plot.screenshot(return_img=True),float(ft[j]))
 if j==snapshot_index:img.save(out/'figures/02_physical_time_snapshot.png',dpi=(300,300))
 if not args.preview:writer.append_data(np.asarray(img))
 if j%48==0:print('FRAME',args.pose,j,len(ft),flush=True)
if not args.preview:writer.close()
plot.close()
record=dict(run_name=args.run_name,flow_sha256=json.loads((root/'runs'/args.run_name/'IDENTITY.json').read_text())['flow_sha256'],
 candidate_id=int(args.pose),text_language='English',video_semantics='PHYSICAL_TIME_SHARED_POPULATION_NOT_AGE_ALIGNED',
 original_camera_style='microbubble/scripts/render_results.py',zoom=1.30,layout=layout,marker_size_display_px=9.1,
 render_backend=gpu_strings,encoder='libx264_CPU_2_THREADS',physical_time_end_s=float(ft[-1]),frames=len(ft),fps=24,
 all_particles_included=True,position_interpolation='Exact held-velocity path from every accepted endpoint',
 static_lines='display decimation <=~1500 exact vertices per track; every accepted state retained in raw files',
 outcome_color='Future final terminal category, not local instantaneous physics',ports_px=ports,preview_only=args.preview)
(out/'data/render_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
if not args.preview:
 (out/'OPEN_RESULTS.html').write_text('''<!doctype html><meta charset="utf-8"><title>BraVa continuous microbubble perfusion</title><style>body{background:#101824;color:#dce2ed;font:17px sans-serif;max-width:1400px;margin:30px auto}video,img{width:100%}a{color:#39d8e6}</style><h1>BraVa continuous microbubble perfusion</h1><p>Physical time; source C = 4e9 m^-3, Q = 18 mL/min, dt = 0.5 ms. Stop infusion after at least 20 actual center crossings at every outlet, then observe for 12 seconds. Colors indicate final outcomes. Red solver stops are numerical failures, not physiological trapping.</p><video src="animations/microbubble_physical_time.mp4" controls loop></video><p><a href="figures/01_trajectories.pdf">Trajectory PDF</a> · <a href="data/render_manifest.json">Rendering provenance</a> · <a href="../POPULATION_VALIDATION_REPORT_ZH.md">Chinese validation report</a></p><img src="figures/01_trajectories.png">''')
print(json.dumps(record),flush=True)
