"""English blue-particle animation; inherited BraVa style and exact 1.30 zoom.

Playback-only slowdown uses the saved accepted trajectory endpoints, never a
new particle solve or interpolation across birth/terminal events.
"""
from pathlib import Path
import argparse,json,csv,os,hashlib
os.environ.setdefault('VTK_DEFAULT_OPENGL_WINDOW','vtkEGLRenderWindow')
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from style_reference import projected_bounds,enlarge_animation
ap=argparse.ArgumentParser();ap.add_argument('--run-name',required=True);ap.add_argument('--pose',default='0');ap.add_argument('--preview',action='store_true');ap.add_argument('--playback-speed',type=float,default=.2);args=ap.parse_args()
if not 0 < args.playback_speed <= 1:ap.error('--playback-speed must be in (0, 1]')
root=Path(__file__).resolve().parents[1];base=root.parent;analysis=root/'analysis'/args.run_name
out=root/('candidate_'+args.pose);out.mkdir(exist_ok=True)
for folder in ['figures','animations','data']:(out/folder).mkdir(exist_ok=True)
rows=list(csv.DictReader((analysis/'trajectory_catalog.csv').open()));summary=json.loads((root/'runs'/args.run_name/'SUMMARY.json').read_text())
identity=json.loads((root/'runs'/args.run_name/'IDENTITY.json').read_text())
dt_s=identity['population_contract']['nominal_dt_s']
inlet_flow_m3_s=identity['population_contract']['source']['Q_in_m3_s']
data=np.load(analysis/'visualization_data.npz');geometry=json.loads((base/'inputs/geometry.json').read_text())
T=np.array(geometry['pose_transforms'][args.pose]['source_mm_to_candidate_mm']);origin=np.array(geometry['origin_source_mm'])
def pose(x):return (np.asarray(x)*1e3+origin)@T[:3,:3].T+T[:3,3]
config=json.loads((root/'data/environment.json').read_text())
if not Path(config['flow']).exists():config=json.loads((root/'data/environment.local.json').read_text())
# Authoritative physical surface, no geometry substitution.
wallpath=base/'mesh/SV_MESH/mesh-surfaces/WALL.vtp'
surface=pv.read(config['flow']).extract_surface(algorithm='dataset_surface');surface.points=pose(surface.points)
center=np.array(surface.center);radius=np.linalg.norm(np.ptp(surface.points,axis=0))/2
direction=np.array([-.72,.69,.11]);direction/=np.linalg.norm(direction)
BACKGROUND='#101824';BUBBLE_COLOR='#3399ff'
COLORS={k:BUBBLE_COLOR for k in ['OUTLET_01','OUTLET_02','OUTLET_03','SUPPORTED_STATIONARY','OBSERVATION_CENSORED','SOLVER_FAILURE']}
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
inlet_area_m2=float(pv.read(base/'mesh/SV_MESH/mesh-surfaces/INLET.vtp').area)
inlet_mean_speed_m_s=inlet_flow_m3_s/inlet_area_m2
ports={}
import vtk
for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
 p=pv.read(base/'mesh/SV_MESH/mesh-surfaces'/(role+'.vtp'))
 xyz=pose(np.array(p.center));c=vtk.vtkCoordinate();c.SetCoordinateSystemToWorld();c.SetValue(*xyz)
 xy=c.GetComputedDoubleDisplayValue(plot.renderer);ports[role]=[float(xy[0]),float(1080-xy[1])]
def overlay(pixels,t=None,static=False):
 img=Image.fromarray(pixels);d=ImageDraw.Draw(img)
 title='Continuous microbubble perfusion'
 time_text='Complete recorded trajectories' if static else f'Physical time: {t:.3f} s'
 header=[title,time_text,f'dt = {dt_s*1000:g} ms',
         f'Inlet mean speed: {inlet_mean_speed_m_s:.3f} m/s ({inlet_flow_m3_s*60e6:g} mL/min)']
 for j,line in enumerate(header):
  assert titlefont.getlength(line)<=1840
  d.text((40,22+40*j),line,font=titlefont,fill='white')
 assert d.textbbox((40,142),header[-1],font=titlefont)[3]<layout['enlarged_bounds_px'][0][1]
 # Same neutral side labels, rigidly transformed official cap identity.
 for role,(x,y) in (ports.items() if static else []):
  name='INLET' if role=='INLET' else 'O'+str(int(role[-2:]));left=x<960;tx=55 if left else 1780
  # At far frame edges put labels in a vertical margin, not on the cap.
  label_offset=-50 if x<180 or x>1700 else 0
  ty=float(np.clip(y+label_offset,135,965))
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
original_ft=data['frame_time_s'];original_fx=data['frame_x_m'];original_frames=len(original_ft);fps=24
frame_count=int(round(original_frames/args.playback_speed))
ft=np.linspace(original_ft[0],original_ft[-1],frame_count)
fx_m=np.full((frame_count,len(rows),3),np.nan)
reference_error_m=0.;track_hashes={}
assert np.array_equal(data['particle_id'],np.array([int(r['particle_id']) for r in rows]))
for i,row in enumerate(rows):
 track_path=analysis/'tracks'/f"mb_{int(row['particle_id']):06d}.npz"
 with np.load(track_path) as track:
  ts=track['physical_time_s'];xs=track['x_m']
  assert np.all(np.diff(ts)>0)
  assert abs(ts[0]-birth[i])<1e-12 and abs(ts[-1]-end[i])<1e-10
  live=(ft>=ts[0])&(ft<ts[-1]);old_live=(original_ft>=ts[0])&(original_ft<ts[-1])
  assert np.array_equal(np.isfinite(original_fx[:,i,0]),old_live)
  for component in range(3):
   fx_m[live,i,component]=np.interp(ft[live],ts,xs[:,component])
   old_x=np.interp(original_ft[old_live],ts,xs[:,component])
   if len(old_x):reference_error_m=max(reference_error_m,float(np.max(abs(old_x-original_fx[old_live,i,component]))))
 with track_path.open('rb') as f:track_hashes[track_path.name]=hashlib.file_digest(f,'sha256').hexdigest()
 if i%250==0:print('RESAMPLED_TRACK',i,len(rows),flush=True)
assert reference_error_m<1e-12,reference_error_m
(out/'data/animation_track_hashes.json').write_text(json.dumps(track_hashes,indent=2)+'\n')
np.savez_compressed(out/'data/animation_samples.npz',frame_time_s=ft,frame_x_m=fx_m,particle_id=data['particle_id'])
fx=pose(fx_m);snapshot_index=int(np.argmin(abs(ft-summary['infusion_stop_time_s'])));indices=[snapshot_index] if args.preview else range(len(ft))
if not args.preview:
 import imageio.v2 as imageio
 os.environ['IMAGEIO_FFMPEG_EXE']='/usr/bin/ffmpeg'
 writer=imageio.get_writer(out/'animations/microbubble_physical_time.mp4',fps=fps,codec='libx264',quality=8,macro_block_size=1,ffmpeg_params=['-threads','2'])
for j in indices:
 for key,ids in groups.items():clouds[key].points=fx[j,ids]
 plot.render();img=overlay(plot.screenshot(return_img=True),float(ft[j]))
 if j==snapshot_index:img.save(out/'figures/02_physical_time_snapshot.png',dpi=(300,300))
 if not args.preview:writer.append_data(np.asarray(img))
 if j%48==0:print('FRAME',args.pose,j,len(ft),flush=True)
if not args.preview:writer.close()
plot.close()
record=dict(run_name=args.run_name,flow_sha256=identity['flow_sha256'],
 production_git_commit=identity['git_commit'],population_contract_sha256=identity['population_contract_sha256'],
 source_identity_sha256=identity['source_identity_sha256'],source_distribution_identity=identity['source_distribution_identity'],
 master_seed=identity['master_seed'],realization_id=identity['realization_id'],
 renderer_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 style_script_sha256=hashlib.sha256((Path(__file__).parent/'style_reference.py').read_bytes()).hexdigest(),
 pyvista_version=pv.__version__,vtk_version=vtk.vtkVersion.GetVTKVersion(),
 candidate_id=int(args.pose),text_language='English',video_semantics='PHYSICAL_TIME_SHARED_POPULATION_NOT_AGE_ALIGNED',
 original_camera_style='microbubble/scripts/render_results.py',zoom=1.30,layout=layout,marker_size_display_px=9.1,
 render_backend=gpu_strings,encoder='libx264_CPU_2_THREADS',physical_time_end_s=float(ft[-1]),frames=len(ft),fps=fps,
 playback_speed_relative_to_previous=args.playback_speed,previous_frames=original_frames,previous_duration_s=original_frames/fps,
 video_duration_s=len(ft)/fps,physical_time_per_video_second=float(ft[-1]-ft[0])/((len(ft)-1)/fps),
 original_animation_sample_max_error_m=reference_error_m,uniform_particle_color=BUBBLE_COLOR,
 header_layout='SINGLE_COLUMN_TITLE_PHYSICAL_TIME_DT_INLET_MEAN_SPEED',footer=None,
 inlet_mean_speed_m_s=inlet_mean_speed_m_s,inlet_area_m2=inlet_area_m2,inlet_flow_m3_s=inlet_flow_m3_s,nominal_dt_s=dt_s,
 all_particles_included=True,position_interpolation='Exact held-velocity path from every accepted endpoint',
 static_lines='display decimation <=~1500 exact vertices per track; every accepted state retained in raw files',
 outcome_color='DISABLED: all particles blue; original terminal categories retained in trajectory_catalog.csv',ports_px=ports,preview_only=args.preview)
(out/'data/render_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
if not args.preview:
 (out/'OPEN_RESULTS.html').write_text(f'''<!doctype html><meta charset="utf-8"><title>BraVa continuous microbubble perfusion</title><style>body{{background:#101824;color:#dce2ed;font:17px sans-serif;max-width:1400px;margin:30px auto}}video,img{{width:100%}}a{{color:#39d8e6}}</style><h1>BraVa continuous microbubble perfusion</h1><p>Candidate {args.pose}. Playback: {args.playback_speed:g}x the previous video speed; physical timestamps refer to the original recorded trajectories. All particles are blue. Original camera, vessel transparency, marker size and 1.30 zoom are retained.</p><video src="animations/microbubble_physical_time.mp4?v=blue-slow-20260929" controls loop></video><p><a href="figures/01_trajectories.pdf">Trajectory PDF</a> · <a href="data/render_manifest.json">Rendering provenance</a> · <a href="../POPULATION_VALIDATION_REPORT_ZH.md">Original simulation report and terminal outcomes</a></p><img src="figures/01_trajectories.png?v=blue-slow-20260929">''')
print(json.dumps(record),flush=True)
