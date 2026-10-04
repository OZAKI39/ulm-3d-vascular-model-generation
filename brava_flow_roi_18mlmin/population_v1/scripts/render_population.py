"""English blue-particle animation; inherited BraVa style, maximum safe fit.

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
ap=argparse.ArgumentParser();ap.add_argument('--run-name',required=True);ap.add_argument('--pose',default='0');ap.add_argument('--preview',action='store_true');ap.add_argument('--playback-speed',type=float,default=.2);ap.add_argument('--hide-audited-anomalies',action='store_true');args=ap.parse_args()
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
BACKGROUND='#000000';BUBBLE_COLOR='#3399ff'
COLORS={k:BUBBLE_COLOR for k in ['OUTLET_01','OUTLET_02','OUTLET_03','SUPPORTED_STATIONARY','OBSERVATION_CENSORED','SOLVER_FAILURE']}
labels={'OUTLET_01':'O1','OUTLET_02':'O2','OUTLET_03':'O3','SUPPORTED_STATIONARY':'Supported stop','OBSERVATION_CENSORED':'Censored','SOLVER_FAILURE':'Solver stop'}
excluded=set();failure_ids=set();risk_ids=set();filter_sources={}
if args.hide_audited_anomalies:
 assert args.run_name=='realization_0000_contactfix','Existing audit applies only to this recorded run'
 audit_identity=json.loads((root/'cfd_wall_audit/AUDIT_IDENTITY.json').read_text())
 assert audit_identity['sha256']['flow']==identity['flow_sha256']
 risk_path=root/'cfd_wall_audit/acceleration_overlap.csv'
 risk_rows=list(csv.DictReader(risk_path.open()))
 risk_ids={int(r['particle_id']) for r in risk_rows if r['window_any_eta_gt_0p5']=='True'}
 failure_ids={int(r['particle_id']) for r in rows if r['status']=='SOLVER_FAILURE'}
 excluded=risk_ids|failure_ids
 assert excluded <= {int(r['particle_id']) for r in rows}
 for path in [risk_path,analysis/'trajectory_catalog.csv',root/'cfd_wall_audit/AUDIT_IDENTITY.json']:
  with path.open('rb') as f:filter_sources[str(path.relative_to(root))]=hashlib.file_digest(f,'sha256').hexdigest()
visibility=[]
for row in rows:
 pid=int(row['particle_id']);reasons=[]
 if pid in failure_ids:reasons.append('SOLVER_FAILURE')
 if pid in risk_ids:reasons.append('SPEEDUP_WINDOW_ETA_GT_0P5_NUMERICAL_RISK')
 visibility.append(dict(particle_id=pid,displayed=pid not in excluded,status=row['status'],outlet=row['outlet'],hidden_reason=';'.join(reasons)))
with (out/'data/display_selection.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(visibility[0]));w.writeheader();w.writerows(visibility)
keys=[r['outlet'] or r['status'] for r in rows]
groups={k:np.array([i for i,q in enumerate(keys) if q==k and int(rows[i]['particle_id']) not in excluded],dtype=int) for k in COLORS}
groups={k:v for k,v in groups.items() if len(v)}
displayed_ids=[int(r['particle_id']) for r in rows if int(r['particle_id']) not in excluded]
assert sum(len(v) for v in groups.values())==len(displayed_ids)
print('DISPLAY_SELECTION',len(displayed_ids),'shown;',len(excluded),'hidden',flush=True)
font=root/'assets/NotoSansCJKsc-Regular.otf';titlefont=ImageFont.truetype(str(font),30);smallfont=ImageFont.truetype(str(font),23)
plot=pv.Plotter(off_screen=True,window_size=(1920,1080));plot.set_background(BACKGROUND)
plot.add_mesh(surface,color='#b2becf',opacity=.16,smooth_shading=True)
plot.camera.position=center+direction*radius*4;plot.camera.focal_point=center;plot.camera.up=(0,0,1)
plot.camera.parallel_projection=True;plot.camera.parallel_scale=radius*.78
plot.show(auto_close=False);bounds=projected_bounds(plot,surface.points);span=bounds[1]-bounds[0]
plot.camera.parallel_scale*=max(span[0]/1750,span[1]/840)*1.30
layout=enlarge_animation(plot,surface)
# Reuse the requested human-spatial RIGHT-panel annotation implementation.
from spatial_axes_reference import add_spatial_reference_axes,axis_annotations
previous_layout=layout
reference_bounds=tuple(surface.bounds)
reference_corners=pv.Box(bounds=reference_bounds).points
coordinate_axes,orientation_labels=add_spatial_reference_axes(plot,reference_bounds)
safe_bounds=np.array([[40.,208.],[1880.,1035.]])
forward=np.asarray(plot.camera.direction);forward/=np.linalg.norm(forward)
right=np.cross(forward,np.asarray(plot.camera.up));right/=np.linalg.norm(right)
up=np.cross(right,forward)
def fit_content_bounds():
 coordinate_axes.update()
 bbox=projected_bounds(plot,reference_corners)
 labels=axis_annotations(coordinate_axes)
 for label in labels:
  q=label['bbox_px'];bbox[0]=np.minimum(bbox[0],q[:2]);bbox[1]=np.maximum(bbox[1],q[2:])
 return bbox
# Text offsets are fixed in pixels. Iterate scale and pan until the physical
# frame plus its actual projected labels tightly fit below the header.
for _ in range(30):
 bbox=fit_content_bounds();factor=float(np.max((bbox[1]-bbox[0])/(safe_bounds[1]-safe_bounds[0])))
 plot.camera.parallel_scale*=factor
 bbox=fit_content_bounds();delta=safe_bounds.mean(axis=0)-bbox.mean(axis=0)
 shift=(-delta[0]*right+delta[1]*up)*(2*plot.camera.parallel_scale/plot.window_size[1])
 plot.camera.position=np.asarray(plot.camera.position)+shift
 plot.camera.focal_point=np.asarray(plot.camera.focal_point)+shift
 content_bounds=fit_content_bounds()
 if np.max(np.abs(content_bounds.mean(axis=0)-safe_bounds.mean(axis=0)))<1e-6 and abs(factor-1)<1e-9:break
assert np.all(content_bounds[0]>=safe_bounds[0]-1e-4) and np.all(content_bounds[1]<=safe_bounds[1]+1e-4),content_bounds
extra_zoom=previous_layout['parallel_scale']/plot.camera.parallel_scale
after=projected_bounds(plot,surface.points)
layout=dict(mode='MAXIMUM_FIT_3D_FRAME_GRID_AND_OUTWARD_LABELS',safe_bounds_px=safe_bounds.tolist(),
            previous_enlarged_bounds_px=previous_layout['enlarged_bounds_px'],enlarged_bounds_px=after.tolist(),
            content_bounds_px=content_bounds.tolist(),projected_frame_bounds_px=projected_bounds(plot,reference_corners).tolist(),
            zoom_relative_to_previous=extra_zoom,zoom_relative_to_original=1.30*extra_zoom,
            parallel_scale=plot.camera.parallel_scale,camera_position=list(plot.camera.position),
            camera_focal_point=list(plot.camera.focal_point),view_direction=forward.tolist())
grid_info=dict(type='SPATIAL_RIGHT_PANEL_3D_SIX_FACE_GRID',units='mm',bounds_candidate_mm=list(reference_bounds),
               frame_color='#FF5C5C',frame_line_width_px=4.,frame_fill_opacity=.025,
               grid_color='#BFC7D5',grid_line_width_px=1.,grid_opacity=.28,
               tick_counts=dict(X=4,Y=4,Z=5),tick_font_family='Arial',tick_font_size=14,title_font_size=16,
               decimal_places=2,single_visible_edge_per_axis=True,ticks_outward=True,
               view_up_retained='Z',reference_view_up='Y; style reused without rotating the original vessel view',
               tick_values_mm={name:np.linspace(reference_bounds[2*i],reference_bounds[2*i+1],count).tolist() for i,(name,count) in enumerate([('X',4),('Y',4),('Z',5)])},
               annotations=axis_annotations(coordinate_axes),
               provenance='vascular_printing/vascular_processing/spatial_display_adapter.py -> compact_annotations.py + utils/rodent_vasculature/interactive.py',
               module_sha256=hashlib.sha256((Path(__file__).parent/'spatial_axes_reference.py').read_bytes()).hexdigest(),
               box_meaning='DISPLAY_BOUNDING_BOX_OF_THIS_CFD_GEOMETRY; NOT_THE_SAMPLING_ROI_FROM_THE_REFERENCE')
caps=plot.ren_win.ReportCapabilities();gpu_strings=[s for s in caps.splitlines() if any(k in s for k in ['vendor string','renderer string','version string'])]
assert any('NVIDIA' in s for s in gpu_strings),gpu_strings
birth=np.array([float(r['birth_time_s']) for r in rows]);end=np.array([float(r['physical_time_s']) for r in rows])
inlet_area_m2=float(pv.read(base/'mesh/SV_MESH/mesh-surfaces/INLET.vtp').area)
inlet_mean_speed_m_s=inlet_flow_m3_s/inlet_area_m2
ports={};port_centers_m={};port_surface_sha256={}
import vtk
for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
 port_path=base/'mesh/SV_MESH/mesh-surfaces'/(role+'.vtp')
 p=pv.read(port_path)
 port_centers_m[role]=list(p.center)
 with port_path.open('rb') as f:port_surface_sha256[role]=hashlib.file_digest(f,'sha256').hexdigest()
 xyz=pose(np.array(p.center));c=vtk.vtkCoordinate();c.SetCoordinateSystemToWorld();c.SetValue(*xyz)
 xy=c.GetComputedDoubleDisplayValue(plot.renderer);ports[role]=[float(xy[0]),float(1080-xy[1])]
# Short neutral leaders in nearby whitespace, connected to official port caps.
label_offsets={'INLET':(-136,74),'OUTLET_01':(120,64),
               'OUTLET_02':(-136,-78),'OUTLET_03':(-150,-60)}
port_annotations={}
for role,(x,y) in ports.items():
 name='Inlet' if role=='INLET' else 'Outlet '+str(int(role[-2:]))
 center_px=np.array([x,y])+label_offsets[role]
 bbox=ImageDraw.Draw(Image.new('RGB',(1,1))).textbbox(tuple(center_px),name,font=smallfont,anchor='mm')
 assert bbox[0]>=20 and bbox[2]<=1900 and bbox[1]>=200 and bbox[3]<=1056,(role,bbox)
 anchor=np.clip([x,y],np.array(bbox[:2])-6,np.array(bbox[2:])+6)
 port_annotations[role]=dict(text=name,text_center_px=center_px.tolist(),text_bbox_px=list(bbox),
                             leader_start_px=anchor.tolist(),cap_px=[x,y])
def overlay(pixels,t=None,static=False):
 img=Image.fromarray(pixels);d=ImageDraw.Draw(img)
 title='Continuous microbubble perfusion'
 time_text='Displayed recorded trajectories' if static else f'Physical time: {t:.3f} s'
 header=[title,time_text,f'dt = {dt_s*1000:g} ms',
         f'Inlet mean speed: {inlet_mean_speed_m_s:.3f} m/s ({inlet_flow_m3_s*60e6:g} mL/min)']
 for j,line in enumerate(header):
  assert titlefont.getlength(line)<=1840
  d.text((40,22+40*j),line,font=titlefont,fill='white')
 assert d.textbbox((40,142),header[-1],font=titlefont)[3]<layout['enlarged_bounds_px'][0][1]
 for role,(x,y) in ports.items():
  annotation=port_annotations[role]
  d.text(tuple(annotation['text_center_px']),annotation['text'],font=smallfont,anchor='mm',fill='white')
  anchor=tuple(annotation['leader_start_px']);d.line([anchor,(x,y)],fill='#dce2ed',width=2)
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
coordinate_axes.dispose();orientation_labels.dispose();plot.close()
record=dict(run_name=args.run_name,flow_sha256=identity['flow_sha256'],
 production_git_commit=identity['git_commit'],population_contract_sha256=identity['population_contract_sha256'],
 source_identity_sha256=identity['source_identity_sha256'],source_distribution_identity=identity['source_distribution_identity'],
 master_seed=identity['master_seed'],realization_id=identity['realization_id'],
 renderer_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 style_script_sha256=hashlib.sha256((Path(__file__).parent/'style_reference.py').read_bytes()).hexdigest(),
 pyvista_version=pv.__version__,vtk_version=vtk.vtkVersion.GetVTKVersion(),
 candidate_id=int(args.pose),text_language='English',video_semantics='PHYSICAL_TIME_SHARED_POPULATION_NOT_AGE_ALIGNED',
 original_camera_style='microbubble/scripts/render_results.py',zoom=1.30*extra_zoom,layout=layout,marker_size_display_px=9.1,
 background=BACKGROUND,port_annotations=port_annotations,port_centers_m=port_centers_m,port_surface_sha256=port_surface_sha256,coordinate_grid=grid_info,
 render_backend=gpu_strings,encoder='libx264_CPU_2_THREADS',physical_time_end_s=float(ft[-1]),frames=len(ft),fps=fps,
 playback_speed_relative_to_previous=args.playback_speed,previous_frames=original_frames,previous_duration_s=original_frames/fps,
 video_duration_s=len(ft)/fps,physical_time_per_video_second=float(ft[-1]-ft[0])/((len(ft)-1)/fps),
 original_animation_sample_max_error_m=reference_error_m,uniform_particle_color=BUBBLE_COLOR,
 header_layout='SINGLE_COLUMN_TITLE_PHYSICAL_TIME_DT_INLET_MEAN_SPEED',footer=None,
 inlet_mean_speed_m_s=inlet_mean_speed_m_s,inlet_area_m2=inlet_area_m2,inlet_flow_m3_s=inlet_flow_m3_s,nominal_dt_s=dt_s,
 all_particles_included=not excluded,source_particle_count=len(rows),displayed_particle_count=len(displayed_ids),hidden_particle_count=len(excluded),
 displayed_particle_ids=displayed_ids,hidden_particle_ids=sorted(excluded),display_filter_enabled=args.hide_audited_anomalies,
 display_filter=dict(failure_ids=sorted(failure_ids),numerical_risk_ids=sorted(risk_ids),sources_sha256=filter_sources,
                     scope='ENTIRE_TRACK_HIDDEN_IN_ANIMATION_AND_STATIC_LINES; RAW_DATA_UNCHANGED',
                     limitation='eta > 0.5 is a numerical-risk screen, not proof of abnormal physics; retained tracks are not independently certified normal'),
 position_interpolation='Exact held-velocity path from every accepted endpoint',
 static_lines='display decimation <=~1500 exact vertices per track; every accepted state retained in raw files',
 outcome_color='DISABLED: all particles blue; original terminal categories retained in trajectory_catalog.csv',ports_px=ports,preview_only=args.preview)
(out/'data/render_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
if not args.preview:
 (out/'OPEN_RESULTS.html').write_text(f'''<!doctype html><meta charset="utf-8"><title>BraVa continuous microbubble perfusion</title><style>body{{background:#000000;color:#dce2ed;font:17px sans-serif;max-width:1400px;margin:30px auto}}video,img{{width:100%}}a{{color:#39d8e6}}</style><h1>BraVa continuous microbubble perfusion</h1><p>Candidate {args.pose}. Blue particles; black background; labeled official ports. Vessel fitted to the available frame, with the original viewing direction and transparency. XYZ direction axes; physical 3-D coordinate grid and outward mm ticks in the referenced spatial ROI style. Duration: {len(ft)/fps:g} seconds at 24 fps; physical timestamps refer to the original recorded trajectories.</p><p>Filtered display: {len(displayed_ids)} of {len(rows)} tracks shown; {len(excluded)} hidden by the selected solver-failure / numerical-risk screen. This is a display filter, not a correction or proof of physical accuracy. Raw trajectories and outcomes are unchanged.</p><video src="animations/microbubble_physical_time.mp4?v=black-spatial-axes-filtered-20260929" controls loop></video><p><a href="figures/01_trajectories.pdf">Trajectory PDF</a> · <a href="data/display_selection.csv">Display selection</a> · <a href="data/render_manifest.json">Rendering provenance</a> · <a href="../POPULATION_VALIDATION_REPORT_ZH.md">Original simulation report and terminal outcomes</a></p><img src="figures/01_trajectories.png?v=black-spatial-axes-filtered-20260929">''')
print(json.dumps({k:record[k] for k in ['candidate_id','frames','fps','video_duration_s','displayed_particle_count','hidden_particle_count','preview_only']}),flush=True)
