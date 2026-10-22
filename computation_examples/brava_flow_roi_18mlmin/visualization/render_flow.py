"""Existing enlarged fixed-Z style, with real BraVa data and mm coordinate labels."""
from pathlib import Path
import argparse,json,os,time,hashlib
import numpy as np,pyvista as pv,imageio_ffmpeg
from PIL import Image,ImageDraw,ImageFont
import style as shared
ROOT=Path(__file__).resolve().parents[1]

def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def scalar_overlay(pixels,title,footer,unit,limits):
 im=Image.fromarray(pixels[:,:,:3]);w,h=im.size;d=ImageDraw.Draw(im)
 font=ImageFont.truetype(shared.FONT,round(h*.021));heading=ImageFont.truetype(shared.FONT,round(h*.032))
 d.text((w*.025,h*.025),title,font=heading,fill='#f0f3fa');d.text((w*.025,h*.958),footer,font=font,fill='#c0c8d3')
 x,y,bw,bh=round(w*.922),round(h*.275),round(w*.014),round(h*.49)
 colors=(shared.CMAP(np.linspace(1,0,bh))[:,:3]*255).astype(np.uint8);im.paste(Image.fromarray(np.repeat(colors[:,None,:],bw,axis=1)),(x,y))
 d.text((w*.904,h*.208),unit,font=font,fill='#f0f3fa');lo,hi=limits
 for val in np.linspace(lo,hi,6):
  yy=y+bh*(hi-val)/(hi-lo);d.line((x+bw,yy,x+bw+w*.005,yy),fill='#e9eef4',width=max(1,round(w/1920)));d.text((x+bw+w*.010,yy-h*.012),f'{val:g}',font=font,fill='#f0f3fa')
 return im

class SurfaceView:
 def __init__(self,mesh,surface,orbit,ports,kind,spec,size,annotations=None):
  self.orbit=orbit;self.spec=spec;self.p=pv.Plotter(off_screen=True,window_size=size,shape=(1,2),col_weights=[.9,.1],border=False)
  p=self.p;p.set_background('black',all_renderers=True);p.subplot(0,0);p.renderer.SetViewport(*shared.VIEWPORT);p.render_window.SetMultiSamples(0)
  p.add_mesh(mesh,scalars=spec['scalar'],preference=spec.get('preference','point'),cmap=shared.CMAP,clim=spec['clim'],smooth_shading=kind!='wss_raw',split_sharp_edges=False,ambient=.68,diffuse=.32,specular=.08,show_scalar_bar=False)
  orbit.setup(p);self.axes=shared.CoordinateGrid(p,surface.bounds,orbit.center,'full_vessel',size);self.annotations=annotations or shared.AnnotationPlan(orbit,surface,self.axes,ports);p.show(auto_close=False,interactive=False)
 def frame(self,i,still=False):
  self.p.camera.position=self.orbit.trace[i]['camera_position_mm'];self.p.camera.up=self.orbit.trace[i]['camera_up_unit'];self.p.render()
  return self.annotations.overlay(scalar_overlay(self.p.screenshot(),self.spec['title'],self.spec['footer'],self.spec['unit'],self.spec['clim']),i,still=still)
 def close(self):self.p.close()

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--stills-only',action='store_true');parser.add_argument('--pose',choices=['0','15','both'],default='both');args=parser.parse_args()
 active=json.loads((ROOT/'reports/ACTIVE_FLOW.json').read_text());case=ROOT/active['case'];prep=json.loads((ROOT/'visualization/data/PREPARATION.json').read_text());geometry=json.loads((ROOT/'inputs/geometry.json').read_text())
 assert prep['source_flow_sha256']==active['flow_sha256'], 'Streamline preparation belongs to a different flow'
 for relative,digest in active['postprocess_data_sha256'].items():assert hashlib.sha256((case/relative).read_bytes()).hexdigest()==digest,relative
 for name,digest in prep['data_sha256'].items():assert hashlib.sha256((ROOT/'visualization/data'/name).read_bytes()).hexdigest()==digest,name
 source={'pressure':pv.read(case/'postprocess/pressure_surface_si.vtp'),'wss':pv.read(case/'postprocess/wall_wss_si.vtp'),'lines':pv.read(ROOT/'visualization/data/streamlines_si.vtp'),'samples':pv.read(ROOT/'visualization/data/velocity_samples_si.vtp')}
 vmax=active['measurements']['velocity_max_m_s']*1000;shared.SPEED_MAX=float(np.ceil(vmax/100)*100)
 pmin,pmax=active['measurements']['pressure_range_pa'];climp=(float(min(0,np.floor(pmin/50)*50)),float(np.ceil(pmax/500)*500));wlim=(0.,max(5.,float(np.ceil(active['raw_WSS']['max_Pa']/5)*5)))
 specs={
 'pressure':dict(scalar='Pressure_Pa',title='Full vessel | Pressure field',footer='Steady FEM pressure | BraVa BG001 RMCA | Inlet: 18 mL/min',unit='Pressure (Pa)',clim=climp),
 'wss':dict(scalar='WSS_display_Pa',title='Full vessel | Wall shear stress',footer='FEM wall traction | Area-weighted nodal display of raw face magnitudes',unit='WSS (Pa)',clim=wlim),
 'wss_raw':dict(scalar='WSS_raw_Pa',preference='cell',title='Full vessel | Raw wall shear stress',footer='Raw P1 wall-triangle values | No scalar smoothing or normalization',unit='WSS (Pa)',clim=wlim)}
 for pose in (['0','15'] if args.pose=='both' else [args.pose]):
  out=ROOT/'visualization'/('candidate_'+pose)
  for name in ['figures','animations','inspection','data']:(out/name).mkdir(parents=True,exist_ok=True)
  T=np.array(geometry['pose_transforms'][pose]['source_mm_to_candidate_mm']);R=T[:3,:3];origin=np.array(geometry['origin_source_mm'])
  def transform(mesh):
   m=mesh.copy();m.points=(np.asarray(m.points,dtype=float)*1000+origin)@R.T+T[:3,3]
   for attr in [m.point_data,m.cell_data]:
    for name in ['Velocity','Velocity_m_s','Direction','Tangential_viscous_traction_Pa','Outward_normal']:
     if name in attr:attr[name]=np.asarray(attr[name])@R.T
   return m
  data={k:transform(v) for k,v in source.items()};surface=data['pressure'];ports=[]
  for name,spec in geometry['ports'].items():
   c=pv.read(ROOT/'mesh/SV_MESH/mesh-surfaces'/f"{spec['role']}.vtp").center;c=(np.array(c)*1000+origin)@R.T+T[:3,3]
   ports.append(('Inlet' if name=='I1' else 'Outlet 0'+name[-1],c))
  center=shared.fitted_z_center(surface.points);full_orbit=shared.Turntable(surface.points,[0,0,1],center)
  focus=(np.array(prep['branch_focus_m'])*1000+origin)@R.T+T[:3,3];detail=shared.crop_ellipsoid(surface,focus,np.array([6.,6.,6.]))
  samples=data['samples'];samples['Sample_id']=np.arange(samples.n_points)
  arrow=pv.Arrow(start=(-.5,0,0),direction=(1,0,0),tip_length=.32,tip_radius=.12,shaft_radius=.035,tip_resolution=10,shaft_resolution=8)
  glyph=samples.glyph(orient='Direction',scale=False,factor=.6,geom=arrow);inside=glyph.select_enclosed_points(surface,tolerance=1e-8,check_surface=True)
  invalid=np.unique(glyph['Sample_id'][inside['SelectedPoints']==0]);keep=~np.isin(samples['Sample_id'],invalid);samples=samples.extract_points(keep).extract_surface(algorithm='dataset_surface')
  glyph=samples.glyph(orient='Direction',scale=False,factor=.6,geom=arrow);assert samples.n_points>20;glyph.save(out/'data/velocity_glyphs_mm.vtp')
  local_orbit=shared.Turntable(np.vstack([detail.points,glyph.points]),[0,0,1],focus)
  annotations={};records=[]
  kinds=['streamlines','velocity_vectors','pressure','wss','wss_raw']
  def make(kind,size):
   if kind in specs:return SurfaceView(data['pressure' if kind=='pressure' else 'wss'],surface,full_orbit,ports,kind,specs[kind],size,annotations.get('full'))
   if kind=='streamlines':return shared.View(surface,data['lines'],full_orbit,size=size,kind='full_vessel',count=96,ports=ports,annotations=annotations.get('full'))
   return shared.View(detail,glyph,local_orbit,size=size,kind='branch_vectors',count=samples.n_points,ports=[],annotations=annotations.get('local'))
  for kind in kinds:
   view=make(kind,(3840,2160));key='local' if kind=='velocity_vectors' else 'full';annotations[key]=view.annotations
   image=view.frame(0,still=True);image.save(out/'figures'/f'{kind}_overview_4k.png',dpi=(300,300));image.save(out/'figures'/f'{kind}_overview.pdf',resolution=300)
   gpu=[s for s in view.p.render_window.ReportCapabilities().splitlines() if 'OpenGL vendor' in s or 'OpenGL renderer' in s];view.close()
   if args.stills_only:continue
   view=make(kind,shared.SIZE);path=out/'animations'/f'{kind}_rotation.mp4';tmp=path.with_name(path.stem+'.tmp.mp4')
   encoder=imageio_ffmpeg.write_frames(str(tmp),size=shared.SIZE,fps=shared.FPS,codec='libx264',quality=8,pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=8,output_params=['-preset','medium','-threads','4','-movflags','+faststart']);encoder.send(None)
   try:
    for i in range(shared.FRAMES):
     frame=view.frame(i);encoder.send(np.ascontiguousarray(frame))
     if i in shared.KEYS:frame.save(out/'inspection'/f'{kind}_{i:04d}.png')
     if i%108==0:print(pose,kind,i,flush=True)
   finally:encoder.close();view.close()
   tmp.replace(path);reader=imageio_ffmpeg.read_frames(str(path));meta=next(reader);count=sum(1 for _ in reader);assert meta['size']==shared.SIZE and count==shared.FRAMES
   records.append(dict(kind=kind,file=str(path.relative_to(out)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),frames=count,OpenGL=gpu))
  dump(out/('PREVIEW.json' if args.stills_only else 'MEDIA_VALIDATION.json'),dict(pose=int(pose),source_flow_sha256=active['flow_sha256'],coordinate_unit='mm',physical_pose_transform=T.tolist(),frames=shared.FRAMES,fps=shared.FPS,style='Original enlarged fixed-Z turntable with side-arrow fading labels',speed_range_mm_s=[0,shared.SPEED_MAX],pressure_range_Pa=climp,wss_range_Pa=wlim,raw_values_unchanged=True,media=records,display_glyph_count=samples.n_points,removed_outside_glyphs=len(invalid)))
  html='<!doctype html><meta charset="utf-8"><title>BraVa flow</title><style>body{background:#101824;color:#edf2fa;max-width:1300px;margin:30px auto;font:18px sans-serif}video,img{width:100%}a{color:#50d9ef}</style><h1>BraVa BG001 RMCA | candidate '+pose+'</h1><p>Inlet: 18 mL/min. Common physical solution, rigidly transformed to the selected print pose.</p>'
  for k in kinds:html+=f'<h2>{k.replace("_"," ")}</h2><video controls src="animations/{k}_rotation.mp4"></video><a href="figures/{k}_overview.pdf">PDF</a>'
  (out/'OPEN_RESULTS.html').write_text(html)
if __name__=='__main__':main()
