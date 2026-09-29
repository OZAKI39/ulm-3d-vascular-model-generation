"""Only camera/font/backend preflight; not a trajectory result."""
from pathlib import Path
import os,json
os.environ['VTK_DEFAULT_OPENGL_WINDOW']='vtkEGLRenderWindow'
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont
from style_reference import projected_bounds,enlarge_animation
root=Path(__file__).resolve().parents[1];base=root.parent;g=json.loads((base/'inputs/geometry.json').read_text());config=json.loads((root/'data/environment.json').read_text());mesh=pv.read(config['flow']).extract_surface(algorithm='dataset_surface')
rows=[];out=root/'validation/layout';out.mkdir(exist_ok=True)
for pose in ['0','15']:
 surface=mesh.copy();T=np.array(g['pose_transforms'][pose]['source_mm_to_candidate_mm']);origin=np.array(g['origin_source_mm']);surface.points=(surface.points*1e3+origin)@T[:3,:3].T+T[:3,3]
 c=np.array(surface.center);r=np.linalg.norm(np.ptp(surface.points,axis=0))/2;d=np.array([-.72,.69,.11]);d/=np.linalg.norm(d)
 plot=pv.Plotter(off_screen=True,window_size=(1920,1080));plot.set_background('#101824');plot.add_mesh(surface,color='#b2becf',opacity=.16,smooth_shading=True)
 plot.camera.position=c+d*r*4;plot.camera.focal_point=c;plot.camera.up=(0,0,1);plot.camera.parallel_projection=True;plot.camera.parallel_scale=r*.78
 plot.show(auto_close=False);b=projected_bounds(plot,surface.points);span=b[1]-b[0];plot.camera.parallel_scale*=max(span[0]/1750,span[1]/840)*1.3
 layout=enlarge_animation(plot,surface);plot.render();img=Image.fromarray(plot.screenshot(return_img=True));draw=ImageDraw.Draw(img);font=ImageFont.truetype(str(root/'assets/NotoSansCJKsc-Regular.otf'),30)
 draw.text((40,22),'Layout preflight only | Actual BraVa geometry | Main elements x1.30',fill='white',font=font);img.save(out/f'candidate_{pose}_layout.png')
 gpu=[s for s in plot.ren_win.ReportCapabilities().splitlines() if any(k in s for k in ['vendor string','renderer string','version string'])]
 rows.append(dict(candidate=pose,layout=layout,gpu=gpu));plot.close()
(root/'validation/render_layout_preflight.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows),flush=True)
