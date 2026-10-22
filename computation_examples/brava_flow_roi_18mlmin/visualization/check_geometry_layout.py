"""Real-mesh camera/port preflight only: deliberately contains no flow field."""
from pathlib import Path
import os,json
os.environ.setdefault('VTK_DEFAULT_OPENGL_WINDOW','vtkEGLRenderWindow')
import numpy as np,pyvista as pv
from PIL import Image,ImageDraw,ImageFont
import style as shared
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'visualization/geometry_preflight';out.mkdir(exist_ok=True)
geometry=json.loads((ROOT/'inputs/geometry.json').read_text())
for pose in ['0','15']:
    T=np.array(geometry['pose_transforms'][pose]['source_mm_to_candidate_mm']);R=T[:3,:3]
    assert np.allclose(R.T@R,np.eye(3),atol=1e-14) and np.isclose(np.linalg.det(R),1)
    def x(p):return (np.asarray(p,dtype=float)*1000+geometry['origin_source_mm'])@R.T+T[:3,3]
    surface=pv.read(ROOT/'mesh/SV_MESH/mesh-complete.exterior.vtp');surface.points=x(surface.points)
    center=shared.fitted_z_center(surface.points);orbit=shared.Turntable(surface.points,[0,0,1],center)
    p=pv.Plotter(off_screen=True,window_size=(1920,1080),shape=(1,2),col_weights=[.9,.1],border=False)
    p.set_background('black',all_renderers=True);p.subplot(0,0);p.renderer.SetViewport(*shared.VIEWPORT);p.render_window.SetMultiSamples(0)
    p.add_mesh(surface,color='#b7c4d1',smooth_shading=True,ambient=.68,diffuse=.32,specular=.08)
    ports=[]
    for name,spec in geometry['ports'].items():
        c=pv.read(ROOT/'mesh/SV_MESH/mesh-surfaces'/f"{spec['role']}.vtp").center
        ports.append(('Inlet' if name=='I1' else 'Outlet 0'+name[-1],x(c)))
    orbit.setup(p);axes=shared.CoordinateGrid(p,surface.bounds,orbit.center,'full_vessel',(1920,1080))
    annotations=shared.AnnotationPlan(orbit,surface,axes,ports)
    p.show(auto_close=False,interactive=False)
    font=ImageFont.truetype(shared.FONT,32)
    for i in [0,108,216,324]:
        p.camera.position=orbit.trace[i]['camera_position_mm'];p.camera.up=orbit.trace[i]['camera_up_unit'];p.render()
        im=Image.fromarray(p.screenshot());draw=ImageDraw.Draw(im);draw.text((48,28),f'Geometry and port identity | Candidate {pose} | NO FLOW FIELD',font=font,fill='white')
        annotations.overlay(im,i,still=True).save(out/f'candidate_{pose}_{i:03d}.png',dpi=(300,300))
    result=dict(pose=int(pose),ports=[dict(name=n,position_mm=np.asarray(c).tolist()) for n,c in ports],
                backend=[s for s in p.render_window.ReportCapabilities().splitlines() if 'OpenGL vendor' in s or 'OpenGL renderer' in s],
                purpose='Geometry, coordinate and label preflight. No flow values used.',
                label_max_alpha_step=annotations.max_alpha_step)
    (out/f'candidate_{pose}.json').write_text(json.dumps(result,indent=2)+'\n');p.close()
    print(json.dumps(result),flush=True)
