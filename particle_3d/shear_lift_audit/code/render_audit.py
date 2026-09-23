"""Same frozen-Z presentation as the existing flow views; audit data only."""
from pathlib import Path
import json,argparse
import numpy as np
import pyvista as pv
import imageio_ffmpeg
from PIL import Image,ImageDraw,ImageFont
import shared_rotation as shared

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/remote_results';FIG=ROOT/'figures'
CLIM=(1.,20000.)
def overlay(pixels,title,footer=True,branch=False):
    im=Image.fromarray(pixels[:,:,:3]);w,h=im.size;draw=ImageDraw.Draw(im)
    font=ImageFont.truetype(shared.FONT,round(h*.020));big=ImageFont.truetype(shared.FONT,round(h*.030))
    draw.text((w*.023,h*.018),title,font=big,fill='#f0f3fa')
    if footer:
        draw.text((w*.023,h*.943),'Candidate Saffman magnitude only',font=font,fill='#c4cad3')
        draw.text((w*.023,h*.973),'Not included in trajectory dynamics',font=font,fill='#c4cad3')
    x,y,bw,bh=round(w*.923),round(h*.275),round(w*.013),round(h*.49)
    colors=(shared.CMAP(np.linspace(1,0,bh))[:,:3]*255).astype(np.uint8)
    im.paste(Image.fromarray(np.repeat(colors[:,None,:],bw,axis=1)),(x,y))
    draw.text((w*.901,h*.207),'Shear rate (s⁻¹)',font=font,fill='white')
    for v in 10.**np.arange(np.log10(CLIM[0]),np.log10(CLIM[1])+.01):
        yy=y+bh*(np.log10(CLIM[1])-np.log10(v))/(np.log10(CLIM[1])-np.log10(CLIM[0]))
        draw.line((x+bw,yy,x+bw+w*.004,yy),fill='white',width=max(1,round(w/1920)))
        draw.text((x+bw+w*.008,yy-h*.011),f'{v:g}',font=font,fill='white')
    if branch:
        draw.text((w*.023,h*.09),'CANDIDATE ONLY | Arrows show unit directions',font=font,fill='#e3d0b2')
        draw.text((w*.023,h*.128),'Cyan: aligned slip     Orange: slip × vorticity',font=font,fill='#d6e1ef')
    return im

class ShearView:
    def __init__(self,mesh,ports,title,size=(1920,1080),plan=None,points=None):
        self.title=title;self.points=points;self.p=pv.Plotter(off_screen=True,window_size=size,shape=(1,2),col_weights=[.90,.10],border=False)
        p=self.p;p.set_background('black',all_renderers=True);p.subplot(0,0);p.renderer.SetViewport(*shared.VIEWPORT);p.render_window.SetMultiSamples(0)
        p.add_mesh(mesh,scalars='local_shear_rate_s_inv',preference='cell',cmap=shared.CMAP,clim=CLIM,log_scale=True,
                   show_scalar_bar=False,ambient=.65,diffuse=.35,specular=.05)
        # Sparse real saved positions only, no new trajectory or particle integration.
        if points is not None:
            for point,radius in points:p.add_mesh(pv.Sphere(radius=radius,center=point,theta_resolution=14,phi_resolution=10),color='#fff1ad')
        self.orbit=shared.Turntable(mesh.points,[0.,0.,1.],shared.fitted_z_center(mesh.points))
        self.orbit.setup(p);self.axes=shared.CoordinateGrid(p,mesh.bounds,self.orbit.center,'full_vessel',size)
        self.plan=plan or shared.AnnotationPlan(self.orbit,mesh,self.axes,ports)
        p.show(auto_close=False,interactive=False)
    def frame(self,i,still=False):
        c=self.orbit.trace[i];self.p.camera.position=c['camera_position_um'];self.p.camera.up=c['camera_up_unit'];self.p.render()
        assert np.allclose(self.p.camera.up,[0,0,1])
        im=self.plan.overlay(overlay(self.p.screenshot(),self.title),i,still=still)
        # Screen-space rings mark true saved positions through the opaque surface.
        # These are location markers, not enlarged physical bubble geometry.
        if self.points:
            draw=ImageDraw.Draw(im);scale=im.width/shared.SIZE[0]
            xy=self.plan.project(np.array([point for point,radius in self.points]),i)*scale
            for x,y in xy:
                r=5*scale
                draw.ellipse((x-r,y-r,x+r,y+r),outline='#fff5db',width=max(1,round(2*scale)))
                draw.ellipse((x-scale,y-scale,x+scale,y+scale),fill='#fff5db')
            font=ImageFont.truetype(shared.FONT,round(im.height*.018))
            draw.text((im.width*.023,im.height*.062),'Circles: saved MB locations',font=font,fill='#d4cdbb')
        return im

def main():
    global CLIM
    parser=argparse.ArgumentParser();parser.add_argument('--stills-only',action='store_true');parser.add_argument('--branch-only',action='store_true');args=parser.parse_args()
    FIG.mkdir(exist_ok=True);(ROOT/'animations').mkdir(exist_ok=True);(ROOT/'data/inspection').mkdir(exist_ok=True)
    limits=[]
    for flow in ['FLOW_0P353_MMPS','FLOW_2P0_MMPS']:
        with np.load(DATA/(flow+'_eulerian.npz')) as z:
            g=z['shear_s_inv'];limits.extend([g[g>0].min(),g.max()])
    CLIM=(float(10**np.floor(np.log10(min(limits)))),float(10**np.ceil(np.log10(max(limits)))))
    if args.branch_only:
        old=pv.read(DATA/'FLOW_0P353_MMPS_shear_surface.vtp');old.points=np.asarray(old.points,float)*1e6
        record=json.loads((ROOT/'data/visualization_validation.json').read_text())
        render_branch(old,record);return
    fem=ROOT.parents[1]/'formal_3D_flow_solver/FEM_SimVascular/frozen_reference/SV_MESH/mesh-surfaces'
    ports=[]
    for name in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        cap=pv.read(fem/(name+'.vtp'));ports.append(('Inlet' if name=='INLET' else 'Outlet '+name[-2:],np.asarray(cap.center)*1e6))
    reps=json.loads((DATA/'representatives.json').read_text());points=[]
    for r in reps:
        with np.load(DATA/'representatives'/f"{r['dataset']}_{r['id']:06d}.npz") as z:
            k=len(z['time_s'])//2;points.append((z['evaluation_x_m'][k]*1e6,z['radius_m'][k]*1e6))
    plan=None;record={}
    if args.stills_only and (ROOT/'data/visualization_validation.json').exists():
        previous=json.loads((ROOT/'data/visualization_validation.json').read_text())
        if 'video' in previous:record['video']=previous['video']
    for flow,tag,title in [('FLOW_0P353_MMPS','01_local_shear_rate_map','Full vessel | Local FEM shear | 0.353 mm/s'),
                            ('FLOW_2P0_MMPS','01b_new_flow_shear_rate_map','Full vessel | Local FEM shear | 2.0 mm/s')]:
        mesh=pv.read(DATA/(flow+'_shear_surface.vtp'));mesh.points=np.asarray(mesh.points,float)*1e6
        view=ShearView(mesh,ports,title,(3840,2160),plan=plan,points=points if '0P353' in flow else None)
        view.frame(0,still=True).save(FIG/(tag+'.png'));plan=view.plan
        record[flow]=dict(axis=[0,0,1],center_um=view.orbit.center.tolist(),scale_um=view.orbit.scale,log_color_range=CLIM,
                         sample_display='ELEMENT_LOCAL_GRADIENT_ON_EXTERIOR_FACES; NO_WSS_CONVERSION',camera=view.orbit.trace)
        view.p.close();print('STILL',tag,flush=True)
    old=pv.read(DATA/'FLOW_0P353_MMPS_shear_surface.vtp');old.points=np.asarray(old.points,float)*1e6
    if not args.stills_only:
        view=ShearView(old,ports,'Full vessel | Local FEM shear | 0.353 mm/s',plan=plan,points=points)
        path=ROOT/'animations/shear_lift_audit_rotating.mp4'
        writer=imageio_ffmpeg.write_frames(str(path),shared.SIZE,fps=shared.FPS,codec='libx264',quality=8,macro_block_size=2,output_params=['-movflags','+faststart']);writer.send(None)
        for i in range(shared.FRAMES):
            im=view.frame(i);writer.send(np.asarray(im))
            if i in shared.KEYS:im.save(ROOT/'data/inspection'/f'rotation_{i:04d}.png')
            if i%108==0:print('FRAME',i,flush=True)
        writer.close();view.plan.export(ROOT/'data/rotation_annotations.json');view.p.close()
        reader=imageio_ffmpeg.read_frames(str(path));metadata=next(reader);n=0;distinct=set()
        for raw in reader:n+=1;distinct.add(shared.hashlib.sha256(raw).hexdigest())
        assert n==shared.FRAMES and len(distinct)==n and tuple(metadata['size'])==tuple(shared.SIZE)
        record['video']=dict(frames=n,distinct_frames=len(distinct),fps=shared.FPS,duration_s=n/shared.FPS,size=shared.SIZE,
                             content='CAMERA_ROTATION_ONLY_OVER_FROZEN_FIELD_AND_STATIC_SAVED_POSITIONS',
                             saved_position_marker_role='SCREEN_SPACE_LOCATION_RINGS_NOT_PHYSICAL_SIZE',sha256=shared.sha(path))
    render_branch(old,record)

def render_branch(old,record):
    # Branch direction diagnostic, using resolved near-field samples and fixed-size unit arrows.
    with np.load(DATA/'plot_points.npz') as f:d={k:f[k] for k in f.files}
    center=np.array([92.1310784,48.2637661,112.1937748]);x=d['evaluation_x_m']*1e6
    ids=np.flatnonzero((np.linalg.norm(x-center,axis=1)<10)&np.isfinite(d['candidate_direction_only']).all(axis=1))
    ids=ids[np.argsort(d['CANDIDATE_SAFFMAN_MAGNITUDE'][ids])[::-1]];chosen=[]
    for i in ids:
        if not chosen or np.min(np.linalg.norm(x[chosen]-x[i],axis=1))>1.4:chosen.append(int(i))
        if len(chosen)>=14:break
    local=shared.crop_ellipsoid(old,center,np.array([17.,17.,13.]))
    p=pv.Plotter(off_screen=True,window_size=(2560,1440));p.set_background('black');p.renderer.SetViewport(*shared.VIEWPORT)
    p.add_mesh(local,scalars='local_shear_rate_s_inv',cmap=shared.CMAP,clim=CLIM,log_scale=True,opacity=.26,show_scalar_bar=False)
    for color,key in [('#37d9ed','slip_m_s'),('#ffb35b','candidate_direction_only')]:
        if chosen:
            vectors=d[key][chosen];vectors=vectors/np.linalg.norm(vectors,axis=1)[:,None]
            pts=pv.PolyData(x[chosen]);pts['direction']=vectors
            p.add_mesh(pts.glyph(orient='direction',scale=False,factor=2.3),color=color,ambient=.7)
    p.add_axes(xlabel='X',ylabel='Y',zlabel='Z',color='white')
    bounds=np.asarray(local.bounds).reshape(3,2)
    # Explicit world-coordinate grid avoids VTK empty-title/font artifacts.
    for fixed_axis in range(3):
        other=[j for j in range(3) if j!=fixed_axis]
        for running in other:
            crossed=next(j for j in other if j!=running)
            step=10. if crossed==0 else 5.
            for value in np.arange(np.ceil(bounds[crossed,0]/step)*step,bounds[crossed,1],step):
                a=bounds[:,0].copy();a[fixed_axis]=bounds[fixed_axis,1 if fixed_axis==1 else 0]
                a[crossed]=value;b=a.copy();b[running]=bounds[running,1]
                p.add_mesh(pv.Line(a,b),color='#536579',opacity=.35,line_width=1)
    for axis in range(3):
        a=bounds[:,0].copy()
        if axis==1:a[0]=bounds[0,1]
        b=a.copy();b[axis]=bounds[axis,1]
        p.add_mesh(pv.Line(a,b),color='#8b9cb0',line_width=1)
    p.camera.position=center+np.array([32.,-32.,7.]);p.camera.focal_point=center;p.camera.up=[0,0,1];p.camera.parallel_projection=True;p.camera.parallel_scale=18.
    p.show(auto_close=False,interactive=False)
    im=overlay(p.screenshot(),'Branch detail | Candidate shear-lift audit',branch=True)
    # Use the same font as the full-vessel labels. VTK's default axes font drops µ.
    draw=ImageDraw.Draw(im);font=ImageFont.truetype(shared.FONT,26);titlefont=ImageFont.truetype(shared.FONT,29)
    def project(point):
        p.renderer.SetWorldPoint(*point,1.);p.renderer.WorldToDisplay();q=p.renderer.GetDisplayPoint()
        return np.array([q[0],im.height-q[1]])
    mid=project(center)
    for axis,step in [(0,10.),(1,5.),(2,5.)]:
        start=bounds[:,0].copy()
        if axis==1:start[0]=bounds[0,1]
        end=start.copy();end[axis]=bounds[axis,1]
        pa,pb=project(start),project(end);along=pb-pa;normal=np.array([-along[1],along[0]])/np.linalg.norm(along)
        if normal@((pa+pb)/2-mid)<0:normal=-normal
        for value in np.arange(np.ceil(bounds[axis,0]/step)*step,bounds[axis,1],step):
            point=start.copy();point[axis]=value;q=project(point)
            draw.line([tuple(q),tuple(q+normal*7)],fill='#c0ccd9',width=2)
            draw.text(tuple(q+normal*25),f'{value:g}',font=font,fill='#c0ccd9',anchor='mm')
        draw.text(tuple((pa+pb)/2+normal*(100 if axis==2 else 66)),'XYZ'[axis]+' (µm)',font=titlefont,fill='#d6e1ef',anchor='mm')
    im.save(FIG/'09_branch_point_audit.png');p.close()
    record['branch']=dict(arrows=len(chosen),indices=chosen,glyph_length_um=2.3,glyph_role='NORMALIZED_DIRECTION_ONLY_NOT_FORCE_MAGNITUDE')
    shared.dump(ROOT/'data/visualization_validation.json',record)

if __name__=='__main__':main()
