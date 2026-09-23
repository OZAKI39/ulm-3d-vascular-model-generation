"""Enlarged full-vessel streamlines and branch-detail FEM velocity glyphs.

The center, axis and magnification stay fixed throughout each turntable.
The camera orbits the fixed axis; physical geometry and its coordinate axes stay
in their original coordinate frame, with no transformation of the source data.
"""
import os
os.environ.setdefault('LP_NUM_THREADS', '2')
os.environ.setdefault('VTK_SMP_MAX_THREADS', '2')
import argparse
import csv
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

import imageio_ffmpeg
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from scipy.optimize import minimize
from scipy.spatial import ConvexHull, cKDTree
from scipy.ndimage import minimum_filter1d, convolve1d
from vtkmodules.vtkCommonDataModel import vtkQuadric
from vtkmodules.vtkFiltersCore import vtkClipPolyData
from vtkmodules.vtkFiltersCore import vtkImplicitPolyDataDistance

import matplotlib
from matplotlib.colors import ListedColormap

HERE = Path(__file__).resolve().parent
CASE = HERE / 'input_data'
CMAP = ListedColormap(matplotlib.colormaps['turbo'](np.linspace(.10, .96, 256)))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def axis_rotation(axis, angle_degrees):
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    skew = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    angle = np.deg2rad(angle_degrees)
    return np.eye(3)*np.cos(angle)+(1-np.cos(angle))*np.outer(axis, axis)+np.sin(angle)*skew


OUT = HERE / 'results'
SIZE = (1920, 1080)
FPS, FRAMES = 24, 432
VIEWPORT = (0., .070, .90, .930)
VIEW_HEIGHT = VIEWPORT[3] - VIEWPORT[1]
ASPECT = SIZE[0] * .90 / (SIZE[1] * VIEW_HEIGHT)
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
TITLE = 'Full vessel | FEM streamlines'
DETAIL_RADII_UM = np.array([18.,18.,11.5])
ELEVATION_DEG = 6.
KEYS = [0, 108, 216, 324, 431]






def fitted_z_center(points):
    """Maximize a constant full-vessel zoom for a fixed vertical Z orbit."""
    hull=points[ConvexHull(points).vertices]
    center=(points.min(0)+points.max(0))/2
    phi=np.deg2rad(ELEVATION_DEG)
    def constraints(x):
        q=hull-x[:3]
        radial=np.linalg.norm(q[:,:2],axis=1)
        return np.r_[x[3]-q[:,2]*np.cos(phi)-radial*np.sin(phi),
                     x[3]+q[:,2]*np.cos(phi)-radial*np.sin(phi),
                     x[3]-radial/ASPECT]
    fit=minimize(lambda x:x[3],np.r_[center,80.],method='SLSQP',
        constraints={'type':'ineq','fun':constraints},options={'ftol':1e-10,'maxiter':500})
    assert fit.success,fit.message
    return fit.x[:3]


def velocity_glyphs(field, full_surface, focus):
    """Sample true P1 vectors at interior tetrahedron centroids, then thin in space."""
    tet = field.cells.reshape(-1, 5)
    assert np.all(tet[:, 0] == 4)
    nodes = tet[:, 1:]
    positions = np.asarray(field.points,dtype=float)[nodes].mean(1)*1e6
    velocity = field['Velocity'][nodes].mean(1)
    speed = np.linalg.norm(velocity, axis=1)*1000
    use = np.flatnonzero((np.linalg.norm((positions-focus)/DETAIL_RADII_UM, axis=1) < .92) & (speed > .15))
    distance = vtkImplicitPolyDataDistance()
    distance.SetInput(full_surface)
    clearance = np.array([abs(distance.EvaluateFunction(p)) for p in positions[use]])
    use = use[clearance > .40]
    # Greedy farthest-point selection avoids clusters from irregular tetra sizes.
    candidates = positions[use]
    chosen = [int(np.argmin(np.linalg.norm(candidates-focus, axis=1)))]
    nearest = np.linalg.norm(candidates-candidates[chosen[0]], axis=1)
    while len(chosen) < 350 and nearest.max() > 1.25:
        j = int(nearest.argmax())
        chosen.append(j)
        nearest = np.minimum(nearest, np.linalg.norm(candidates-candidates[j], axis=1))
    ids = use[chosen]
    samples = pv.PolyData(positions[ids])
    samples['Velocity_m_s'] = velocity[ids]
    samples['Direction'] = velocity[ids] / np.linalg.norm(velocity[ids], axis=1)[:, None]
    samples['Speed_mm_s'] = speed[ids]
    samples['Tetra_id'] = ids
    samples['Sample_id'] = np.arange(len(ids))
    arrow = pv.Arrow(start=(-.5,0,0), direction=(1,0,0), scale=1.,
                     tip_length=.32, tip_radius=.12, shaft_radius=.035,
                     tip_resolution=10, shaft_resolution=8)
    glyphs = samples.glyph(orient='Direction', scale=False, factor=1.55, geom=arrow)
    # Suppress only glyphs with any display vertex outside the actual lumen.
    enclosed = glyphs.select_enclosed_points(full_surface, tolerance=1e-8, check_surface=True)
    invalid = np.unique(glyphs['Sample_id'][enclosed['SelectedPoints'] == 0])
    keep = ~np.isin(samples['Sample_id'], invalid)
    samples = samples.extract_points(keep).extract_surface(algorithm='dataset_surface')
    glyphs = samples.glyph(orient='Direction', scale=False, factor=1.55, geom=arrow)
    assert samples.n_points > 40 and glyphs.n_points > 1000
    assert np.allclose(samples['Direction'], samples['Velocity_m_s']/np.linalg.norm(samples['Velocity_m_s'],axis=1)[:,None])
    samples.save(OUT/'data/velocity_glyph_samples_um.vtp')
    glyphs.save(OUT/'data/velocity_glyph_geometry_um.vtp')
    with (OUT/'data/velocity_glyph_samples.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['sample_id','tetra_id','x_um','y_um','z_um','vx_m_s','vy_m_s','vz_m_s','speed_mm_s'])
        writer.writerows([int(i),int(t),*p,*v,float(s)] for i,t,p,v,s in zip(
            samples['Sample_id'],samples['Tetra_id'],samples.points,samples['Velocity_m_s'],samples['Speed_mm_s']))
    return glyphs, samples.n_points, dict(method='P1 velocity at tetrahedron centroid (mean of 4 nodal vectors)',
        length_um=1.55, length_role='Constant display length; direction from actual velocity; color encodes speed',
        minimum_sampling_spacing_um=1.25, rejected_outside_lumen_glyphs=len(invalid),
        actual_vector_data_exported=True, candidate_count=len(use), glyph_count=samples.n_points)


def crop_ellipsoid(mesh, center, radii):
    sphere = vtkQuadric()
    diag=1/np.square(radii)
    sphere.SetCoefficients(*diag,0.,0.,0.,*(-2*center*diag),float(np.sum(center*center*diag)-1))
    clip = vtkClipPolyData()
    clip.SetInputData(mesh)
    clip.SetClipFunction(sphere)
    clip.InsideOutOn()
    clip.Update()
    result = pv.wrap(clip.GetOutput()).copy(deep=True)
    assert result.n_points > 100
    assert np.linalg.norm((result.points-center)/radii,axis=1).max()<=1.0001
    return result


class Turntable:
    def __init__(self, points, axis, center, start_angle=90.):
        self.points = points
        self.axis = np.asarray(axis, dtype=float)
        self.axis /= np.linalg.norm(self.axis)
        self.center = np.asarray(center, dtype=float)
        assert np.array_equal(self.axis,[0.,0.,1.])
        phi=np.deg2rad(ELEVATION_DEG)
        self.direction=np.array([np.cos(phi)/np.sqrt(2),np.cos(phi)/np.sqrt(2),np.sin(phi)])
        self.right=np.cross(self.axis,self.direction)
        self.right/=np.linalg.norm(self.right)
        self.up=np.cross(self.direction,self.right)
        q = points - self.center
        height = q @ self.axis
        radial = np.linalg.norm(q - height[:, None] * self.axis, axis=1)
        vertical=np.abs(height)*np.cos(phi)+radial*np.sin(phi)
        self.scale=float(max(vertical.max(),radial.max()/ASPECT)/.975)
        self.radius = float(np.linalg.norm(q, axis=1).max())
        self.position = self.center + self.direction * 4.5 * self.radius
        self.matrices, self.trace = [], []
        for i in range(FRAMES):
            angle = start_angle + 360 * i / FRAMES
            rotation = axis_rotation(self.axis, angle)
            matrix = np.eye(4)
            matrix[:3, :3] = rotation
            matrix[:3, 3] = self.center - rotation @ self.center
            transformed = q @ rotation.T
            x = transformed @ self.right / (self.scale * ASPECT)
            y = transformed @ self.up / self.scale
            bounds = [float(x.min()), float(x.max()), float(y.min()), float(y.max())]
            assert max(abs(v) for v in bounds) <= .975001
            assert abs(self.right@self.axis)<1e-14, 'Z must project vertically'
            assert np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-13)
            assert np.isclose(np.linalg.det(rotation), 1.)
            assert np.allclose(rotation @ self.axis, self.axis, atol=1e-13)
            self.matrices.append(matrix)
            self.trace.append(dict(frame=i, video_time_s=i/FPS, angle_deg=angle,
                axis_unit=self.axis.tolist(), center_um=self.center.tolist(),
                camera_position_um=(self.center+rotation.T@(self.position-self.center)).tolist(),
                camera_up_unit=self.axis.tolist(),
                parallel_scale_um=self.scale, projected_bounds=bounds, display_matrix=matrix.tolist()))

    def setup(self, plotter):
        plotter.camera.position = self.position
        plotter.camera.focal_point = self.center
        plotter.camera.up = self.axis
        plotter.camera.parallel_projection = True
        plotter.camera.parallel_scale = self.scale
        plotter.camera.clipping_range = (.01, 1000.)


def annotate(pixels, count, kind):
    """Static typography on the rendered image; no opaque text backgrounds."""
    image = Image.fromarray(pixels[:, :, :3])
    w, h = image.size
    draw = ImageDraw.Draw(image)
    title = ImageFont.truetype(FONT, round(h*.032))
    font = ImageFont.truetype(FONT, round(h*.021))
    draw.text((w*.025, h*.025), TITLE if kind=='full_vessel' else 'Local detail | FEM velocity vectors', font=title, fill='#f0f3fa')
    draw.text((w*.025, h*.958), f'{count} '+('streamlines' if kind=='full_vessel' else 'vectors'), font=font, fill='#c0c8d3')
    x, y, bw, bh = round(w*.922), round(h*.275), round(w*.014), round(h*.49)
    colors = (CMAP(np.linspace(1, 0, bh))[:, :3]*255).astype(np.uint8)
    image.paste(Image.fromarray(np.repeat(colors[:, None, :], bw, axis=1)), (x, y))
    draw.text((w*.904, h*.208), 'Speed (mm/s)', font=font, fill='#f0f3fa')
    for value in np.linspace(0, 7.5, 6):
        yy = y + bh * (1 - value/7.5)
        draw.line((x+bw, yy, x+bw+w*.005, yy), fill='#e9eef4', width=max(1, round(w/1920)))
        draw.text((x+bw+w*.010, yy-h*.012), f'{value:g}', font=font, fill='#f0f3fa')
    return image


class AnnotationPlan:
    """Plan a complete periodic rotation before drawing any text.

    Every candidate follows a continuous projection, never a frame-local greedy
    relocation. Switches between candidates use overlapping raised-cosine fades.
    Validity is eroded in time before smoothing, so a fading port label cannot
    enter a vessel or cover another label. The same plan is used at 1080p/4K.
    """
    FADE_RADIUS = 10  # 21 samples / 0.875 s; includes the last-to-first seam.

    def __init__(self, orbit, wall, axes, ports):
        self.orbit = orbit
        self.tracks = []
        self.projected_wall = [self.project(wall.points, i) for i in range(FRAMES)]
        self.trees = [cKDTree(p) for p in self.projected_wall]
        self.kernel = np.hanning(2*self.FADE_RADIUS+3)[1:-1]
        self.kernel /= self.kernel.sum()
        self.anchors = (np.stack([self.project([p[1] for p in ports], i)
                        for i in range(FRAMES)]) if ports else None)
        if ports:
            # The larger terminal branches get space first; all names stay
            # associated with the original cap centers throughout the rotation.
            for j in [0, 3, 1, 2]:
                name = ports[j][0]
                anchor = self.anchors[:,j]
                offsets = np.array([(dx,dy) for dx in [-180,180,-260,260,-360,360,0]
                                    for dy in [0,-70,70,-140,140,-210,210]],float)
                centers = anchor[:,None,:]+offsets[None,:,:]
                width,height = self.text_size(name,25)
                centers[:,:,0] = np.clip(centers[:,:,0],45+width/2,1690-width/2)
                centers[:,:,1] = np.clip(centers[:,:,1],110+height/2,975-height/2)
                clearance = self.clearance(centers,width,height)
                valid = clearance >= 22.
                for other in range(len(ports)):
                    distance = np.maximum(np.abs(centers-self.anchors[:,other,None,:])-
                                          [width/2,height/2],0.)
                    valid &= np.linalg.norm(distance,axis=2)>=55.
                valid &= self.unoccupied(centers,width,height,pad=16.)
                cost = np.linalg.norm(centers-anchor[:,None,:],axis=2)/100.
                cost += .002*np.abs(centers[:,:,1]-anchor[:,None,1])
                track = self.add(name,25,centers,valid,cost,'port',
                                 anchor=anchor,clearance=clearance)
                assert track is not None
                print('Planned clear arrow label:',name,flush=True)
        for axis,ends in enumerate(axes.axis_ends):
            pairs=np.stack([self.project(ends,i) for i in range(FRAMES)])
            direction=pairs[:,1]-pairs[:,0]
            direction/=np.linalg.norm(direction,axis=1)[:,None]
            normal=np.column_stack([-direction[:,1],direction[:,0]])
            positions=[]
            for t in [1.,.8,1.2,.6,.4,.2,-.2,-.4,-.6,-.8]:
                point=axes.origin+(ends[1]-axes.origin)*t
                projected=np.stack([self.project([point],i)[0] for i in range(FRAMES)])
                for offset in [(0,-54),(0,54),(-75,0),(75,0),
                               (-140,0),(140,0),(0,-110),(0,110),
                               (-210,0),(210,0),(0,-180),(0,180)]:
                    positions.append(projected+offset)
            centers=np.stack(positions,axis=1)
            text='XYZ'[axis]+' (µm)'
            width,height=self.text_size(text,24)
            valid=self.in_view(centers,width,height)
            valid &= self.unoccupied(centers,width,height,pad=12.)
            valid &= self.clearance(centers,width,height)>9.
            cost=np.broadcast_to(np.arange(centers.shape[1])*.045,centers.shape[:2]).copy()
            self.add(text,24,centers,valid,cost,'axis_title')
        # Keep the identity and world coordinate of every tick. Crowded labels
        # fade at their own projected location rather than being moved elsewhere.
        for axis,ends in enumerate(axes.axis_ends):
            pairs=np.stack([self.project(ends,i) for i in range(FRAMES)])
            direction=pairs[:,1]-pairs[:,0]
            direction/=np.linalg.norm(direction,axis=1)[:,None]
            normal=np.column_stack([-direction[:,1],direction[:,0]])
            for value in axes.ticks[axis]:
                point=axes.origin.copy();point[axis]=value
                anchors=np.stack([self.project([point],i)[0] for i in range(FRAMES)])
                centers=(anchors+normal*23.)[:,None,:]
                width,height=self.text_size(f'{value:g}',18)
                valid=self.in_view(centers,width,height)
                valid &= self.unoccupied(centers,width,height,pad=7.)
                valid &= self.clearance(centers,width,height)>5.
                safe=minimum_filter1d(valid.astype(float),size=2*self.FADE_RADIUS+1,axis=0,mode='wrap')
                alpha=convolve1d(safe,self.kernel,axis=0,mode='wrap')
                # Tick marks remain on the physical axis even when crowded
                # text fades away; only the viewport edge affects their fade.
                mark_valid=self.in_view(anchors[:,None,:],10.,10.).astype(float)
                mark_safe=minimum_filter1d(mark_valid,size=2*self.FADE_RADIUS+1,axis=0,mode='wrap')
                mark_alpha=convolve1d(mark_safe,self.kernel,axis=0,mode='wrap')[:,0]
                self.tracks.append(dict(text=f'{value:g}',font_size=18,kind='tick',
                    centers=centers,alpha=alpha,width=width,height=height,
                    anchor=anchors,normal=normal,axis=axis,value=float(value),mark_alpha=mark_alpha))
        self.audit()

    @staticmethod
    def text_size(text, size):
        box=ImageFont.truetype(FONT,size).getbbox(text,anchor='mm')
        return float(max(abs(box[0]),abs(box[2]))*2+8),float(max(abs(box[1]),abs(box[3]))*2+8)

    def project(self, points, i):
        matrix=self.orbit.matrices[i]
        q=np.asarray(points)@matrix[:3,:3].T+matrix[:3,3]-self.orbit.center
        x=q@self.orbit.right/(self.orbit.scale*ASPECT)
        y=q@self.orbit.up/self.orbit.scale
        return np.column_stack(((x+1)*.5*SIZE[0]*.90,
            (1-VIEWPORT[1]-(y+1)*.5*VIEW_HEIGHT)*SIZE[1]))

    @staticmethod
    def in_view(centers,width,height):
        return ((centers[:,:,0]-width/2>35)&(centers[:,:,0]+width/2<1690)&
                (centers[:,:,1]-height/2>110)&(centers[:,:,1]+height/2<975))

    def clearance(self,centers,width,height):
        # Dense samples of the entire text rectangle give a conservative bound
        # on distance to all projected surface vertices (spacing <= 6 px).
        xx=np.linspace(-width/2,width/2,max(2,int(np.ceil(width/6))+1))
        yy=np.linspace(-height/2,height/2,max(2,int(np.ceil(height/6))+1))
        offsets=np.stack(np.meshgrid(xx,yy),axis=-1).reshape(-1,2)
        result=[]
        for tree,c in zip(self.trees,centers):
            query=c[:,None,:]+offsets[None,:,:]
            distance=tree.query(query.reshape(-1,2))[0].reshape(len(c),-1)
            result.append(distance.min(1)-3.*np.sqrt(2))
        return np.asarray(result)

    def unoccupied(self,centers,width,height,pad):
        valid=np.ones(centers.shape[:2],dtype=bool)
        for other in self.tracks:
            for j in range(other['centers'].shape[1]):
                delta=np.abs(centers-other['centers'][:,j,None,:])
                overlap=((delta[:,:,0]<(width+other['width'])/2+pad)&
                         (delta[:,:,1]<(height+other['height'])/2+pad)&
                         (other['alpha'][:,j,None]>1e-8))
                valid &= ~overlap
        return valid

    def add(self,text,size,centers,valid,cost,kind,**extra):
        safe=minimum_filter1d(valid.astype(np.uint8),size=2*self.FADE_RADIUS+1,axis=0,mode='wrap').astype(bool)
        if not safe.any(1).all():
            missing=np.flatnonzero(~safe.any(1))
            raise AssertionError(f'No continuous clear label location for {text}: {missing.tolist()}')
        # Penalize changing a label location; changes that remain necessary are
        # crossfaded, with no position interpolation through the vessel.
        costs=np.where(safe,cost,np.inf)
        periodic_costs=np.tile(costs,(3,1))
        previous=periodic_costs[0].copy()
        history=np.zeros(periodic_costs.shape,dtype=int)
        for i in range(1,len(periodic_costs)):
            best=int(previous.argmin())
            stay=previous
            switch=previous[best]+25.
            use_stay=stay<=switch
            history[i]=np.where(use_stay,np.arange(len(stay)),best)
            previous=periodic_costs[i]+np.minimum(stay,switch)
        selected=np.zeros(len(periodic_costs),dtype=int)
        selected[-1]=int(previous.argmin())
        for i in range(len(periodic_costs)-1,0,-1):
            selected[i-1]=history[i,selected[i]]
        selected=selected[FRAMES:2*FRAMES]
        raw=np.zeros(centers.shape[:2])
        raw[np.arange(FRAMES),selected]=1.
        alpha=convolve1d(raw,self.kernel,axis=0,mode='wrap')
        assert np.all(alpha[~valid]<1e-12)
        assert np.allclose(alpha.sum(1),1.)
        width,height=self.text_size(text,size)
        used=np.flatnonzero(alpha.max(0)>1e-8)
        track=dict(text=text,font_size=size,kind=kind,centers=centers[:,used],
                   alpha=alpha[:,used],width=width,height=height,**extra)
        if 'clearance' in track:
            track['clearance']=track['clearance'][:,used]
        self.tracks.append(track)
        return track

    def audit(self):
        self.max_alpha_step=max(float(np.abs(np.roll(t['alpha'],-1,axis=0)-t['alpha']).max()) for t in self.tracks)
        self.max_position_step=max(float(np.linalg.norm(np.roll(t['centers'],-1,axis=0)-t['centers'],axis=2).max()) for t in self.tracks)
        assert self.max_alpha_step<=self.kernel.max()+1e-12
        assert self.max_position_step<32.,self.max_position_step
        self.minimum_port_clearance=min((float(t['clearance'][t['alpha']>1e-8].min())
            for t in self.tracks if t['kind']=='port'),default=None)
        self.record=dict(text_background_rectangles=0,layout='Continuous projected candidates with periodic crossfades',
            transition_kernel='Raised cosine; periodic convolution, including last-to-first frame',
            transition_samples=len(self.kernel),transition_duration_s=len(self.kernel)/FPS,
            maximum_alpha_change_per_frame=self.max_alpha_step,
            maximum_candidate_position_step_px=self.max_position_step,
            minimum_port_text_to_projected_surface_clearance_px=self.minimum_port_clearance,
            port_arrow_tip='Original boundary cap center',fixed_text_remains_constant=True,
            full_cycle_planned_before_rendering=True,all_pass=True)

    def export(self,path):
        frames=[]
        for i in range(FRAMES):
            labels=[]
            for number,t in enumerate(self.tracks):
                for j,alpha in enumerate(t['alpha'][i]):
                    center=t['centers'][i,j]
                    row=dict(id=f'{number}:{j}',text=t['text'],kind=t['kind'],
                        center_px=center.tolist(),bbox_px=[float(center[0]-t['width']/2),float(center[1]-t['height']/2),
                        float(center[0]+t['width']/2),float(center[1]+t['height']/2)],opacity=float(alpha))
                    if t['kind']=='port':
                        row['arrow_tip_px']=t['anchor'][i].tolist()
                        row['surface_clearance_px']=float(t['clearance'][i,j])
                    elif t['kind']=='tick':
                        row.update(axis='XYZ'[t['axis']],value_um=t['value'],
                                   tick_center_px=t['anchor'][i].tolist(),
                                   tick_opacity=float(t['mark_alpha'][i]))
                    labels.append(row)
            frames.append(dict(frame=i,labels=labels))
        dump(path,dict(reference_size_px=SIZE,validation=self.record,frames=frames))

    def overlay(self,image,i,still=False):
        scale=image.height/SIZE[1]
        layer=Image.new('RGBA',image.size,(0,0,0,0))
        draw=ImageDraw.Draw(layer)
        for t in self.tracks:
            font=ImageFont.truetype(FONT,round(t['font_size']*scale))
            if t['kind']=='tick':
                anchor=t['anchor'][i];normal=t['normal'][i]
                mark_opacity=round(float(t['mark_alpha'][i])*255)
                draw.line([tuple((anchor-normal*5)*scale),tuple((anchor+normal*5)*scale)],
                          fill=(170,183,202,mark_opacity),width=max(1,round(scale)))
            opacity_values=t['alpha'][i].copy()
            if still:
                if t['kind']=='tick':
                    opacity_values=(opacity_values>=.5).astype(float)
                else:
                    opacity_values[:]=0.
                    opacity_values[int(t['alpha'][i].argmax())]=1.
            for j,alpha in enumerate(opacity_values):
                opacity=round(float(alpha)*255)
                if opacity==0:
                    continue
                center=t['centers'][i,j]
                color=(*( (239,244,252) if t['kind']=='port' else (218,228,242) if t['kind']=='axis_title' else (177,191,211)),opacity)
                if t['kind']=='port':
                    anchor=t['anchor'][i]
                    delta=anchor-center
                    # Start at the edge of the transparent text footprint.
                    fraction=1./max(abs(delta[0])/(t['width']/2+6),abs(delta[1])/(t['height']/2+6))
                    tail=center+delta*fraction
                    unit=(anchor-tail)/np.linalg.norm(anchor-tail)
                    perpendicular=np.array([-unit[1],unit[0]])
                    base=anchor-unit*12.
                    draw.line([tuple(tail*scale),tuple(base*scale)],fill=(212,225,240,opacity),width=max(1,round(1.4*scale)))
                    draw.polygon([tuple(anchor*scale),tuple((base+perpendicular*4.5)*scale),
                                  tuple((base-perpendicular*4.5)*scale)],fill=(227,237,250,opacity))
                draw.text(tuple(center*scale),t['text'],font=font,anchor='mm',fill=color)
        return Image.alpha_composite(image.convert('RGBA'),layer).convert('RGB')


class CoordinateGrid:
    """Graduated physical axes through a reference point, not a bounding cage.

    Three Cartesian reference planes use the original absolute coordinates.
    Tick labels are screen-sized overlays at exact projected physical positions.
    The field of view can therefore remain enlarged without clipped box labels.
    """
    def __init__(self, plotter, bounds, center, kind, size):
        self.size = size
        self.step = 20. if kind=='full_vessel' else 5.
        self.origin = np.round(center/self.step)*self.step
        self.ranges = np.asarray(bounds).reshape(3,2)
        self.ticks = [np.arange(np.ceil(lo/self.step)*self.step, hi+.001, self.step)
                      for lo,hi in self.ranges]
        segments = []
        # Coordinate-plane grids at the reference point, in original µm units.
        for normal in range(3):
            plane = [axis for axis in range(3) if axis!=normal]
            for fixed, varying in [plane,plane[::-1]]:
                for value in self.ticks[fixed]:
                    pair = np.tile(self.origin,(2,1))
                    pair[:,fixed] = value
                    pair[:,varying] = self.ranges[varying]
                    segments.extend(pair)
        grid=plotter.add_lines(np.asarray(segments),color='#374252',width=.7*size[1]/1080,connected=False)
        grid.prop.opacity=.44
        self.axis_ends = []
        for axis in range(3):
            pair = np.tile(self.origin,(2,1))
            pair[:,axis] = self.ranges[axis]
            self.axis_ends.append(pair)
        axes=plotter.add_lines(np.concatenate(self.axis_ends),color='#8492a6',width=1.*size[1]/1080,connected=False)
        axes.prop.opacity=.80



class View:
    def __init__(self, wall, lines, orbit, size=SIZE, kind='full_vessel', count=None, ports=None, annotations=None):
        self.orbit = orbit
        self.kind, self.ports, self.size = kind, ports, size
        self.count = count if count is not None else len(np.unique(lines.cell_data['Line_id']))
        self.p = pv.Plotter(off_screen=True, window_size=size, shape=(1, 2), col_weights=[.90, .10], border=False)
        p = self.p
        p.set_background('black', all_renderers=True)
        p.subplot(0, 0)
        p.renderer.SetViewport(*VIEWPORT)
        p.enable_depth_peeling(number_of_peels=8, occlusion_ratio=0)
        p.render_window.SetMultiSamples(0)
        self.actors = [p.add_mesh(wall, color='#b7c4d1', opacity=.17 if kind=='full_vessel' else .15, smooth_shading=True,
                                ambient=.65, diffuse=.35, show_scalar_bar=False)]
        self.actors.append(p.add_mesh(lines, scalars='Speed_mm_s', cmap=CMAP, clim=(0., 7.5),
            line_width=1.55*size[1]/1080, opacity=.92, lighting=False, show_scalar_bar=False))
        orbit.setup(p)
        self.axes = CoordinateGrid(p,wall.bounds,orbit.center,kind,size)
        self.annotations = annotations or AnnotationPlan(orbit,wall,self.axes,ports)
        p.show(auto_close=False, interactive=False)

    def frame(self, i, still=False):
        self.p.camera.position = self.orbit.trace[i]['camera_position_um']
        self.p.camera.up = self.orbit.trace[i]['camera_up_unit']
        self.p.render()
        assert np.allclose(self.p.camera.position, self.orbit.trace[i]['camera_position_um'])
        assert np.allclose(self.p.camera.focal_point, self.orbit.center)
        assert self.p.camera.parallel_scale == self.orbit.scale
        frame=annotate(self.p.screenshot(),self.count,self.kind)
        return self.annotations.overlay(frame,i,still=still)

    def close(self):
        self.p.close()


def validate_video(path, key):
    reader = imageio_ffmpeg.read_frames(str(path))
    metadata = next(reader)
    assert metadata['size'] == SIZE and metadata['fps'] == FPS
    assert abs(metadata['duration'] - FRAMES/FPS) < .05
    frames, distinct = {}, set()
    count = 0
    for i, raw in enumerate(reader):
        count += 1
        distinct.add(hashlib.sha256(raw).hexdigest())
        if i in KEYS:
            frame = np.frombuffer(raw, np.uint8).reshape(SIZE[1], SIZE[0], 3).copy()
            frames[i] = frame
            Image.fromarray(frame).save(OUT / 'inspection' / f'{key}_{i:04d}.png')
            # No obsolete wording beneath the scalar bar.
            assert frame[900:1000, 1735:1910].mean() < 1.
            assert frame[1020:1078, 470:1910].mean() < 1.
    assert count == FRAMES and len(distinct) == FRAMES
    sheet = Image.new('RGB', (3840, 2160), 'black')
    for k, i in enumerate(KEYS[:4]):
        sheet.paste(Image.fromarray(frames[i]), (k % 2 * 1920, k // 2 * 1080))
    sheet.save(OUT / 'figures' / f'{key}_rotation_views.png')
    return dict(file=str(path.relative_to(OUT)), sha256=sha(path), decoded_frames=count,
                distinct_frames=len(distinct), fps=FPS, size=SIZE, duration_s=FRAMES/FPS,
                footer_only_primitive_count=True, no_extra_colorbar_text=True, all_pass=True)


def main():
    global CASE, OUT
    parser = argparse.ArgumentParser()
    parser.add_argument('--stills-only', action='store_true')
    parser.add_argument('--case', type=Path, default=CASE, help='Input case directory (default: bundled input_data)')
    parser.add_argument('--output', type=Path, default=OUT, help='Output directory (default: results)')
    args = parser.parse_args()
    CASE = args.case.expanduser().resolve()
    OUT = args.output.expanduser().resolve()
    start = time.time()
    for folder in ['animations', 'figures', 'inspection', 'data']:
        (OUT / folder).mkdir(parents=True, exist_ok=True)
    sources = ['streamlines/data/streamlines_si.vtp', 'field_diagnostics/data/wall_wss_si.vtp',
               'field_diagnostics/data/pressure_surface_si.vtp', 'frozen_flow/steady_flow_mean_2p0_mmps.vtu']
    hashes = {name: sha(CASE/name) for name in sources}
    compute = json.loads((CASE/'streamlines/COMPUTE_VALIDATION.json').read_text())
    assert compute['all_pass']
    assert hashes[sources[0]] == compute['outputs_sha256']['data/streamlines_si.vtp']
    assert hashes[sources[-1]] == compute['source_field_sha256']
    wall = pv.read(CASE/sources[1])
    full_surface = pv.read(CASE/sources[2])
    lines = pv.read(CASE/sources[0])
    assert lines.n_cells == compute['selected_count'] == 96
    original_speed = lines['Speed_mm_s'].copy()
    for mesh in [wall, full_surface, lines]:
        mesh.points = mesh.points.astype(float) * 1e6
    assert np.array_equal(original_speed, lines['Speed_mm_s'])
    previous = json.loads((CASE/'field_diagnostics/MEDIA_VALIDATION.json').read_text())
    axis = np.array([0.,0.,1.])
    center = fitted_z_center(np.vstack([full_surface.points, lines.points]))
    focus = np.array(json.loads((CASE/'reports/render_manifest.json').read_text())['branch_focus_um'])
    detail_wall = crop_ellipsoid(wall, focus, DETAIL_RADII_UM)
    detail_wall.save(OUT/'data/detail_wall_display_um.vtp')
    field = pv.read(CASE/sources[-1])
    glyphs, glyph_count, glyph_record = velocity_glyphs(field, full_surface, focus)
    ports = []
    for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        cap = pv.read(CASE/'SV_MESH/mesh-surfaces'/f'{role}.vtp')
        cap_center = np.array(cap.center)*1e6
        cap_points = cap.points*1e6
        _,_,vt = np.linalg.svd(cap_points-cap_points.mean(0),full_matrices=False)
        normal = vt[-1]
        if normal@(cap_center-center)<0:
            normal *= -1
        ports.append(('Inlet' if role=='INLET' else 'Outlet '+role[-2:], cap_center, cap_center+normal*6.))
    records, videos = {}, []
    for key, surface, selected, points, pivot, count in [
        ('full_vessel', wall, lines, full_surface.points, center, 96),
        ('branch_vectors', detail_wall, glyphs, detail_wall.points,
         fitted_z_center(np.vstack([detail_wall.points,glyphs.points])), glyph_count),
    ]:
        # Include lines as well as wall in the exact display-framing audit.
        orbit = Turntable(np.vstack([points, selected.points]), axis, pivot)
        dump(OUT/f'{key}_camera.json', orbit.trace)
        record = dict(primitive_count=count, primitive='streamlines' if key=='full_vessel' else 'vectors',
            fixed_parallel_scale_um=orbit.scale, axis_unit=axis.tolist(), center_um=pivot.tolist(),
            pixels_per_um=SIZE[1]*VIEW_HEIGHT/(2*orbit.scale),
            magnification_vs_previous_full=previous['fixed_parallel_scale_um']/orbit.scale*VIEW_HEIGHT/.835,
            magnification_vs_previous_velocity_detail=24./orbit.scale*VIEW_HEIGHT/.835 if key=='branch_vectors' else None,
            all_visible_geometry_inside_view_all_frames=True,
            display_crop_radii_xyz_um=DETAIL_RADII_UM.tolist() if key=='branch_vectors' else None,
            physical_axes_titles=['X (µm)','Y (µm)','Z (µm)'], grid_and_ticks=True,
            z_axis_vertical_all_frames=True, camera_elevation_deg=ELEVATION_DEG)
        if key=='branch_vectors':
            record['glyphs'] = glyph_record
        else:
            record['outlet_line_counts'] = dict(Counter(int(v) for v in lines.cell_data['Outlet_id']))
            record['boundary_labels'] = [p[0] for p in ports]
        print(key, json.dumps(record), flush=True)
        view = View(surface, selected, orbit, size=(3840, 2160), kind=key, count=count, ports=ports if key=='full_vessel' else None)
        annotations=view.annotations
        annotations.export(OUT/f'{key}_annotations.json')
        record['annotations']=annotations.record
        print(key, 'annotations', json.dumps(annotations.record), flush=True)
        view.frame(0,still=True).save(OUT/'figures'/f'{key}_4k.png')
        record['opengl'] = [l for l in view.p.render_window.ReportCapabilities().splitlines()
                            if any(t in l for t in ['OpenGL vendor', 'OpenGL renderer', 'OpenGL version'])]
        view.close()
        view = View(surface, selected, orbit, kind=key, count=count, ports=ports if key=='full_vessel' else None, annotations=annotations)
        record['coordinate_grid_reference_um']=view.axes.origin.tolist()
        record['coordinate_tick_step_um']=view.axes.step
        if args.stills_only:
            for i in KEYS[:4]:
                view.frame(i).save(OUT/'inspection'/f'{key}_preview_{i:04d}.png')
            view.close()
        else:
            path = OUT/'animations'/('streamlines_full_vessel_enlarged.mp4' if key=='full_vessel' else 'velocity_vectors_branch_detail_enlarged.mp4')
            temp = path.with_name(path.stem + '.tmp.mp4')
            writer = imageio_ffmpeg.write_frames(str(temp), size=SIZE, fps=FPS, codec='libx264',
                quality=8, pix_fmt_in='rgb24', pix_fmt_out='yuv420p', macro_block_size=8,
                output_params=['-preset', 'medium', '-movflags', '+faststart'])
            writer.send(None)
            try:
                for i in range(FRAMES):
                    writer.send(np.ascontiguousarray(view.frame(i)))
                    if i % 72 == 0:
                        print(key, i, '/', FRAMES, flush=True)
            finally:
                writer.close()
                view.close()
            os.replace(temp, path)
            videos.append(validate_video(path, key))
        records[key] = record
    assert {name: sha(CASE/name) for name in sources} == hashes
    dump(OUT/('PREVIEW_VALIDATION.json' if args.stills_only else 'MEDIA_VALIDATION.json'),
         dict(all_pass=True, source_sha256=hashes, source_unchanged=True, new_streamline_integrations=0,
              title=TITLE, detail_title='Local detail | FEM velocity vectors', footer_format='{count} streamlines / {count} vectors', colorbar_title='Speed (mm/s)',
              colorbar_range=[0., 7.5], colorbar_extra_text=False, background='black',
              text_background='Transparent; no opaque rectangles behind any text',
              text_transitions='Periodic raised-cosine crossfades, 0.875 s; stable titles and numeric colorbar',
              still_annotation_policy='Dominant clear label position at full opacity, for clean static figures',
              camera='Orthographic camera orbit about one fixed axis and center; constant scale and angular speed',
              full_axis_center_method='Constant-scale fitting for orbit around a fixed axis parallel to original Z through the view center',
              z_axis_vertical_all_frames=True, camera_up=[0,0,1],
              animation_role='360-degree display rotation of a steady flow field',
              source_coordinates_modified=False, scalar_values_modified=False,
              detail='Laterally widened ellipsoidal display crop around the existing velocity-detail focus; clipping interpolates boundary points only',
              frames=FRAMES, fps=FPS, views=records, videos=videos,
              script_sha256=sha(Path(__file__)), elapsed_seconds=time.time()-start))
    if not args.stills_only:
        (OUT/'OPEN_RESULTS.html').write_text('''<!doctype html>
<html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>FEM 全血管流线与局部速度矢量</title>
<style>body{max-width:1280px;margin:32px auto;padding:0 20px;background:#0b0d11;color:#edf2f7;font:17px/1.6 system-ui}
video,img{display:block;width:100%;background:black}a{color:#90caff}section{margin:32px 0}h1{font-size:27px}h2{font-size:21px}</style>
<h1>FEM 全血管流线与局部速度矢量</h1>
<p>透明文字 · 端口箭头标注 · 文字换位渐变 · Z 轴竖直，绕 Z 轴旋转 · 1920 × 1080 · 24 fps · 每段 18 秒</p>
<section><h2>全血管</h2><video controls loop preload="metadata" poster="figures/full_vessel_4k.png"
src="animations/streamlines_full_vessel_enlarged.mp4"></video>
<p><a href="animations/streamlines_full_vessel_enlarged.mp4" download>下载 MP4</a> ·
<a href="figures/full_vessel_4k.png" download>下载 4K 静态图</a></p></section>
<section><h2>局部分叉：速度矢量 / Glyph</h2><video controls loop preload="metadata" poster="figures/branch_vectors_4k.png"
src="animations/velocity_vectors_branch_detail_enlarged.mp4"></video>
<p><a href="animations/velocity_vectors_branch_detail_enlarged.mp4" download>下载 MP4</a> ·
<a href="figures/branch_vectors_4k.png" download>下载 4K 静态图</a></p></section>
<p>局部箭头朝向表示 FEM 速度方向，颜色表示速度大小；箭头长度统一。刻度对应原始物理坐标，单位 µm。</p>
</html>''', encoding='utf-8')
        (OUT/'README_ZH.md').write_text(f'''# Z 轴竖直的全血管流线与局部速度矢量

打开 OPEN_RESULTS.html 可直接预览两段视频，并下载 MP4 和 4K 静态图。

- 全血管：animations/streamlines_full_vessel_enlarged.mp4；相对于上一版固定轴全血管流线，单位物理长度的显示比例放大 {records['full_vessel']['magnification_vs_previous_full']:.3f} 倍。
- 局部分叉：animations/velocity_vectors_branch_detail_enlarged.mp4；相对于之前速度局部图，显示比例放大 {records['branch_vectors']['magnification_vs_previous_velocity_detail']:.3f} 倍。围绕同一焦点的显示区域向 X/Y 两侧扩宽，椭球半轴分别为 18、18、11.5 µm；横向宽度由上一预览的 27 µm 扩至 36 µm，竖向取景聚焦分叉。
- 两段动画均为 1920 × 1080、24 fps、432 帧、18 秒（展示旋转时长）。相机绕固定的竖直 Z 方向旋转，旋转中心和缩放比例恒定；6° 仰角恒定，Z 轴在画面中始终竖直。
- 全血管顶部为 `{TITLE}`，底部仅为 `96 streamlines`；局部顶部为 `Local detail | FEM velocity vectors`，底部仅为 `{glyph_count} vectors`。颜色条仅包含 `Speed (mm/s)` 与 0–7.5 的数字刻度。
- 坐标轴为原始物理坐标 X、Y、Z（单位 µm），不是重新归一化的显示坐标。三个带刻度的参考坐标轴通过记录中的 reference 点；参考面网格与几何共用物理坐标，网格不表示血管壁或流场网格单元。全血管主刻度间隔 20 µm，局部 5 µm。
- 全血管保留 96 条已有流线，三个出口分别有 16、56、24 条，并标注 Inlet / Outlet 01 / Outlet 02 / Outlet 03。数量是展示分配，不等于流量比例。
- 局部 {glyph_count} 个箭头来自四面体重心处的 P1 速度：四个节点速度矢量的平均。采用最远点空间抽样，最小采样间距 1.25 µm；不按速度大小重新指定方向。颜色为速度模长，单位 mm/s；箭头长度统一为 1.55 µm，仅用于表达方向。已检查箭头几何处于原血管内，实际采样数据已导出 CSV/VTP。
- 端口文字与血管投影至少保留 22 px 间距，以箭头指向原始边界中心；所有文字均为透明背景。需要换位或隐现的文字使用 0.875 秒周期交叉淡化，包含视频末帧回到首帧的衔接；标题、颜色条等固定文字保持不变。两个 annotations.json 保存逐帧位置、透明度与检查记录。
- 沿用已有流线数据，无新增积分；完整流场和原始流线文件均通过 SHA256 不变检查。
- 全血管原始几何及局部裁剪后的几何均完成逐帧视野边界检查；两段视频均完成全帧解码与重复帧检查。

MEDIA_VALIDATION.json 和两个 camera.json 保存验证与相机记录；data 下为局部显示用几何，坐标单位 µm。
''', encoding='utf-8')
        files = sorted(p for p in OUT.rglob('*') if p.is_file() and p.name != 'SHA256SUMS.txt')
        (OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(OUT)}\n' for p in files))
    print('PREVIEWS_COMPLETE' if args.stills_only else 'RENDER_COMPLETE', flush=True)


if __name__ == '__main__':
    main()
