"""Sampled geometric access, separate from slicer support removability.

FCL checks capsule probes against the actual accepted vascular surface and frame.
VTK signed distance plus a deterministic grid flood-fill detects isolated sampled
free-space pockets. Neither approximation certifies removal of real supports.
"""
import numpy as np
from scipy import ndimage
import trimesh
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from . import sacrificial_fixture as legacy
from .surface_continuity_qc import surface


def vascular_regions(inputs,cfg):
    base=inputs['base'];previous=inputs['folder']/'resolved_config.yaml'
    settings=legacy.load_config(previous)['input']
    graph=legacy.read_source(base/settings['fitted_swc']).graph
    mapping=legacy.read_csv(base/settings['edge_mapping'])
    transform=np.array(__import__('json').loads((base/settings['transform']).read_text())['transform_4x4'])
    groups={}
    for row in mapping:
        a,b=int(row['parent_id']),int(row['child_id'])
        groups.setdefault(row['original_branch'],[]).append((a,b))
    regions=[]
    for branch,edges in sorted(groups.items()):
        # The frozen mapping supplies identity; no new branch or skeleton model.
        samples=[]
        for a,b in edges:
            ca,cb=graph.nodes[a]['coords'],graph.nodes[b]['coords']
            length=np.linalg.norm(cb[:3]-ca[:3]);n=max(1,int(np.ceil(length/2.)))
            samples.extend(ca*(1-u)+cb*u for u in np.linspace(0,1,n,endpoint=False))
        samples=np.array(samples)
        distances=np.r_[0,np.cumsum(np.linalg.norm(np.diff(samples[:,:3],axis=0),axis=1))]
        bins=np.floor(distances/cfg['support_access']['region_length_mm']).astype(int)
        for group in np.unique(bins):
            subset=samples[bins==group];point=subset[len(subset)//2]
            regions.append(dict(region_id=f'{branch}_{group:02}',nearest_vascular_branch=branch,
                center_mm=transform[:3,:3]@point[:3]+transform[:3,3],radius_mm=float(point[3])))
    return regions


def capsule(start,end,radius,sections):
    vector=np.asarray(end)-start;length=np.linalg.norm(vector)
    if length<1e-8:raise ValueError('Zero-length access probe')
    shape=trimesh.creation.capsule(height=length,radius=radius,count=[sections,sections])
    # trimesh capsule is centered at the origin, with cap centers at +/-h/2.
    transform=trimesh.geometry.align_vectors([0,0,1],vector/length)
    transform[:3,3]=(np.asarray(start)+end)/2
    shape.apply_transform(transform)
    return shape


def line_visible(frame,core,start,end):
    direction=np.asarray(end)-start;length=np.linalg.norm(direction)
    if length<1e-8:return False
    direction/=length
    for mesh in (frame,core):
        hits,_,_=mesh.ray.intersects_location(np.array([start]),np.array([direction]),multiple_hits=True)
        if len(hits):
            distances=(hits-start)@direction
            if np.any((distances>1e-5)&(distances<length-1e-5)):return False
    return True


def free_space_components(design,core,cfg):
    bounds=design['outer_bounds'];step=cfg['support_access']['free_space_grid_mm']
    dims=np.ceil((bounds[1]-bounds[0])/step).astype(int)+1
    coordinates=[np.linspace(bounds[0,k],bounds[1,k],dims[k]) for k in range(3)]
    points=np.stack(np.meshgrid(*coordinates,indexing='ij'),axis=-1)
    probe=cfg['support_removal_probe']['diameter_mm']/2
    occupied=np.zeros(tuple(dims),bool)
    for member in design['members']:
        lo,hi=member['bounds'];gap=np.maximum(np.maximum(lo-points,points-hi),0.)
        occupied |= np.linalg.norm(gap,axis=-1)<=probe+1e-8
    distance=vtk.vtkImplicitPolyDataDistance();distance.SetInput(surface(core))
    sample=vtk.vtkSampleFunction();sample.SetImplicitFunction(distance)
    sample.SetModelBounds(*bounds.T.ravel());sample.SetSampleDimensions(*[int(n) for n in dims]);sample.ComputeNormalsOff();sample.Update()
    distances=vtk_to_numpy(sample.GetOutput().GetPointData().GetScalars()).reshape(tuple(dims),order='F')
    occupied|=distances<=probe
    labels,count=ndimage.label(~occupied,ndimage.generate_binary_structure(3,1))
    axis=design['open_axis'];exit_labels=set()
    for index in (0,-1):
        selection=[slice(None)]*3;selection[axis]=index
        exit_labels.update(np.unique(labels[tuple(selection)]).tolist())
    exit_labels.discard(0)
    return dict(coordinates=coordinates,labels=labels,exit_labels=exit_labels,
        component_count=count,grid_spacing_mm=[float(x[1]-x[0]) for x in coordinates],
        occupied_fraction=float(occupied.mean()))


def point_component(point,grid):
    indices=tuple(int(np.argmin(abs(axis-point[k]))) for k,axis in enumerate(grid['coordinates']))
    label=int(grid['labels'][indices])
    if label==0:
        # A valid probe endpoint can round onto an occupied sample at a curved
        # vessel boundary. Search only within half the cell diagonal, preserving
        # the stated grid uncertainty; do not jump through a whole wall.
        axes=grid['coordinates'];limit=np.linalg.norm([x[1]-x[0] for x in axes])/2
        nearby=[]
        import itertools
        for cell in itertools.product(*[range(max(0,i-1),min(len(axes[k]),i+2)) for k,i in enumerate(indices)]):
            value=int(grid['labels'][cell])
            distance=np.linalg.norm(np.array([axes[k][cell[k]] for k in range(3)])-point)
            if value and distance<=limit+1e-8:nearby.append((distance,value))
        if nearby:label=min(nearby)[1]
    return label,label!=0 and label in grid['exit_labels']


def evaluate_access(design,core,regions,cfg):
    manager=trimesh.collision.CollisionManager();manager.add_object('frame',design['mesh']);manager.add_object('vascular',core)
    grid=free_space_components(design,core,cfg);axis=design['open_axis'];bounds=design['outer_bounds']
    radius=cfg['support_removal_probe']['diameter_mm']/2;standoff=cfg['support_removal_probe']['surface_standoff_mm']
    rows=[];paths=[]
    for region in regions:
        center=np.asarray(region['center_mm']);r=region['radius_mm'];row=dict(region_id=region['region_id'],nearest_vascular_branch=region['nearest_vascular_branch'])
        connectivity=[];labels=[]
        for side,letter in [(0,'A'),(1,'B')]:
            sign=-1 if side==0 else 1
            best=None;ray_trials=[];other=[k for k in range(3) if k!=axis]
            offsets=cfg['support_access']['entry_offsets_mm']
            for offset0 in offsets:
                for offset1 in offsets:
                    start=center.copy();start[axis]=bounds[side,axis]+sign*(radius+standoff)
                    for k,offset in zip(other,(offset0,offset1)):
                        start[k]=np.clip(center[k]+offset,design['inner_bounds'][0,k]+radius,design['inner_bounds'][1,k]-radius)
                    heading=start-center;heading/=np.linalg.norm(heading)
                    hits,_,_=core.ray.intersects_location(center[None,:],heading[None,:],multiple_hits=True)
                    distances=(hits-center)@heading if len(hits) else np.array([])
                    distances=distances[distances>1e-5]
                    if not len(distances):continue
                    surface_offset=float(distances.min())
                    if surface_offset>max(3*r,r+1.):continue
                    surface_point=center+surface_offset*heading
                    end=surface_point+(radius+standoff)*heading
                    visible_end=surface_point+standoff*heading
                    visible=line_visible(design['mesh'],core,start,visible_end)
                    collision=True
                    if visible:
                        shape=capsule(start,end,radius,cfg['support_removal_probe']['capsule_sections'])
                        collision=manager.in_collision_single(shape)
                    label,reachable=point_component(end,grid)
                    candidate=dict(region_id=region['region_id'],face=design['open_faces'][side],start=start.tolist(),
                        end=end.tolist(),probe_accessible=not collision,visible=visible,
                        grid_label=label,grid_reachable=reachable,entry_offsets_mm=[offset0,offset1])
                    ray_trials.append(candidate)
                    if best is None or (candidate['probe_accessible'],visible,reachable)>(best['probe_accessible'],best['visible'],best['grid_reachable']):best=candidate
            if best is None:
                best=dict(region_id=region['region_id'],face=design['open_faces'][side],start=center.tolist(),end=center.tolist(),
                    probe_accessible=False,visible=False,grid_label=0,grid_reachable=False,entry_offsets_mm=[])
            row['visible_from_open_face_'+letter]=bool(any(v['visible'] for v in ray_trials))
            row['probe_accessible_'+letter]=bool(any(v['probe_accessible'] for v in ray_trials))
            row['tested_paths_'+letter]=len(ray_trials)
            labels.append(best['grid_label']);connectivity.append(any(v['grid_reachable'] for v in ray_trials))
            paths.append(best)
        los=row['visible_from_open_face_A'] or row['visible_from_open_face_B']
        probe=row['probe_accessible_A'] or row['probe_accessible_B']
        # Only a nonzero isolated free-space label supports a trapped-pocket
        # proxy. A blocked voxel or missing straight ray alone does not prove it.
        trapped=bool(any(labels) and not any(connectivity) and not los and not probe)
        row.update(grid_exit_reachable=bool(any(connectivity)),trapped_support_proxy=trapped,
            status='TRAPPED_SUPPORT_PROXY' if trapped else ('GEOMETRY_SUPPORT_ACCESS_PASS' if los or probe else 'GEOMETRY_SUPPORT_ACCESS_RISK'))
        rows.append(row)
    count=len(rows);accessible=sum(r['status']=='GEOMETRY_SUPPORT_ACCESS_PASS' for r in rows)
    probes=sum(r['probe_accessible_A'] or r['probe_accessible_B'] for r in rows)
    trapped=sum(r['trapped_support_proxy'] for r in rows)
    result=dict(status='GEOMETRY_SUPPORT_ACCESS_PASS' if accessible==count and not trapped else 'GEOMETRY_SUPPORT_ACCESS_RISK',
        evidence_level='GEOMETRIC_SUPPORT_ACCESS_ONLY',region_count=count,accessible_region_count=accessible,
        inaccessible_region_count=count-accessible,probe_accessible_region_count=probes,
        accessible_fraction=accessible/count if count else 0.,probe_accessible_fraction=probes/count if count else 0.,
        trapped_proxy_region_count=trapped,grid_spacing_mm=grid['grid_spacing_mm'],free_space_component_count=grid['component_count'],
        actual_slicer_support_removal_verified=False,
        limitation='Sampled sightlines, capsule probes and voxel free-space connectivity are accessibility proxies; actual support break-up and removal are not certified.')
    return result,rows,paths
