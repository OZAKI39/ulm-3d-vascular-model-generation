"""Source-sized blind magnet cups on opposite exterior sides of casting molds.

Read the actual trimmed Brep edges in the supplied Rhino model. Its four
touching solids form a rectangular pad and circular blind pocket; directly
unioning their independent render meshes creates microscopic seam tunnels.
Consolidate that same analytic shape before meshing, without scaling the slot.
Construction curves and points are retained in the reference 3dm, not printed.
"""
from pathlib import Path
import shutil

import cadquery as cq
import numpy as np
import rhino3dm as rhino
import trimesh

from . import y_channel_casting as y
from . import straight_channel_casting as s


def validate_config(cfg):
    p=cfg['magnet_slots']
    if cfg.get('fixture_mode')!='bottomless_external_dovetail':
        raise ValueError('MAGNETS_REQUIRE_BOTTOMLESS_FIXTURE')
    if p['sides']!=['Y_MINUS','Y_PLUS']:
        raise ValueError('MAGNETS_REQUIRE_OPPOSITE_Y_SIDES')
    if not .1<=p['wall_overlap_mm']<=.4 or not 0.5<=p['fixture_clearance_mm']<=3:
        raise ValueError('INVALID_MAGNET_ATTACHMENT_CLEARANCE')
    if y.digest(y.ROOT/p['source'])!=p['source_sha256']:
        raise ValueError('MAGNET_SOURCE_SHA256_MISMATCH')


def read_source(cfg):
    validate_config(cfg)
    path=y.ROOT/cfg['magnet_slots']['source'];doc=rhino.File3dm.Read(str(path))
    if doc is None or doc.Settings.ModelUnitSystem!=rhino.UnitSystem.Millimeters:
        raise ValueError('MAGNET_SOURCE_REQUIRES_MM')
    solids=[o.Geometry for o in doc.Objects if isinstance(o.Geometry,rhino.Brep)]
    if len(solids)!=4 or not all(b.IsSolid and b.IsValid for b in solids):
        raise ValueError('UNSUPPORTED_MAGNET_SOURCE_SOLIDS')
    # Trimmed edge/vertex coordinates, never the untrimmed surface bounding box.
    vertices=np.array([[v.Location.X,v.Location.Y,v.Location.Z] for b in solids for v in b.Vertices])
    lo,hi=vertices.min(0),vertices.max(0)
    circles=[c for b in solids for e in b.Edges if (c:=e.TryGetCircle()) is not None]
    radii=np.array([c.Radius for c in circles]);radius=float(radii.min())
    inner=[c for c in circles if abs(c.Radius-radius)<1e-5]
    centers=np.array([[c.Center.X,c.Center.Y,c.Center.Z] for c in inner])
    center=centers[0,:2];levels=np.unique(np.round(centers[:,2]-lo[2],6))
    if len(levels)!=3 or np.max(np.linalg.norm(centers[:,:2]-center,axis=1))>1e-5:
        raise ValueError('UNSUPPORTED_MAGNET_BLIND_POCKET')
    size=hi-lo;floor=float(levels[1]);height=float(size[2])
    if not (np.allclose(size,[8.,8.60641266,3.],atol=2e-5) and abs(radius-3.2)<1e-5 and abs(floor-1.)<1e-5):
        raise ValueError('MAGNET_SOURCE_DIMENSIONS_CHANGED')
    cxy=center-lo[:2]
    # These source shapes are straight extrusions; verify intermediate edge
    # samples remain in the rectangular envelope and the circular recess.
    for b in solids:
        for e in b.Edges:
            pts=np.array([[q.X,q.Y,q.Z] for q in
                          (e.PointAt(t) for t in np.linspace(e.Domain.T0,e.Domain.T1,257))])
            if ((pts<lo-1e-5)|(pts>hi+1e-5)).any():
                raise ValueError('MAGNET_NONRECTANGULAR_OUTER_PROFILE')
            if ((np.linalg.norm(pts[:,:2]-center,axis=1)<radius-1e-5)&(pts[:,2]>lo[2]+floor+1e-5)).any():
                raise ValueError('MAGNET_SOURCE_OCCUPIES_POCKET')
    solid=cq.Solid.makeBox(*map(float,size)).cut(cq.Solid.makeCylinder(radius,height-floor+.1,
        cq.Vector(*map(float,cxy),floor),cq.Vector(0,0,1))).clean()
    mesh=s.cad_mesh(solid,cfg)
    ideal=float(np.prod(size)-np.pi*radius**2*(height-floor))
    if not mesh.is_volume or mesh.euler_number!=2 or abs(mesh.volume-ideal)>.03:
        raise ValueError('MAGNET_CLEAN_SOLID_INVALID')
    meta=dict(source=cfg['magnet_slots']['source'],source_sha256=y.digest(path),units='mm',
        source_solid_count=4,ignored_construction_object_count=len(doc.Objects)-4,
        source_bounds_mm=[lo,hi],pad_size_mm=size,pocket_center_local_xy_mm=cxy,
        pocket_radius_mm=radius,pocket_diameter_mm=2*radius,pocket_floor_mm=floor,
        pocket_depth_mm=height-floor,source_dimensions_rescaled=False,
        consolidation='Trimmed Brep edges define one rectangular pad minus the original circular blind pocket; mesh seams removed',
        analytic_volume_mm3=ideal,mesh_volume_mm3=float(mesh.volume),single_volume=True)
    return mesh,meta


def transform_for(side,x,outer,meta,overlap):
    """Rigidly rotate the pad; no mirror/scale of the supplied pocket."""
    sign=1. if side=='Y_PLUS' else -1.
    T=np.diag([sign,sign,1.,1.])
    T[0,3]=x-sign*meta['pad_size_mm'][0]/2
    T[1,3]=outer[1 if sign>0 else 0,1]-sign*overlap
    return T


def moved(m,T):
    result=m.copy();result.apply_transform(T);return result


def add_slots(cfg,frame,other_parts,inner,out):
    """Choose a common X for opposing cups close to the chamber midpoint."""
    slot,meta=read_source(cfg);p=cfg['magnet_slots']
    outer=np.asarray(inner)+[[-cfg['wall_mm']]*2,[cfg['wall_mm']]*2]
    half=meta['pad_size_mm'][0]/2;clearance=p['fixture_clearance_mm']
    midpoint=float(inner[:,0].mean())
    candidates=np.arange(outer[0,0]+half+1,outer[1,0]-half-1+.001,.25)
    candidates=sorted(set([midpoint,*map(float,candidates)]),key=lambda x:(abs(x-midpoint),x))
    # Exclude the deliberate attachment strip from the baseline-frame test.
    wall_box=y.box([np.r_[outer[0],-.1],np.r_[outer[1],30]])
    external_fixture=trimesh.boolean.difference([frame,wall_box],engine='manifold')
    obstacles=[external_fixture,*other_parts.values()]
    chosen=None
    for x in candidates:
        placed=[];records=[]
        for side in p['sides']:
            T=transform_for(side,x,outer,meta,p['wall_overlap_mm']);m=moved(slot,T)
            # Keep both cups away from the external sockets and retainers,
            # including overhead access. Contact with the plain wall is desired.
            b=m.bounds.copy();b[:,:2]+=[[-clearance]*2,[clearance]*2];b[1,2]=22
            clear=y.box(b)
            if any(s.volume_intersection(ob,clear)>1e-4 for ob in obstacles):break
            placed.append(m)
            center=(T@np.r_[meta['pocket_center_local_xy_mm'],meta['pocket_floor_mm'],1])[:3]
            records.append(dict(side=side,slot_to_casting_transform=T,pocket_floor_center_mm=center,
                pocket_axis=[0,0,1],slot_bounds_mm=m.bounds,attachment_overlap_mm=p['wall_overlap_mm'],
                clearance_to_external_fixtures_mm=clearance))
        if len(placed)==2:
            chosen=(placed,records);break
    if chosen is None:raise ValueError('NO_CLEAR_OPPOSING_MAGNET_SLOT_POSITION')
    placed,records=chosen;result=y.union([frame,*placed])
    qc=check_slots(result,frame,other_parts,inner,cfg,meta,records)
    if not qc['passed']:raise ValueError('MAGNET_GEOMETRY_QC_FAILED: '+str(qc))
    ref=out/'reference';ref.mkdir(parents=True,exist_ok=True)
    shutil.copy2(y.ROOT/p['source'],ref/'magnet_slot.3dm')
    slot.export(ref/'magnet_slot_consolidated_reference.stl')
    report=dict(enabled=True,count=2,frame_part_only=True,source=meta,placements=records,
        common_x_mm=x,cavity_midpoint_x_mm=midpoint,placement_search_count=len(candidates),
        placement_rule='Opposing Y walls; closest feasible common X to cavity midpoint; >=1 mm fixture access envelope',
        qc=qc,physical_magnetic_retention_tested=False)
    y.save(out/'magnet_slots_qc.json',report)
    return result,report


def check_slots(frame,baseline,others,inner,cfg,meta,records):
    # Read actual mesh ray hits: pocket remains empty above its 1 mm floor;
    # floor is continuous, side walls at the source diameter, access is clear.
    rows=[];r=meta['pocket_radius_mm'];floor=meta['pocket_floor_mm'];height=meta['pad_size_mm'][2]
    for record in records:
        c=np.array(record['pocket_floor_center_mm']);origin=c+[0,0,10]
        hit=frame.ray.intersects_location([origin],[[0,0,-1]])[0]
        levels=np.unique(np.round(hit[:,2],5));levels.sort()
        widths=[]
        for z in [floor+.4,height-.3]:
            widths_z=[]
            for d in [[1,0,0],[0,1,0]]:
                distances=[]
                for sign in [-1,1]:
                    vec=np.array(d)*sign;orig=np.r_[c[:2],z]
                    hits=frame.ray.intersects_location([orig],[vec])[0]
                    positive=(hits-orig)@vec;positive=positive[positive>1e-5]
                    distances.append(float(positive.min()) if len(positive) else float('inf'))
                widths_z.append(sum(distances))
            widths.append(widths_z)
        cavity=y.cylinder(c[:2],r-.01,floor+.01,22)
        blocked=s.volume_intersection(frame,cavity)+sum(s.volume_intersection(o,cavity) for o in others.values())
        solid_floor=y.cylinder(c[:2],r-.01,.01,floor-.01)
        floor_missing=max(0.,float(solid_floor.volume)-s.volume_intersection(frame,solid_floor))
        passed=(len(levels)==2 and np.allclose(levels,[0,floor],atol=1e-4)
                and np.allclose(widths,2*r,atol=.004) and blocked<1e-4 and floor_missing<1e-4)
        rows.append(dict(side=record['side'],vertical_ray_levels_mm=levels,pocket_diameters_xy_mm=widths,
            pocket_access_blocked_volume_mm3=blocked,floor_missing_volume_mm3=floor_missing,passed=bool(passed)))
    original_missing=max(0.,float(baseline.volume)-s.volume_intersection(frame,baseline))
    interior=y.box([np.r_[inner[0]+1e-4,-.1],np.r_[inner[1]-1e-4,22]])
    intrusion=s.volume_intersection(frame,interior)
    gates=dict(two_clear_blind_pockets=len(rows)==2 and all(r['passed'] for r in rows),
        single_integrated_frame=frame.is_volume and len(frame.split())==1,
        no_new_mesh_tunnels=frame.euler_number==baseline.euler_number,
        original_frame_preserved=original_missing<1e-3,no_chamber_intrusion=intrusion<1e-4,
        frame_height_preserved=abs(frame.bounds[1,2]-baseline.bounds[1,2])<1e-5)
    return dict(passed=all(gates.values()),gates=gates,pockets=rows,
        original_frame_missing_volume_mm3=original_missing,internal_abs_volume_mm3=intrusion,
        added_abs_volume_mm3=float(frame.volume-baseline.volume))


def native_toolpath_qc(summary,records,segments):
    """Check actual positive extrusion leaves each blind pocket open."""
    from .bambu_casting_mold import sample_segments
    report=summary['magnet_slots'];meta=report['source'];r=meta['pocket_radius_mm']
    row=next(r for r in records if 'box' in r['name']);T=np.asarray(row['casting_to_plate_transform'])
    points=sample_segments(segments,.15);results=[]
    for slot in report['placements']:
        center=(T@np.r_[slot['pocket_floor_center_mm'],1])[:3]
        radial=np.linalg.norm(points[:,:2]-center[:2],axis=1)
        floor=(radial<r-.6)&(points[:,2]>.05)&(points[:,2]<=center[2]+.05)
        filled=(radial<r-.55)&(points[:,2]>center[2]+.25)&(points[:,2]<3.01)
        wall=(abs(radial-r)<.65)&(points[:,2]>center[2]+.25)&(points[:,2]<3.01)
        layers=np.unique(np.round(points[wall,2],2))
        results.append(dict(side=slot['side'],floor_extrusion_samples=int(floor.sum()),
            undesired_pocket_extrusion_samples=int(filled.sum()),pocket_wall_layer_tops_mm=layers,
            passed=bool(floor.sum()>10 and not filled.any() and len(layers)>=7)))
    return dict(passed=all(r['passed'] for r in results),pockets=results,
        scope='Native positive model extrusion; source-size recess, floor and wall layers; no physical force test')
