"""Layer-friendly ruled curtains, arch windows and explicit geometric cuts.

All dimensions are millimetres. These are geometry hypotheses, not calibrated
ABS breakaway or bridge parameters. No anatomical mesh is modified here.
"""
from __future__ import annotations
import numpy as np
import trimesh
from shapely.geometry import Polygon, GeometryCollection, box as rectangle
from .support_v2_primitives import quantize, loft
from .support_detection import section_polygon


def union(meshes):
    meshes = [m for m in meshes if len(m.faces)]
    if not meshes:
        return trimesh.Trimesh()
    return trimesh.boolean.union(meshes, engine='manifold') if len(meshes) > 1 else meshes[0].copy()


def sweep(xy, bottom, top, widths, lower_widths=None):
    """Closed ruled strip; curved XY spine, horizontal crosswise stations.

    Separate stations avoid a forest of cylinder primitives and retain a
    continuous longitudinal extrusion route. Callers split folded projections.
    """
    xy = np.asarray(xy, float)
    n = len(xy)
    if n < 2:
        raise ValueError('A_RIB_NEEDS_TWO_DISTINCT_STATIONS')
    tangent = np.gradient(xy, axis=0)
    norms = np.linalg.norm(tangent, axis=1)
    if np.any(norms < 1e-8):
        raise ValueError('FOLDED_OR_VERTICAL_RIB_PROJECTION')
    normal = np.c_[-tangent[:, 1], tangent[:, 0]] / norms[:, None]
    widths = np.broadcast_to(widths, (n,))
    lower = widths if lower_widths is None else np.broadcast_to(lower_widths, (n,))
    bottom, top = np.broadcast_arrays(np.asarray(bottom), np.asarray(top), np.zeros(n))[:2]
    if np.any(top <= bottom):
        raise ValueError('NON_POSITIVE_RIB_HEIGHT')
    vertices = []
    for i in range(n):
        for width, z in [(lower[i], bottom[i]), (widths[i], top[i])]:
            vertices.extend([np.r_[xy[i]-normal[i]*width/2, z], np.r_[xy[i]+normal[i]*width/2, z]])
    faces = [[0, 1, 3], [0, 3, 2], [4*n-4, 4*n-2, 4*n-1], [4*n-4, 4*n-1, 4*n-3]]
    for i in range(n-1):
        a, b = 4*i, 4*(i+1)
        for j, k in [(0, 2), (2, 3), (3, 1), (1, 0)]:
            faces.extend([[a+j, b+j, b+k], [a+j, b+k, a+k]])
    m = trimesh.Trimesh(vertices, faces, process=True)
    trimesh.repair.fix_normals(m)
    if not m.is_volume:
        raise ValueError('RIB_SWEEP_NOT_CLOSED')
    return m


def window_cutter(center, tangent, width, bottom, apex, depth, slope=.70):
    """An open-bottom, pointed arch. Sides close at <=0.7 mm/mm.

    The roof crown remains material, not an uncalibrated horizontal bridge.
    Width is the clear rectangular opening below the arch shoulders.
    """
    center = np.asarray(center, float)
    tangent = np.asarray(tangent, float)[:2]
    tangent /= np.linalg.norm(tangent)
    normal = np.array([-tangent[1], tangent[0]])
    shoulder = apex-width/2/slope
    if shoulder <= bottom:
        raise ValueError('NO_HEIGHT_FOR_CLEAR_ARCH_WINDOW')
    polygon = Polygon([[-width/2, bottom], [width/2, bottom], [width/2, shoulder], [0, apex], [-width/2, shoulder]])
    local = trimesh.creation.extrude_polygon(polygon, depth)
    # local x=along rib, y=height, z=across rib
    m = np.eye(4)
    m[:3, 0] = [*tangent, 0]
    m[:3, 1] = [0, 0, 1]
    m[:3, 2] = [*normal, 0]
    m[:3, 3] = [*(center[:2]-normal*depth/2), 0]
    local.apply_transform(m)
    return local, shoulder


def plane_frame(origin, normal):
    m = trimesh.geometry.align_vectors([0, 0, 1], normal)
    m[:3, 3] = origin
    return m


def cut_section(support, intended, origin, normal, aperture, tolerance=1e-5):
    """Audit a finite tool plane, including neighbouring unweakened members.

    The intended material comes from the isolated, pre-union neck primitive;
    it is never redefined to include whatever the completed scaffold contains.
    """
    inverse = np.linalg.inv(plane_frame(origin, normal))
    def section(mesh):
        q = mesh.copy(); q.apply_transform(inverse)
        path = q.section_multiplane([0, 0, 0], [0, 0, 1], [0])[0]
        return section_polygon(path).intersection(aperture)
    actual, allowed = section(support), section(intended)
    extra = actual.difference(allowed.buffer(1e-6))
    area = float(extra.area)
    # Propagate binary-STL float32 coordinate precision into a boundary-area
    # bound. This is recorded per plane, not an arbitrary relaxed percentage.
    ulp=float(np.spacing(np.float32(max(1.,np.max(np.abs(origin))))))
    area_tolerance=max(tolerance,allowed.length*2*ulp+np.pi*(2*ulp)**2)
    return dict(passed=area <= area_tolerance, intended_neck_cross_section_area_mm2=float(allowed.area),
                other_support_cross_section_area_mm2=area, area_tolerance_mm2=area_tolerance,float32_coordinate_ulp_mm=ulp,
                status='GEOMETRIC_CUT_SECTION_VALID' if area <= area_tolerance else 'CUT_SECTION_CROSSES_UNWEAKENED_SUPPORT')


def curtain(xy, top, base, profile, lines=2, roof_width=None, segment_length=5., window_width=3.):
    """Build one interval curtain, segmented by arched tool windows.

    Root notch narrows the wall at the cut height; roof growth uses a <=0.65
    outward slope. All horizontal thicknesses use integer extrusion widths.
    """
    w, h = profile['line_width'], profile['layer_height']
    xy = np.asarray(xy, float); top = np.asarray(top, float)
    base = np.broadcast_to(base, top.shape).copy()
    width = lines*w
    neck_width=min(lines,int(profile.get('cut_neck_lines',1)))*w
    roof_width = quantize(roof_width or (lines+2)*w, w)
    transition = quantize(max(roof_width-width, 0)/2/.60 + 3*h, h)
    roof_bottom = top-transition
    # A full-length, one-line root neck is explicitly separated from the body.
    neck_z = base+5*h
    neck_top = base+7*h
    flare_end = neck_top+quantize(max(width-neck_width, 0)/2/.60+h, h)
    if np.any(roof_bottom <= flare_end+.2):
        raise ValueError('INSUFFICIENT_HEIGHT_FOR_RIB_ROOT_AND_ROOF')
    foot = sweep(xy, np.maximum(base-.03, 0), base+3*h, width+2*w)
    neck = sweep(xy, base+3*h, neck_top, neck_width)
    flare = sweep(xy, neck_top-.02, flare_end, width, neck_width)
    body = sweep(xy, flare_end-.02, roof_bottom+.02, width)
    roof = sweep(xy, roof_bottom, top, roof_width, width)
    assembly = union([foot, neck, flare, body, roof])
    arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
    length = arc[-1]
    windows, cuts, window_meshes = [], [], []
    # Search length is a software heuristic. Curvature and available height
    # are handled by the caller; a short interval needs no artificial window.
    count = int(np.floor((length+window_width)/(segment_length+window_width)))
    centers = np.linspace(0, length, count+1)[1:-1] if count > 1 else []
    for s in centers:
        p = np.array([np.interp(s, arc, xy[:, i]) for i in range(2)])
        i = min(np.searchsorted(arc, s), len(arc)-2)
        direction = xy[i+1]-xy[i]
        if np.linalg.norm(direction) < 1e-6:
            continue
        apex = float(np.interp(s, arc, top)-3*h)
        low = float(np.interp(s, arc, base)-1)
        try:
            cutter, shoulder = window_cutter(p, direction, window_width, low, apex, roof_width*4+2)
        except ValueError:
            continue
        window_meshes.append(cutter)
        windows.append(dict(arc_mm=float(s), xyz=[*p, apex], clear_width_mm=window_width,
                            clear_shoulder_z_mm=shoulder, arch_closure_slope=.70))
        # Finite cross-cut of the thin crown above the arch apex.
        tangent = direction/np.linalg.norm(direction)
        origin = np.r_[p, apex+1.5*h]
        frame = plane_frame(origin, np.r_[tangent, 0])
        cutter2 = trimesh.creation.box([roof_width*3+2, 6*h, .06])
        # A global plane-frame rectangle encompasses the entire local crown.
        cutter2 = trimesh.creation.box([roof_width*3+2, roof_width*3+2, .06], transform=frame)
        cuts.append(dict(kind='WINDOW_CROWN', origin=origin, normal=np.r_[tangent, 0], cutter=cutter2,
                         aperture=rectangle(-roof_width*1.5-1,-roof_width*1.5-1,roof_width*1.5+1,roof_width*1.5+1),
                         aperture_halfwidth=roof_width*1.5+1, probe_target=[*p, shoulder-1.5]))
    if window_meshes:
        holes = union(window_meshes)
        def subtract(m): return trimesh.boolean.difference([m, holes], engine='manifold')
        assembly, foot, neck, body, roof = map(subtract, [assembly, foot, neck, union([flare,body]), roof])
    else:
        body = union([flare, body])
    # Root cuts are local along each contiguous post, not an infinite plane.
    for part in neck.split():
        center = part.bounds.mean(0); z = center[2]
        size = part.extents[:2]+2*w
        cutter = loft(center-[0,0,.03], center+[0,0,.03], size)
        cuts.append(dict(kind='ROOT_NECK', origin=center, normal=[0,0,1], cutter=cutter,
                         aperture=rectangle(-size[0]/2,-size[1]/2,size[0]/2,size[1]/2),
                         aperture_halfwidth=float(max(size)/2), probe_target=center.tolist()))
    retained=[foot]
    for cut in cuts:
        cut['intended'] = trimesh.boolean.intersection([assembly,cut['cutter']], engine='manifold')
        if cut['kind']=='ROOT_NECK':
            lo,hi=cut['cutter'].bounds;low=lo.copy();low[2]=base.min()-1
            below=trimesh.creation.box(hi-low,transform=trimesh.transformations.translation_matrix((hi+low)/2))
            # Include the lower half of the designed neck with its foot.
            below.vertices[:,2]=np.minimum(below.vertices[:,2],cut['origin'][2]-.03)
            retained.append(trimesh.boolean.intersection([neck,below],engine='manifold'))
    return dict(mesh=assembly, ribs=body, roofs=roof, necks=neck, feet=foot, cuts=cuts, windows=windows,
                approved_retained_stubs=union(retained),
                segment_count=len(windows)+1, width_mm=width, roof_width_mm=roof_width,
                segment_length_mm=segment_length, path_xy=xy, top=top, base=base)
