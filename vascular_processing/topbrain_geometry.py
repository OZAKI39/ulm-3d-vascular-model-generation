"""Original TopBrain mask surface in affine world coordinates; no tree conversion."""
import numpy as np
from scipy import ndimage as ndi
from .topbrain_qc import crop, physical_points

def mask_surface(mask, case):
    """Native segmentation geometry; never a replacement for a VascularMD surface."""
    import pyvista as pv
    from skimage.measure import marching_cubes
    slices, origin = crop(mask)
    local = np.pad(mask[slices],1)
    vertices, faces, _, _ = marching_cubes(local.astype(np.uint8),level=.5,allow_degenerate=False)
    indices = vertices + origin - 1
    world = physical_points(indices,case.affine_mm)
    mesh = pv.PolyData(world,np.column_stack((np.full(len(faces),3),faces)).ravel())
    mesh.field_data["geometry_source"] = np.array(["TopBrain native mask marching cubes; unsmoothed; NIfTI affine"])
    # EDT is a voxel estimate. Surface colouring uses a nearby interior voxel only,
    # and never changes vertices or the strict mask.
    distances = ndi.distance_transform_edt(local,sampling=case.spacing_mm)
    nearest = np.rint(vertices).astype(int)
    samples = ndi.maximum_filter(distances,size=3)[tuple(nearest.T)]
    mesh.point_data["display_diameter_estimate_mm"] = 2*samples
    return mesh
