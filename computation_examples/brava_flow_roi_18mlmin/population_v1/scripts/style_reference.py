"""Camera layout copied unchanged from ../microbubble/scripts/render_results.py."""
import numpy as np
ANIMATION_ZOOM=1.30

def projected_bounds(plot,points):
    matrix=plot.camera.GetCompositeProjectionTransformMatrix(plot.renderer.GetTiledAspectRatio(),-1,1)
    matrix=np.array([[matrix.GetElement(i,j) for j in range(4)] for i in range(4)])
    clip=np.column_stack((points,np.ones(len(points))))@matrix.T
    ndc=clip[:,:2]/clip[:,3,None]
    width,height=plot.window_size
    pixels=np.column_stack(((ndc[:,0]+1)*width/2,(1-ndc[:,1])*height/2))
    return np.array([pixels.min(axis=0),pixels.max(axis=0)])

def enlarge_animation(plot,surface):
    """Exact 1.3x orthographic display zoom; retain direction and physical data."""
    old_scale=plot.camera.parallel_scale
    before=projected_bounds(plot,surface.points)
    plot.camera.parallel_scale=old_scale/ANIMATION_ZOOM
    enlarged=projected_bounds(plot,surface.points)
    # Center within the existing header/footer whitespace, without changing view direction.
    target=np.array([960.,558.])
    delta=target-enlarged.mean(axis=0)
    forward=np.asarray(plot.camera.direction);forward/=np.linalg.norm(forward)
    right=np.cross(forward,np.asarray(plot.camera.up));right/=np.linalg.norm(right)
    up=np.cross(right,forward)
    shift=(-delta[0]*right+delta[1]*up)*(2*plot.camera.parallel_scale/plot.window_size[1])
    plot.camera.position=np.asarray(plot.camera.position)+shift
    plot.camera.focal_point=np.asarray(plot.camera.focal_point)+shift
    after=projected_bounds(plot,surface.points)
    ratios=(after[1]-after[0])/(before[1]-before[0])
    assert np.allclose(ratios,ANIMATION_ZOOM,rtol=0,atol=1e-10),ratios
    assert np.all(after[0]>=[40,108]) and np.all(after[1]<=[1880,1008]),after
    return dict(original_bounds_px=before.tolist(),enlarged_bounds_px=after.tolist(),
                projected_size_ratios=ratios.tolist(),original_parallel_scale=old_scale,
                parallel_scale=plot.camera.parallel_scale,display_translation_px=delta.tolist(),
                camera_position=list(plot.camera.position),camera_focal_point=list(plot.camera.focal_point))
