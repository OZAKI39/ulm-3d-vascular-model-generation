"""Independent barycentric-gradient, surface mapping, source and video checks."""
from pathlib import Path
import hashlib,json
import numpy as np
import pyvista as pv
import imageio_ffmpeg
from PIL import Image
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    c=json.loads((HERE/'COMPUTE_VALIDATION.json').read_text())
    r=json.loads((HERE/'MEDIA_VALIDATION.json').read_text())
    assert c['all_pass'] and r['all_pass']
    for name,digest in c['source_lock'].items():assert sha(ROOT/name)==digest,name
    for name,digest in c['files'].items():assert sha(HERE/name)==digest,name
    assert r['compute_record_sha256']==sha(HERE/'COMPUTE_VALIDATION.json')
    assert r['script_sha256']==sha(HERE/'render_shear_rate.py')
    assert c['script_sha256']==sha(HERE/'prepare_shear_rate.py')
    a=np.load(ROOT/'input_data/frozen_flow/flow_arrays_si.npz')
    x=a['points_m'];t=a['tetra'];u=a['velocity_m_s']
    # Independent shape-function derivatives from cross products, not the
    # production core's batched linear solves.
    edges=x[t[:,1:]]-x[t[:,:1]];aa,bb,cc=edges[:,0],edges[:,1],edges[:,2]
    det=np.einsum('ni,ni->n',aa,np.cross(bb,cc));assert np.all(det>0)
    shape_grad=np.stack([np.cross(bb,cc),np.cross(cc,aa),np.cross(aa,bb)],axis=1)/det[:,None,None]
    g=np.einsum('nai,naj->nij',u[t[:,1:]]-u[t[:,:1]],shape_grad)
    expected=np.sqrt(2*(g[:,0,0]**2+g[:,1,1]**2+g[:,2,2]**2)
        +(g[:,0,1]+g[:,1,0])**2+(g[:,0,2]+g[:,2,0])**2+(g[:,1,2]+g[:,2,1])**2)
    grid=pv.read(HERE/'data/strain_shear_rate_volume_si.vtu')
    actual=grid.cell_data['ShearRate_s_inv']
    assert grid.n_cells==len(t)==371402 and np.array_equal(grid.points,x)
    assert np.array_equal(grid['Velocity'],u)
    assert np.allclose(actual,expected,rtol=1e-12,atol=1e-9)
    d=.5*(g+g.swapaxes(1,2));stored=grid.cell_data['StrainRateTensor_s_inv'].reshape(-1,3,3)
    assert np.allclose(d,stored,rtol=1e-11,atol=1e-9)
    assert np.allclose(actual,np.sqrt(2)*grid.cell_data['StrainRateFrobenius_s_inv'],rtol=1e-14,atol=1e-10)
    surface=pv.read(HERE/'data/shear_rate_exterior_si.vtp')
    tri=surface.point_data['GlobalNodeID_zero_based'][surface.faces.reshape(-1,4)[:,1:]]
    assert np.array_equal(tri,a['boundary_triangles'])
    owners=surface.cell_data['Parent_tetra_zero_based'];parent=t[owners]
    assert np.all(np.any(tri[:,:,None]==parent[:,None,:],axis=2))
    assert np.array_equal(surface.cell_data['ShearRate_s_inv'],actual[owners])
    assert np.array_equal(surface.cell_data['Boundary_tag'],a['facet_tags'])
    assert surface.n_cells==45704
    assert actual.min()>=0 and actual.max()<=r['colorbar_range_s_inv'][1]
    assert r['mapping']['mode']=='UseCellFieldData' and r['mapping']['scalar']=='ShearRate_s_inv'
    assert not r['mapping']['interpolate_scalars_before_mapping']
    assert not r['cell_to_point_averaging'] and not r['scalar_clipping']
    for stem in ['camera','annotations']:
        assert (HERE/f'shear_rate_{stem}.json').read_bytes()==(ROOT/f'results/surface_fields/wss_{stem}.json').read_bytes()
    movie=HERE/r['video']['file'];assert sha(movie)==r['video']['sha256']
    reader=imageio_ffmpeg.read_frames(str(movie),pix_fmt='rgb24');meta=next(reader);count=0;distinct=set()
    for frame in reader:count+=1;distinct.add(hashlib.sha256(frame).hexdigest())
    assert count==len(distinct)==432 and tuple(meta['size'])==(1920,1080)
    assert meta['duration']==18 and meta['fps']==24
    images=[]
    for p in sorted((HERE/'figures').glob('*_4k.png')):
        with Image.open(p) as im:assert im.size==(3840,2160) and min(im.info['dpi'])>299
        images.append(p.name)
    assert len(images)==4
    result=dict(all_pass=True,protected_source_files_verified=len(c['source_lock']),
        independent_gradient_method='Barycentric shape gradients via edge cross products',
        shear_rate_max_abs_difference_s_inv=float(np.max(abs(actual-expected))),
        strain_tensor_max_abs_difference_s_inv=float(np.max(abs(stored-d))),
        all_371402_tetrahedra_checked=True,all_45704_exterior_parents_checked=True,
        no_nodal_averaging_or_scalar_clipping=True,unit='s^-1',
        same_camera_and_labels_as_existing_flow_visualization=True,
        video=dict(**meta,decoded_frames=count,distinct_frames=len(distinct)),figures_4k=images,
        rendering_backend=r['OpenGL'],encoder=r['encoder'])
    (HERE/'LOCAL_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
