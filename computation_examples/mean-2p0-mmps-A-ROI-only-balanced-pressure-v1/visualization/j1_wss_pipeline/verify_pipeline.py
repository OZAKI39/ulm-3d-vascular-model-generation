"""Independent local numeric and complete media verification; no CFD execution."""
from pathlib import Path
import hashlib,json,csv
import numpy as np
import pyvista as pv
import imageio_ffmpeg
from PIL import Image
from scipy.spatial.distance import pdist
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    c=json.loads((HERE/'COMPUTE_VALIDATION.json').read_text());m=json.loads((HERE/'MEDIA_VALIDATION.json').read_text())
    assert c['all_pass'] and m['all_pass'] and len(m['videos'])==5
    for name,h in c['source_lock'].items():assert sha(ROOT/name)==h,name
    for name,h in c['files'].items():assert sha(HERE/name)==h,name
    a=np.load(ROOT/'input_data/frozen_flow/flow_arrays_si.npz');wall=pv.read(HERE/'data/J1_wss_pipeline_si.vtp')
    ids=wall.point_data['GlobalNodeID_zero_based'];assert np.array_equal(wall.points,a['points_m'][ids])
    faces=wall.faces.reshape(-1,4)[:,1:];tri=ids[faces];allids=wall.cell_data['Global_boundary_facet_zero_based']
    assert np.all(a['facet_tags'][allids]==1)
    assert np.array_equal(np.sort(tri,axis=1),np.sort(a['boundary_triangles'][allids],axis=1))
    owner=wall.cell_data['Parent_tetra_zero_based'];tet=a['tetra'][owner];x=a['points_m'][tet];u=a['velocity_m_s'][tet]
    # Barycentric edge cross products, independent of production batched solve.
    e=x[:,1:]-x[:,:1];det=np.einsum('ij,ij->i',e[:,0],np.cross(e[:,1],e[:,2]));assert np.all(det>0)
    gradshape=np.stack([np.cross(e[:,1],e[:,2]),np.cross(e[:,2],e[:,0]),np.cross(e[:,0],e[:,1])],axis=1)/det[:,None,None]
    g=np.einsum('nki,nkj->nij',u[:,1:]-u[:,:1],gradshape)
    grad_error=np.max(abs(g.reshape(-1,9)-wall.cell_data['VelocityGradient_s_inv']))
    assert grad_error<1e-8
    xyz=a['points_m'][tri];cross=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);area=np.linalg.norm(cross,axis=1)/2
    n=cross/(2*area[:,None]);direction=xyz.mean(1)-x.mean(1);n[(n*direction).sum(1)<0]*=-1
    assert np.max(abs(n-wall.cell_data['Outward_normal']))<1e-12
    t=c['mu_Pa_s']*np.einsum('nij,nj->ni',g+g.swapaxes(1,2),n)
    projection=np.eye(3)-n[:,:,None]*n[:,None,:];tau=np.einsum('nij,nj->ni',projection,t)
    err=float(np.max(abs(tau-wall.cell_data['WSSVector_Pa'])));assert err<1e-10
    w=np.linalg.norm(tau,axis=1);werr=float(np.max(abs(w-wall.cell_data['WSS_raw_Pa'])));assert werr<1e-10
    assert np.max(abs((tau*n).sum(1)))<1e-11
    node_n=wall.point_data['NodeOutwardNormal'];assert np.max(abs(np.linalg.norm(node_n,axis=1)-1))<1e-12
    full=pv.read(ROOT/'input_data/field_diagnostics/data/wall_wss_si.vtp')
    index=np.searchsorted(full.point_data['GlobalNodeID_zero_based'],ids)
    assert np.array_equal(full.point_data['GlobalNodeID_zero_based'][index],ids)
    assert np.array_equal(wall.point_data['WSS_display_Pa'],full.point_data['WSS_display_Pa'][index])
    center=xyz.mean(1)*1e6;r=np.linalg.norm(center-c['J1_um'],axis=1);assert np.max(r)<c['display_radius_um']
    rows=list(csv.DictReader((HERE/'data/J1_summary.csv').open()))
    for radius in [5.,11.]:
        mask=r<radius;row=next(s for s in rows if s['region']==f'J1_r{radius:g}um' and s['quantity']=='WSS_raw')
        assert int(row['facets'])==mask.sum()
        assert abs(float(row['mean'])-np.average(w[mask],weights=area[mask]))<1e-10
    samples=np.load(HERE/'data/glyph_samples.npz');fi=samples['face_indices'];ni=samples['node_indices']
    min_face=float(pdist(center[fi]).min());min_node=float(pdist(wall.points[ni]*1e6).min())
    assert min_face>=c['glyph_spacing_um']-1e-9 and min_node>=c['glyph_spacing_um']-1e-9
    cameras=json.loads((HERE/'J1_camera.json').read_text());assert len(cameras)==432
    for row in cameras:
        assert row['axis_unit']==[0.,0.,1.] and row['camera_up_unit']==[0.,0.,1.]
        assert row['center_um']==c['J1_um']
        assert max(abs(v) for v in row['projected_bounds'])<=.975001
    labels=json.loads((HERE/'J1_annotations.json').read_text());assert labels['validation']['all_pass']
    for frame in labels['frames']:
        for label in frame['labels']:
            if label['kind']=='port' and label['opacity']>1e-8:assert label['surface_clearance_px']>=22
    videos=[]
    for v in m['videos']:
        p=HERE/v['file'];assert sha(p)==v['sha256'];rr=imageio_ffmpeg.read_frames(str(p));meta=next(rr)
        assert meta['size']==(1920,1080) and meta['fps']==24 and abs(meta['duration']-18)<.05
        count=0;unique=set()
        for raw in rr:count+=1;unique.add(hashlib.sha256(raw).hexdigest())
        assert count==432 and len(unique)==432
        videos.append(dict(file=v['file'],decoded_frames=count,distinct_frames=len(unique),all_pass=True))
    pngs=[]
    for p in sorted((HERE/'figures').glob('*_4k.png')):
        with Image.open(p) as im:assert im.size==(3840,2160)
        pngs.append(p.name)
    result=dict(all_pass=True,independent_gradient_max_abs_error_s_inv=float(grad_error),
        independent_traction_max_abs_error_Pa=err,independent_WSS_max_abs_error_Pa=werr,
        all_source_files_unchanged=True,only_wall_facets=True,real_geometry_and_parent_mapping_verified=True,
        nodal_WSS_identical_to_existing_display=True,unit_node_normals_verified=True,
        source_facet_normal_projection_verified=True,minimum_face_glyph_spacing_um=min_face,
        minimum_node_glyph_spacing_um=min_node,same_camera_framing_verified=True,annotation_clearance_all_frames=True,
        videos=videos,images_4k=pngs,script_sha256=sha(Path(__file__)))
    (HERE/'LOCAL_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
