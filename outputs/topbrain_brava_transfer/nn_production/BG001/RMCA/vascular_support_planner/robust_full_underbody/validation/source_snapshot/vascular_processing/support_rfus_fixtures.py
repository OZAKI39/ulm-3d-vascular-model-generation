"""Permanent RFUS geometry fixtures and an actual-slice horizontal tube test."""
from pathlib import Path
import numpy as np
import trimesh
from .support_rfus_geometry import make_segment
from .support_v3_geometry import union
from .support_inputs import save

NAMES=['full_support_straight_tube','full_support_curved_tube','full_support_S_tube','full_support_Y_branch','full_support_external_port','vertical_overlap_vessels','curtain_lateral_offset','twin_curtain','segmented_curtain','staggered_seam','removal_handle','handle_collision','base_neck_isolation','curtain_removal','trapped_curtain']


def tube(xyz,radius=.6):
    import pyvista as pv
    p=pv.PolyData(np.asarray(xyz,dtype=float));p.lines=np.r_[len(xyz),np.arange(len(xyz))]
    m=p.tube(radius=radius,n_sides=32,capping=True).triangulate()
    result=trimesh.Trimesh(m.points,np.asarray(m.faces).reshape(-1,4)[:,1:],process=True)
    trimesh.repair.fix_normals(result);return result


def branch(xyz,is_port=False,name='SYNTHETIC'):
    xyz=np.asarray(xyz);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(xyz,axis=0),axis=1))];t=np.gradient(xyz,arc,axis=0);t/=np.linalg.norm(t,axis=1)[:,None]
    return dict(branch_id=name,xyz=xyz,arc=arc,radius=np.full(len(xyz),.6),tangent=t,curvature=np.linalg.norm(np.gradient(t,arc,axis=0),axis=1),is_port=is_port)


def fixture(name,profile):
    if name not in NAMES:raise ValueError(name)
    x=np.linspace(90,120,181);y=np.full(len(x),100.)
    if name=='full_support_curved_tube':y+=3*np.sin(np.linspace(0,np.pi,len(x)))
    if name=='full_support_S_tube':y+=2*np.sin(np.linspace(0,2*np.pi,len(x)))
    xyz=np.c_[x,y,np.full(len(x),10.)];b=branch(xyz,name=='full_support_external_port');v=tube(xyz);n=np.gradient(xyz[:,:2],axis=0);n=np.c_[-n[:,1],n[:,0]]/np.linalg.norm(n,axis=1)[:,None]
    top=profile['first_layer_height']+np.floor((9.4-profile['roof_gap_mm']-profile['first_layer_height'])/profile['layer_height'])*profile['layer_height']
    offset=profile['line_width'] if name=='curtain_lateral_offset' else 1.1 if name=='twin_curtain' else 0.
    r=make_segment(xyz[:,:2],np.full(len(x),top),0.,n,profile,2,2,1,offset,name=='twin_curtain');obstacles=[];branches=[b]
    if name in ('segmented_curtain','staggered_seam'):
        cut=trimesh.creation.box([profile['line_width']/2,4,20]);cut.apply_translation([105,100,10])
        if name=='staggered_seam':cut.apply_transform(trimesh.transformations.rotation_matrix(np.radians(30),[0,0,1],[105,100,0]))
        r['mesh']=trimesh.boolean.difference([r['mesh'],cut],engine='manifold')
    if name in ('removal_handle','handle_collision'):
        from .support_v2_primitives import loft
        handle=union([loft([104,100,3],[104,102,8],.66),loft([104,102,7.99],[104,102,14],.66)])
        r['handles']=handle;r['mesh']=union([r['mesh'],handle])
        if name=='handle_collision':obstacles=[tube([[100,102,13],[110,102,13]],1.)]
    if name=='full_support_Y_branch':
        daughter=np.c_[np.linspace(105,120,61),np.linspace(100,110,61),np.full(61,10.)];v=union([v,tube(daughter)])
        nn=np.tile(np.array([-10.,15.])/np.linalg.norm([10,15]),(61,1));other=make_segment(daughter[:,:2],np.full(61,top),0,nn,profile)
        # All three anatomical arcs belong to the verification target.
        branches=[branch(xyz[:91],name='Y_PARENT'),branch(xyz[90:],name='Y_DAUGHTER_1'),branch(daughter,name='Y_DAUGHTER_2')]
        r['daughter']=other;r['mesh']=union([r['mesh'],other['mesh']])
        slit=trimesh.creation.box([profile['line_width']/2,30,20]);slit.apply_translation([105,100,10])
        r['mesh']=trimesh.boolean.difference([r['mesh'],slit],engine='manifold')
    if name=='vertical_overlap_vessels':obstacles=[tube([[104,96,5],[104,104,5]])]
    if name=='trapped_curtain':
        roof=trimesh.creation.box([35,8,1]);roof.apply_translation([105,100,10.5]);obstacles=[roof]
    return dict(support=r,vascular=v,branches=branches,obstacles=obstacles)


def generate(out,profile):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[]
    for name in NAMES:
        f=fixture(name,profile);f['support']['mesh'].export(out/(name+'.stl'));f['vascular'].export(out/(name+'_vascular.stl'))
        if f['obstacles']:union(f['obstacles']).export(out/(name+'_obstacles.stl'))
        rows.append(dict(name=name,closed=bool(f['support']['mesh'].is_volume),volume_mm3=float(f['support']['mesh'].volume)))
    save(out/'manifest.json',dict(fixtures=rows,profile=profile,physical_validation=False));return rows


def actual_straight_test(c):
    from .support_slicer_adapter import discover_engine,engine_profiles
    from .vascular_support_pipeline import sliced
    from .support_qc import parse_actual_paths
    from .support_rfus_travel import events,actual_sections
    from .support_rfus_validation import underbody_coverage
    from .support_rfus_geometry import corridor
    from shapely.geometry import GeometryCollection
    from .support_inputs import table
    out=c['out']/'validation/actual_straight_tube';out.mkdir(parents=True,exist_ok=True);f=fixture('full_support_straight_tube',c['profile'])
    # Extend only the fixture roof ends by one line, retaining the 30 mm target.
    r=f['support'];xy=r['path_xy'].copy();xy[0,0]-=c['profile']['line_width'];xy[-1,0]+=c['profile']['line_width']
    # A 2 mm test foundation and the 2-line neck variant retain the saved
    # profile's first-layer width; never override the machine profile.
    f['support']=make_segment(xy,r['top'],np.full(len(xy),2.),r['lateral'],c['profile'],neck_lines=2)
    foundation=trimesh.creation.box([34,6,2]);foundation.apply_translation([105,100,1])
    # A separate witness column keeps global layers nonempty across the
    # intentional 0.20 mm roof gap, like the real model's untouched box walls.
    # It is outside the tube projection and is never credited in coverage.
    witness=trimesh.creation.box([1,1,12]);witness.apply_translation([92,102,6])
    path=out/'horizontal_tube_and_full_curtain.stl';union([foundation,witness,f['vascular'],f['support']['mesh']]).export(path)
    engine=discover_engine('Bambu',c['discovery'],out);profiles=engine_profiles(c['discovery'],'Bambu',{'enable_support':'0'})
    result=sliced(engine,profiles,path,path,out/'bambu',900)
    if not result.get('slicer_success'):raise RuntimeError('REAL_STRAIGHT_TUBE_FIXTURE_SLICE_FAILED')
    p=events(result['archive']);fp=actual_sections(p,f['support']['mesh'],c['profile'])
    bs=[corridor(b,f['vascular'],np.array([[80,80,0],[130,130,20]]),c['profile'],c['cfg']) for b in f['branches']]
    actual=underbody_coverage(bs,fp,{},c['profile']);table(out/'actual_underbody_coverage.csv',actual['rows'])
    # Sparse control removes middle 20 mm only from the measured footprint,
    # never credits a nominal roof or the vascular extrusion itself.
    from shapely.geometry import box
    sparse={z:p.difference(box(95,90,115,110)) for z,p in fp.items()};negative=underbody_coverage(bs,sparse,{},c['profile'])
    # Also slice a real sparse control, not just a removed-footprint unit check.
    cutter=trimesh.creation.box([20,4,20]);cutter.apply_translation([105,100,10])
    sparse_support=trimesh.boolean.difference([f['support']['mesh'],cutter],engine='manifold')
    sparse_path=out/'horizontal_tube_sparse_control.stl'
    union([foundation,witness,f['vascular'],sparse_support]).export(sparse_path)
    sparse_slice=sliced(engine,profiles,sparse_path,sparse_path,out/'bambu_sparse_control',900)
    if not sparse_slice.get('slicer_success'):raise RuntimeError('REAL_SPARSE_TUBE_CONTROL_SLICE_FAILED')
    sparse_fp=actual_sections(events(sparse_slice['archive']),sparse_support,c['profile'])
    real_negative=underbody_coverage(bs,sparse_fp,{},c['profile'])
    table(out/'sparse_actual_underbody_coverage.csv',real_negative['rows'])
    summary=dict(fixture_layer_witness=dict(center=[92,102,6],extents=[1,1,12],credited_as_support=False),fixture_foundation_mm=[34,6,2],neck_lines=2,profile_overridden=False,archive=str(result['archive']),real_bambu_slice_success=True,full_underbody=actual,sparse_negative_control={k:v for k,v in real_negative.items() if k!='rows'},sparse_archive=str(sparse_slice['archive']),real_sparse_control=True,removed_footprint_control={k:v for k,v in negative.items() if k!='rows'},passed=actual['all_samples_supported'] and real_negative['max_uncovered_underbody_arc_mm']>19.,claim='Actual full and sparse slice coverage fixtures; no physical strength/removal claim')
    save(out/'result.json',summary);return summary
