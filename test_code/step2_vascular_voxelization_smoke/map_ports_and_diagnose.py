#!/usr/bin/python3
"""Analyze actual Palabos flags, create a one-link cap opening, export VTK."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree
from prepare_port_contract import read_stl, projected_inside, write_json, PORT_IDS

SIX = ndimage.generate_binary_structure(3, 1)
DIRECTIONS = np.concatenate([np.eye(3,dtype=int), -np.eye(3,dtype=int)])

def segment_hits_triangles(start, end, triangles):
    """Moller-Trumbore finite segment intersection; returns original patch-local ID."""
    hit = np.full(len(start), -1, dtype=int)
    direction = end-start
    for i, (a,b,c) in enumerate(triangles):
        e1=b-a; e2=c-a; h=np.cross(direction,e2); determinant=h@e1
        valid=abs(determinant)>1e-12
        inv=np.zeros(len(start)); inv[valid]=1/determinant[valid]
        s=start-a; u=np.einsum('ij,ij->i',s,h)*inv
        q=np.cross(s,e1); v=np.einsum('ij,ij->i',direction,q)*inv
        t=(q@e2)*inv
        yes=valid & (u>=-1e-10) & (v>=-1e-10) & (u+v<=1+1e-10) & (t>=0) & (t<=1)
        hit[(hit<0)&yes]=i
    return hit

def image_output(path, data, name, origin, dx):
    from vtkmodules.vtkCommonDataModel import vtkImageData
    from vtkmodules.vtkIOXML import vtkXMLImageDataWriter
    from vtkmodules.util.numpy_support import numpy_to_vtk
    image=vtkImageData(); nz,ny,nx=data.shape
    image.SetDimensions(nx+1,ny+1,nz+1)
    # Integer Palabos node i is the CENTER of cell i in diagnostic image.
    image.SetOrigin(*(origin-0.5*dx)); image.SetSpacing(dx,dx,dx)
    array=numpy_to_vtk(data.ravel(),deep=False); array.SetName(name)
    image.GetCellData().SetScalars(array)
    writer=vtkXMLImageDataWriter(); writer.SetFileName(str(path)); writer.SetInputData(image)
    writer.SetDataModeToBinary(); writer.SetCompressorTypeToZLib()
    assert writer.Write()==1

def cap_output(path, triangles, ids, labels):
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkCellArray
    from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
    from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    points=vtkPoints();points.SetData(numpy_to_vtk(triangles.reshape(-1,3),deep=True))
    connectivity=np.arange(3*len(triangles),dtype=np.int64)
    offsets=np.arange(0,len(connectivity)+1,3,dtype=np.int64)
    cells=vtkCellArray();cells.SetData(numpy_to_vtkIdTypeArray(offsets,deep=True),numpy_to_vtkIdTypeArray(connectivity,deep=True))
    data=vtkPolyData();data.SetPoints(points);data.SetPolys(cells)
    for name,values in [('STLTriangleID',ids),('PortLabel',labels)]:
        a=numpy_to_vtk(np.asarray(values),deep=True);a.SetName(name);data.GetCellData().AddArray(a)
    w=vtkXMLPolyDataWriter();w.SetFileName(str(path));w.SetInputData(data);w.SetDataModeToBinary();assert w.Write()==1

def center_output(path, ports):
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkCellArray
    from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
    from vtkmodules.util.numpy_support import numpy_to_vtk
    pts=vtkPoints();pts.SetData(numpy_to_vtk(np.array([ports[k]['center'] for k in PORT_IDS]),deep=True))
    cells=vtkCellArray()
    for i in range(4):cells.InsertNextCell(1);cells.InsertCellPoint(i)
    data=vtkPolyData();data.SetPoints(pts);data.SetVerts(cells)
    for name,values in [('PortLabel',np.arange(1,5,dtype=np.int32)),('OutwardNormal',np.array([ports[k]['normal'] for k in PORT_IDS]))]:
        a=numpy_to_vtk(values,deep=True);a.SetName(name);data.GetPointData().AddArray(a)
        if name=='OutwardNormal':data.GetPointData().SetVectors(a)
    w=vtkXMLPolyDataWriter();w.SetFileName(str(path));w.SetInputData(data);w.SetDataModeToBinary();assert w.Write()==1

def boundary_counts(data):
    return {k:int(np.count_nonzero(v)) for k,v in zip(
        ['x_min','x_max','y_min','y_max','z_min','z_max'],
        [data[:,:,0],data[:,:,-1],data[:,0,:],data[:,-1,:],data[0,:,:],data[-1,:,:]])}

def diagnose(run):
    settings=json.loads((run/'run_settings.json').read_text())
    geom=json.loads((run/'diagnostics/palabos_geometry.json').read_text())
    caps=json.loads((run/'diagnostics/cap_identification.json').read_text()); ports=caps['ports']
    tri=read_stl(settings['input_stl']); dx=geom['effective_dx_m']; origin=np.array(geom['physical_origin_m'])
    bbox=np.array(settings['physical_bbox_m']); expected_origin=bbox[0]-(geom['margin']+geom['extra_layer'])*dx
    vertices=np.fromfile(run/'diagnostics/palabos_vertices_preinflate.f64',dtype='<f8').reshape(-1,3)
    original=np.unique(tri.reshape(-1,3),axis=0); expected_lu=(original-origin)/dx
    residual_forward=float(cKDTree(vertices).query(expected_lu)[0].max())
    residual_reverse=float(cKDTree(expected_lu).query(vertices)[0].max())
    transform_pass=bool(np.max(abs(origin-expected_origin))<1e-18 and
        abs(dx-settings['effective_dx_m_expected'])<1e-20 and
        max(residual_forward,residual_reverse)<1e-8 and len(vertices)==len(original))
    transform=dict(status='PASS' if transform_pass else 'UNVERIFIED',
        formula='p_LU = (p_m - physical_origin_m) / effective_dx_m',
        scale_per_m=1/dx,translation_lu=(-origin/dx).tolist(),physical_origin_m=origin.tolist(),
        bbox_min_m=bbox[0].tolist(),margin=geom['margin'],extra_layer=geom['extra_layer'],
        original_unique_vertices=len(original),palabos_vertices=len(vertices),
        max_bidirectional_vertex_residual_lu=max(residual_forward,residual_reverse),
        normals='unchanged by positive isotropic affine scale',area_scale=1/dx**2,
        vti_cell_data_origin_m=(origin-0.5*dx).tolist(),vti_cell_center_formula='origin_m + [i,j,k] * dx',
        inflation_lu=geom['inflate_lu'],inflation_m=geom['inflate_lu']*dx)
    write_json(run/'diagnostics/transform_check.json',transform)
    assert transform_pass, 'PHYSICAL_TO_LATTICE_TRANSFORM=UNVERIFIED; no opening'
    shape=tuple(reversed(geom['lattice_shape']))
    closed=np.fromfile(run/'diagnostics/closed_flag_matrix.u8',dtype=np.uint8).reshape(shape)
    native=np.memmap(run/'diagnostics/palabos_native_flags.u8',dtype=np.uint8,mode='r',shape=shape)
    assert np.count_nonzero(closed)==geom['closed_fluid_voxels']
    assert np.array_equal(closed,((native==3)|(native==4)).astype(np.uint8))
    image_output(run/'closed_flag_matrix.vti',closed,'ClosedFluid',origin,dx)
    cc,ncc=ndimage.label(closed,structure=SIX)
    counts=np.bincount(cc.ravel()); component_sizes=sorted(counts[1:].tolist(),reverse=True);del cc
    bounds=boundary_counts(closed)
    # Counterfactual ONLY: exact global, one-rank processing order of official CopyFromNeighbor.
    xcounter=[]
    for dest,source in [(0,1),(1,2),(shape[2]-2,shape[2]-3),(shape[2]-1,shape[2]-3)]:
        difference=closed[:,:,dest]!=closed[:,:,source]
        xcounter.append(dict(destination_x=dest,source_x=source,changed_voxels=int(difference.sum()),
            would_add_fluid_voxels=int(((closed[:,:,dest]==0)&(closed[:,:,source]==1)).sum())))
    topology=dict(connectivity='6 face neighbors (strict; no diagonal-only connections)',
        closed_connected_components=int(ncc),closed_component_sizes=component_sizes,
        closed_domain_boundary_fluid_counts=bounds,official_x_copy_counterfactual=xcounter,
        closed_x_end_copy_executed=False,fluid_timesteps_run=0)
    write_json(run/'diagnostics/topology_check.json',topology)
    if ncc != 1:
        write_json(run/'diagnostics/voxelization_summary.json',dict(settings,**geom,
            closed_connected_components=int(ncc),step2_auto_check='FAIL',step2_status='PARTIAL',
            reason='Closed lumen disconnected; four-port opening not attempted'))
        raise RuntimeError('CLOSED_LUMEN_CONNECTIVITY=FAIL; no opening')
    assert not any(bounds.values()), 'Closed fluid reaches bounding box; inspect before opening'
    assert caps['cap_triangle_sets_disjoint'] and all(p['identification_pass'] for p in ports.values())
    opened=closed.copy(); labels=np.zeros(shape,dtype=np.uint8); overlap=0; all_cap_ids=[]; all_cap_labels=[]
    for name in PORT_IDS:
        p=ports[name]; label=p['label']; ids=np.array(p['triangle_ids']); cap=(tri[ids]-origin)/dx
        c=(np.array(p['center'])-origin)/dx; n=np.array(p['normal'])
        lower=np.maximum(np.floor(cap.min(axis=(0,1))).astype(int)-2,0)
        upper=np.minimum(np.ceil(cap.max(axis=(0,1))).astype(int)+2,np.array(geom['lattice_shape'])-1)
        grid=np.stack(np.meshgrid(*(np.arange(a,b+1) for a,b in zip(lower,upper)),indexing='ij'),axis=-1).reshape(-1,3)
        s=(grid-c)@n
        candidate=(closed[tuple(grid[:,::-1].T)]==0)&(s>0)&(s<=max(abs(n))+1e-10)
        grid=grid[candidate];grid=grid[projected_inside(grid,cap,n)]
        witnesses={}
        for direction in DIRECTIONS:
            inside=grid+direction
            valid=((inside>=0)&(inside<np.array(geom['lattice_shape']))).all(axis=1)
            indices=np.flatnonzero(valid)
            indices=indices[closed[tuple(inside[indices,::-1].T)]==1]
            indices=indices[((inside[indices]-c)@n)<=0]
            hits=segment_hits_triangles(inside[indices],grid[indices],cap)
            for index,hit in zip(indices[hits>=0],hits[hits>=0]):
                key=tuple(grid[index]); witnesses.setdefault(key,(inside[index],int(ids[hit])))
        changed=np.array(list(witnesses),dtype=int).reshape(-1,3)
        assert len(changed),f'{name}: empty one-link opening; no enlargement attempted'
        lookup=tuple(changed[:,::-1].T)
        overlap+=int(np.count_nonzero(labels[lookup])); labels[lookup]=label;opened[lookup]=1
        nearest_i=int(np.argmin(np.linalg.norm(changed-c,axis=1)));nearest=changed[nearest_i]
        distances=np.linalg.norm(changed-c,axis=1);signed=(changed-c)@n
        # Direct face-adjacency witnesses are stronger than a component-count-only check.
        link_records=[]
        for node,(neighbor,tid) in witnesses.items():
            link_records.append(dict(opened_voxel=list(map(int,node)),closed_lumen_neighbor=neighbor.tolist(),stl_triangle_id=tid))
        write_json(run/'diagnostics'/f'{name}_opening_witnesses.json',link_records)
        p.update(lattice_center=c.tolist(),lattice_normal=n.tolist(),opened_port_voxels=len(changed),
            nearest_opened_voxel_lu=nearest.tolist(),nearest_opened_voxel_m=(origin+nearest*dx).tolist(),
            distance_in_voxels=float(distances[nearest_i]),
            outside_signed_distance_range_lu=[float(signed.min()),float(signed.max())],
            all_projected_centers_inside_actual_cap=bool(projected_inside(changed,cap,n).all()),
            direct_closed_lumen_neighbor_count=len(link_records),exact_cap_crossing_link_count=len(link_records),
            added_layers='one exterior lattice link; no reservoir',
            mapping_pass=bool(distances[nearest_i]<=settings['center_mapping_max_distance_lu']))
        all_cap_ids.extend(ids.tolist());all_cap_labels.extend([label]*len(ids))
        print(name,'new_voxels',len(changed),'nearest_distance_LU',p['distance_in_voxels'],flush=True)
    changed=(opened!=closed)
    unlabeled=int(np.count_nonzero(changed&(labels==0)))
    labeled_unchanged=int(np.count_nonzero((labels!=0)&~changed))
    negative=int(np.count_nonzero(opened<closed))
    cc,opened_ncc=ndimage.label(opened,structure=SIX)
    patch_components={name:np.unique(cc[labels==ports[name]['label']]).tolist() for name in PORT_IDS};del cc
    opening_locality=bool(unlabeled==0 and labeled_unchanged==0 and negative==0 and overlap==0 and
        all(p['all_projected_centers_inside_actual_cap'] and p['mapping_pass'] and p['exact_cap_crossing_link_count']==p['opened_port_voxels'] for p in ports.values()))
    topology.update(opened_connected_components=int(opened_ncc),ports_connected_component_ids=patch_components,
        all_ports_connected_to_same_closed_lumen=bool(opened_ncc==1 and all(v==[1] for v in patch_components.values())),
        opened_domain_boundary_fluid_counts=boundary_counts(opened),
        added_fluid_voxels=int(changed.sum()),removed_fluid_voxels=negative,
        unlabeled_changed_voxels=unlabeled,labeled_unchanged_voxels=labeled_unchanged,port_overlap_voxels=overlap,
        port_opening_locality='PASS' if opening_locality else 'FAIL',
        wall_changes_outside_four_recovered_caps=unlabeled,
        unrelated_wall_breach_check='All changed nodes have exact cap-intersection and direct closed-lumen witnesses; no other nodes change')
    topology['opened_domain_boundary_port_label_counts']={k:np.bincount(a.ravel(),minlength=5).tolist() for k,a in [
        ('x_min',labels[:,:,0]),('x_max',labels[:,:,-1]),('y_min',labels[:,0,:]),
        ('y_max',labels[:,-1,:]),('z_min',labels[0,:,:]),('z_max',labels[-1,:,:])]}
    write_json(run/'diagnostics/topology_check.json',topology)
    image_output(run/'opened_flag_matrix.vti',opened,'OpenedFluid',origin,dx)
    image_output(run/'port_label_field.vti',labels,'PortLabel',origin,dx)
    cap_output(run/'cap_triangles.vtp',tri[all_cap_ids],all_cap_ids,all_cap_labels)
    center_output(run/'mapped_centers_normals.vtp',ports)
    opened.tofile(run/'diagnostics/opened_flag_matrix.u8');labels.tofile(run/'diagnostics/port_label_field.u8')
    auto=bool(opening_locality and opened_ncc==1 and ncc==1 and transform_pass and
              caps['cap_triangle_sets_disjoint'] and not any(bounds.values()))
    summary=dict(settings,**geom,closed_connected_components=int(ncc),opened_fluid_voxels=int(opened.sum()),
        ports=ports,unintended_x_end_opening=False,unlabeled_changed_voxels=unlabeled,port_overlap_voxels=overlap,
        port_patch_disjointness='PASS' if overlap==0 and caps['cap_triangle_sets_disjoint'] else 'FAIL',
        port_opening_locality=topology['port_opening_locality'],
        step2_auto_check='PASS' if auto else 'FAIL',step2_status='AUTO_PASS_HUMAN_PENDING' if auto else 'PARTIAL',
        human_paraview_review='PENDING',physiological_port_roles='ASSUMED',
        vtk_coordinate_unit='m',vtk_data_association='CellData; cubes centered at integer lattice nodes')
    write_json(run/'diagnostics/voxelization_summary.json',summary)
    assert auto,'STEP2_AUTO_CHECK=FAIL'
    print('GEOMETRY_MAPPING_AUTO_CHECK=PASS; independent audit and human ParaView review remain',flush=True)

def export_existing(run):
    """Export-only path; preserves original Palabos flags and all opening witnesses."""
    s=json.loads((run/'diagnostics/voxelization_summary.json').read_text())
    shape=tuple(reversed(s['lattice_shape']));origin=np.array(s['physical_origin_m']);dx=s['effective_dx_m']
    for file,field in [('closed_flag_matrix','ClosedFluid'),('opened_flag_matrix','OpenedFluid'),('port_label_field','PortLabel')]:
        data=np.memmap(run/'diagnostics'/(file+'.u8'),dtype=np.uint8,mode='r',shape=shape)
        image_output(run/(file+'.vti'),data,field,origin,dx)
    ids=[];labels=[]
    for name in PORT_IDS:
        p=s['ports'][name];ids.extend(p['triangle_ids']);labels.extend([p['label']]*len(p['triangle_ids']))
    tri=read_stl(s['input_stl']);cap_output(run/'cap_triangles.vtp',tri[ids],ids,labels)
    center_output(run/'mapped_centers_normals.vtp',s['ports'])
    print('VTK_EXPORT_ONLY_COMPLETE: compressed inline binary; raw flags unchanged',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);ap.add_argument('--export-only',action='store_true')
    a=ap.parse_args();export_existing(a.run) if a.export_only else diagnose(a.run)
