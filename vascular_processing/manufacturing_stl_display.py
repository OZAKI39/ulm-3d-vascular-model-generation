"""The saved print STL as a data-only alternative in the current compact UI."""
from dataclasses import replace
from pathlib import Path
import csv
import json

import numpy as np
from scipy.spatial import cKDTree

from .compact_display_adapter import CompactBraVaViewer, compact_record
from .manufacturing_roi import topology_paths
from .print_orientation import transform_points
from .swc_export import read_source
from .topbrain_qc import sha256, write_json
from utils.rodent_vasculature import interactive as ui


def diameter_on_surface(points_mm, smooth_graph):
    """Interpolate fitted SWC radii on the exact closest centerline segment.

    This attaches display scalars only. STL vertices/triangles are never fitted,
    remeshed or replaced. STL itself does not contain a radius array.
    """
    import pyvista as pv
    from vtkmodules.vtkCommonCore import reference
    from vtkmodules.vtkCommonDataModel import vtkStaticCellLocator
    ids=list(smooth_graph); index={n:i for i,n in enumerate(ids)}
    xyzr=np.array([smooth_graph.nodes[n]['coords'] for n in ids])
    edges=np.array([[index[a],index[b]] for a,b in smooth_graph.edges],dtype=np.int64)
    lines=pv.PolyData(xyzr[:,:3],lines=np.c_[np.full(len(edges),2),edges].ravel())
    locator=vtkStaticCellLocator();locator.SetDataSet(lines);locator.BuildLocator()
    closest=[0.,0.,0.];cell=reference(0);sub=reference(0);distance2=reference(0.)
    result=np.empty(len(points_mm)); distances=np.empty(len(points_mm))
    for i,point in enumerate(points_mm):
        locator.FindClosestPoint(point,closest,cell,sub,distance2)
        a,b=edges[int(cell)]; delta=xyzr[b,:3]-xyzr[a,:3]
        weight=np.clip(np.dot(np.asarray(closest)-xyzr[a,:3],delta)/np.dot(delta,delta),0.,1.)
        result[i]=2*((1-weight)*xyzr[a,3]+weight*xyzr[b,3])*1000.
        distances[i]=float(distance2)**.5
    return result,dict(method='Closest native fitted SWC segment, linearly interpolated radius; display scalar only',
        minimum_diameter_um=float(result.min()),maximum_diameter_um=float(result.max()),
        maximum_surface_to_centerline_mm=float(distances.max()))


def load_print_stl(path):
    import pyvista as pv
    path=Path(path).resolve();directory=path.parent
    manifest_path=directory/'manufacturing_manifest.json';transform_path=directory/'print_transform.json'
    if not path.is_file() or not manifest_path.is_file() or not transform_path.is_file():
        raise ValueError('PRINT_STL_PROVENANCE_MISSING: STL requires its manufacturing manifest and print_transform.json')
    m=json.loads(manifest_path.read_text());artifact=json.loads((directory/'artifact_manifest.json').read_text())
    if path!=Path(m['selected_stl']).resolve():
        raise ValueError('PRINT_STL_NOT_THE_MANIFEST_SELECTED_CANDIDATE')
    compensated=Path(m['exports']['compensated']['path'])
    mapping_path=compensated.with_name(compensated.stem+'_mapping.csv')
    smooth_path=Path(m['models']['compensated']['outputs']['swc'])
    native_path=Path(m['orientation_input']['path'])
    for file in [path,manifest_path,transform_path,compensated,mapping_path,smooth_path,native_path]:
        relative=str(file.relative_to(directory))
        if sha256(file)!=artifact['files'].get(relative):
            raise ValueError('PRINT_STL_ARTIFACT_CHANGED: '+relative)
    transform=json.loads(transform_path.read_text());matrix=np.asarray(transform['transform_4x4'],float)
    if matrix.shape!=(4,4) or not np.isfinite(matrix).all():raise ValueError('INVALID_PRINT_TRANSFORM')
    np.testing.assert_allclose(matrix[3],[0,0,0,1],atol=1e-12)
    np.testing.assert_allclose(matrix[:3,:3].T@matrix[:3,:3],np.eye(3),atol=1e-10)
    np.testing.assert_allclose(np.linalg.det(matrix[:3,:3]),1.,atol=1e-10)
    np.testing.assert_allclose(matrix,m['selected_orientation']['transform_4x4'],atol=1e-12)
    printed=pv.read(path)
    if not printed.n_cells or not printed.is_all_triangles:raise ValueError('PRINT_STL_INVALID_TRIANGLES')
    points=np.asarray(printed.points,dtype=float)
    original_points=transform_points(points,np.linalg.inv(matrix))
    # Independent stored native surface proves that inverse-transform is the
    # correct coordinate alignment; no approximate registration is performed.
    native=pv.read(native_path)
    error=max(float(cKDTree(native.points).query(original_points)[0].max()),
              float(cKDTree(original_points).query(native.points)[0].max()))
    if error>2e-4:raise ValueError('PRINT_STL_NATIVE_ALIGNMENT_MISMATCH: '+str(error))
    graph=read_source(compensated).graph
    source=json.loads((directory.parent/'roi_manifest.json').read_text())['source']
    if sha256(Path(source['raw_source']))!=source['raw_sha256']:raise ValueError('BRAVA_SOURCE_CHANGED')
    raw=read_source(Path(source['raw_source'])).graph
    with mapping_path.open() as stream:
        selected=[int(row['original_swc_id']) for row in csv.DictReader(stream)]
    stats={**m['manufacturing_topology'], 'branch_count':len(topology_paths(graph)),
           'bifurcation_count':sum(graph.out_degree(n)>1 for n in graph)}
    item=dict(exports={'compensated':m['exports']['compensated']},selected_original_ids=selected,
        stats=stats,model=m['models']['compensated'],surface_vtk=m['models']['compensated']['outputs']['surface_vtk'],
        display_detail='Saved print STL | original anatomical coordinates',display_status='FULL_CONTEXT_REFERENCE | '+m['status'])
    record=compact_record('PRINT_STL',item,raw,0)
    smooth_graph=read_source(smooth_path).graph
    if set(smooth_graph)!=set(graph) or set(smooth_graph.edges)!=set(graph.edges):
        raise ValueError('FITTED_CENTERLINE_TOPOLOGY_MISMATCH')
    # Replace overlay data too: original pre-fit curves can lie outside the STL.
    # Source IDs still identify the exact original edges highlighted on the left;
    # positions/radii for the right-hand overlay come from the paired native fit.
    fitted=np.array([smooth_graph.nodes[n]['coords'] for n in graph])
    positions=fitted[:,:3]*1000.;radii=fitted[:,3]*1000.
    ports=tuple(replace(port,intersection_position_um=tuple(positions[port.local_node_id]),
                        radius_at_cut_um=float(radii[port.local_node_id])) for port in record.cut_ports)
    record=replace(record,local_node_positions_um=positions,local_node_radius_um=radii,
        local_edge_points_um=positions[record.local_edges],local_edge_radius_um=radii[record.local_edges],
        cut_ports=ports,radius_features={f'r{q}':float(np.percentile(radii,q)) for q in (10,25,50,75,90)})
    surface=printed.copy(deep=True);surface.points=original_points*1000.
    surface.point_data[ui.DIAMETER_SCALAR_NAME],radius_info=diameter_on_surface(original_points,smooth_graph)
    # Smooth shading is the existing style. Normals change no point or triangle.
    faces=surface.faces.copy()
    surface=surface.compute_normals(point_normals=True,cell_normals=False,split_vertices=False,
        consistent_normals=False,auto_orient_normals=False,inplace=False)
    np.testing.assert_array_equal(surface.faces,faces)
    np.testing.assert_array_equal(surface.points,original_points*1000.)
    low=np.minimum(np.asarray(record.bbox_min_um),surface.points.min(axis=0))
    high=np.maximum(np.asarray(record.bbox_max_um),surface.points.max(axis=0))
    record=replace(record,bbox_min_um=tuple(low),bbox_max_um=tuple(high),bbox_center_um=tuple((low+high)/2),bbox_size_um=tuple(high-low))
    data=dict(default_print_candidate='PRINT_STL',candidates={'PRINT_STL':item},context_roi_description='Selected print STL ROI')
    provenance=dict(stl=str(path),stl_sha256=sha256(path),stl_points=printed.n_points,stl_triangles=printed.n_cells,
        stored_print_transform=matrix.tolist(),display_transform=np.linalg.inv(matrix).tolist(),
        inverse_transform_for_display_only=True,native_alignment_max_error_mm=error,
        connectivity_preserved=True,stl_file_modified=False,diameter_colour_source=str(smooth_path),diameter_scalar=radius_info,
        overlay_centerline_source=str(smooth_path),left_highlight_mapping='Original pre-fit source edge IDs; right overlay uses paired fitted coordinates',
        full_context_edges=raw.number_of_edges(),selected_source_edges=record.edge_count)
    return data,raw,[record],surface,provenance


class ManufacturingSTLViewer(CompactBraVaViewer):
    """Swap the existing tube actor's dataset, preserving all rendering properties."""
    output_subdirectory='manufacturing_stl_visualization'

    def __init__(self, path, *, initial_roi=None, show=True, output_dir=None):
        self.stl_path=Path(path).resolve()
        super().__init__(self.stl_path.parent,initial_roi=initial_roi,show=show,view='strict',output_dir=output_dir)

    def load_dataset(self):
        data,raw,records,self.stl_surface,self.stl_provenance=load_print_stl(self.stl_path)
        return data,raw,records

    def select(self,index):
        super().select(index)
        self.plotter.subplot(0,1)
        actors=[actor for actor in self.plotter.renderer.actors.values()
            if getattr(actor,'mapper',None) is not None and hasattr(actor.mapper,'dataset')
            and actor.mapper.scalar_visibility and ui.DIAMETER_SCALAR_NAME in actor.mapper.dataset.point_data]
        if len(actors)!=1:raise ValueError('EXPECTED_EXISTING_DIAMETER_TUBE_ACTOR')
        self.stl_actor=actors[0]
        self.stl_actor.mapper.dataset=self.stl_surface
        values=self.stl_surface.point_data[ui.DIAMETER_SCALAR_NAME]
        self.stl_actor.mapper.scalar_range=(float(values.min()),float(values.max()))
        self.plotter.subplot(0,0);self.plotter.render()

    def run_window(self, *, show=True, smoke_seconds=0):
        report=super().run_window(show=show,smoke_seconds=smoke_seconds)
        report.update(source='manufacturing-stl',stl_provenance=self.stl_provenance,
            roi_surface='Actual saved STL triangles; inverse print transform in memory only',
            differences=['Data-only STL mode; shared current rendering properties retained'],
            displayed_stl_triangles=self.stl_surface.n_cells)
        write_json(self.run/'brava_strict_ui_compatibility.json',report)
        return report
