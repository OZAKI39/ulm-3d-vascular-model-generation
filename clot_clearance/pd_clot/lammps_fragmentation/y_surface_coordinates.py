"""Geometry-only surface contract and bounded assembly checks; never apply loads."""
import itertools
import json
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
from .config import dump
from .y_bifurcation_sources import OUT, RESULTS
from .y_mural_clot_geometry import distances, profile


def validate_and_map(sources,g,c,solid,top,attach,*,out=OUT):
    OUT=out
    n=g['clot_count'];X=g['X'][:n];A=g['X'][n:];dx=c['dx_m']
    assert np.allclose(g['X']/dx,np.rint(g['X']/dx),atol=1e-10,rtol=0)
    assert len(np.unique(g['X'],axis=0))==len(g['X']) and np.isfinite(g['X']).all()
    ids=sources.mesh.find_containing_cell(X);aid=sources.mesh.find_containing_cell(A)
    assert (ids>=0).all() and (aid<0).all()
    # Verify 27 support points per retained full voxel, no resolution campaign.
    offsets=np.array(list(itertools.product([-.5,0,.5],repeat=3)))*dx
    corner_inside=sources.mesh.find_containing_cell((X[:,None,:]+offsets).reshape(-1,3))>=0
    assert corner_inside.all()
    sources.sample((X[:,None,:]+offsets).reshape(-1,3))  # strict barycentric membership, no VTK-only acceptance
    wa=distances(sources.exterior,attach.points)
    cent=attach.cell_centers().points
    wac=distances(sources.exterior,cent)
    topdist=distances(sources.exterior,top.points)
    wall_nodes_distance=-distances(sources.exterior,X)
    anchor_distance=distances(sources.exterior,A)
    assert np.max(np.abs(wa))<2e-9 and np.max(np.abs(wac))<2e-9
    assert topdist.max()<2e-9 and anchor_distance.min()>0
    # Closed reference surface must have one body, with no port overlap.
    con=solid.connectivity();assert len(np.unique(con['RegionId']))==1
    footcon=attach.connectivity();assert len(np.unique(footcon['RegionId']))==1
    gaps={name:float(cKDTree(face.points).query(solid.points)[0].min()) for name,face in sources.faces.items() if name!='WALL'}
    assert min(gaps.values())>.005 and solid.bounds[1]<-.005
    types=g['bond_types'];b=g['bonds'];length=np.linalg.norm(g['X'][b[:,1]]-g['X'][b[:,0]],axis=1)
    assert np.all(b[types==1]<n) and np.all(b[types==3]>=n)
    assert np.all(b[types==2,0]<n) and np.all(b[types==2,1]>=n)
    assert np.unique(b,axis=0).shape==b.shape and length.min()>0 and length.max()<=c['delta_m']*(1+1e-10)
    degree=np.bincount(b.ravel(),minlength=len(g['X']));assert degree.min()>0
    report=dict(status='PASS',actual_footpoint_m=c['actual_wall_footpoint_m'],nominal_footpoint_correction_m=c['nominal_footpoint_correction_m'],
                attachment_vertex_gap_max_m=float(np.abs(wa).max()),attachment_vertex_gap_median_m=float(np.median(np.abs(wa))),
                attachment_face_centroid_gap_max_m=float(np.abs(wac).max()),maximum_surface_penetration_m=float(max(topdist.max(),0)),
                particle_center_wall_distance_median_m=float(np.median(wall_nodes_distance)),particle_center_wall_distance_min_m=float(wall_nodes_distance.min()),
                anchor_center_wall_distance_range_m=[float(anchor_distance.min()),float(anchor_distance.max())],
                material_centers_inside=int((ids>=0).sum()),anchor_centers_outside=int((aid<0).sum()),
                voxel_support_samples_checked=len(corner_inside),voxel_support_samples_inside=int(corner_inside.sum()),
                attachment_patch_components=1,attachment_patch_area_m2=float(attach.area),
                min_bond_degree=int(degree.min()),max_bond_degree=int(degree.max()),
                bond_length_range_m=[float(length.min()),float(length.max())],invalid_type_pairs=0,duplicate_nodes=0,duplicate_bonds=0,
                bond_segment_policy='BULK: exact reference solid boundary intersection exclusion. INTERFACE: crosses actual attachment patch and no exposed face. WALL: seven samples outside lumen and within footprint.',
                voxel_note='Node volumes are equal PD quadrature weights. Wall-adjacent cubes are conservatively retained only if their corners are in native lumen; smooth fitted surface and particle stair-step differ.',
                surface_tolerance_m=2e-9,CLOT_PRESENT_CFD_SOLVED='NOT_RUN')
    dump(OUT/'GEOMETRY_AND_WALL_FIT_CHECK.json',report)
    # Slice area measured from actual saved solid by triangulated planar contour.
    def section_area(mesh):
        cut=mesh.slice(normal=(1,0,0),origin=(c['x0_m'],0,0)).clean()
        import vtk
        f=vtk.vtkContourTriangulator();f.SetInputData(cut);f.Update()
        return float(pv.wrap(f.GetOutput()).area)
    clot_area=section_area(solid);lumen_area=section_area(sources.exterior)
    col=int(np.rint(c['x0_m']/dx));discrete_area=int((np.rint(X[:,0]/dx).astype(int)==col).sum())*dx**2
    clearance=dict(status='PASS',original_domain='UNOBSTRUCTED_Y_LUMEN',assembled_domain='Y_LUMEN_WITH_STATIC_MURAL_SOLID_OVERLAY',
                   CLOT_PRESENT_CFD_SOLVED='NOT_RUN',minimum_distance_to_port_surface_vertices_m=gaps,
                   center_section_x_m=c['x0_m'],center_section_clot_area_m2=clot_area,center_section_lumen_area_m2=lumen_area,
                   center_section_area_occupancy=clot_area/lumen_area,discrete_center_column_area_m2=discrete_area,
                   discrete_center_column_area_occupancy=discrete_area/lumen_area,
                   occupancy_definition='Intersection of actual fitted reference solid with x=x0 plane divided by native lumen cross-section area at same plane; separate lattice-column dx^2 estimate',
                   opposite_wall_clearance_lower_bound_m=float(solid.bounds[4]+c['Rv_m']-.00003),
                   boundary_sample_max_penetration_m=report['maximum_surface_penetration_m'],nonocclusive=True)
    assert 0<clot_area/lumen_area<.3
    dump(OUT/'LUMEN_INTERSECTION_AND_CLEARANCE.json',clearance)
    surface_ids=np.flatnonzero(g['surface']);sx=g['X'][surface_ids]
    wall=sources.faces['WALL'].compute_normals(cell_normals=True,point_normals=False,consistent_normals=True,split_vertices=False)
    nearest=wall.find_closest_cell(sx);outward=np.asarray(wall.cell_normals[nearest],float)
    radial=sx.copy();radial[:,0]=0;radial/=np.linalg.norm(radial,axis=1)[:,None]
    outward[np.einsum('ij,ij->i',outward,radial)<0]*=-1
    outward/=np.linalg.norm(outward,axis=1)[:,None]
    axial=np.tile([1.,0.,0.],(len(sx),1));axial-=np.einsum('ij,ij->i',axial,outward)[:,None]*outward
    axial/=np.linalg.norm(axial,axis=1)[:,None];circ=np.cross(outward,axial)
    u,p=sources.sample(sx);direction=u/np.maximum(np.linalg.norm(u,axis=1)[:,None],1e-30)
    assert np.allclose(np.linalg.norm(g['normals'][surface_ids],axis=1),1) and np.all(g['areas'][surface_ids]>0)
    assert np.isclose(g['areas'].sum(),top.area,rtol=1e-12)
    assert np.allclose(np.einsum('ij,ij->i',axial,outward),0,atol=1e-12)
    np.savez_compressed(OUT/'SURFACE_NODES_AND_NORMALS.npz',surface_node_id=surface_ids+1,position_xyz_m=sx,
                        normal_xyz=g['normals'][surface_ids],A_i_m2=g['areas'][surface_ids],wall_distance_m=-distances(sources.exterior,sx),
                        axial_tangent_xyz=axial,circumferential_tangent_xyz=circ,wall_to_lumen_normal_xyz=-outward,
                        local_flow_direction=direction,background_velocity_m_s=u,background_pressure_Pa=p,
                        geometry_id=c['geometry_id'],source_geometry_id=sources.geometry_id,source_mesh_id=sources.mesh_id,
                        source_flow_id=sources.flow_id,source_frame=c['coordinate_frame'],source_units='SI',
                        semantics='INITIAL_EXPOSED_SURFACE_ONLY; NO APPLIED LOAD; NO TRACTION CERTIFICATION')
    schema=dict(schema_version=1,status='READY',applied_load=False,CLOT_PRESENT_CFD_SOLVED='NOT_RUN',
                provider_status={'BACKGROUND_Y_FLOW':'REQUIRES_CLOT_PRESENT_CFD_FOR_TRACTION','ANALYTIC_STREAMING':'SCHEMA_ONLY_NOT_APPLIED',
                                 'BUBBLE_STREAMING':'NOT_YET','CFD_STREAMING':'NOT_YET'},
                required=['load_source','surface_node_id','position_xyz_m','time_s','source_geometry_id','source_flow_id','coordinate_frame','source_units','validity_domain'],
                load_alternatives=[['force_xyz_N'],['traction_xyz_Pa','nodal_area_m2']],conversion='force_xyz_N = traction_xyz_Pa * nodal_area_m2; use initial surface node IDs in existing LAMMPS addforce provider',
                coordinate_frame=c['coordinate_frame'],transform_to_source=np.eye(4),geometry_id=c['geometry_id'],source_geometry_id=sources.geometry_id,source_flow_id=sources.flow_id,
                surface_area_method=g['info']['surface_area_rule'],exposed_area_sum_m2=float(g['areas'].sum()),attachment_area_excluded_m2=float(attach.area),
                local_tangents='Project global axial e_x to nearest actual wall triangle tangent; circumferential = outward_wall_normal cross axial; wall-to-lumen = negative wall normal',
                initial_exposed_surface_supported=True,attachment_surface_loadable=False,new_fracture_surface_update='NOT_SUPPORTED',
                node_tag_overlap='A volume node can be near both boundaries at a thin rim. Exposed and attachment triangle areas are disjoint; only exposed triangle areas enter A_i.',
                validity_domain='One frozen initial mural geometry; no motion or fracture surface tracking. Provider must reject mismatched identities, non-SI units, nonfinite or outside-domain values.',
                CURVED_INTERFACE_STRENGTH_VALIDATION='NOT_ESTABLISHED')
    dump(OUT/'SURFACE_LOAD_INTERFACE_SCHEMA.json',schema)
    return report,clearance
