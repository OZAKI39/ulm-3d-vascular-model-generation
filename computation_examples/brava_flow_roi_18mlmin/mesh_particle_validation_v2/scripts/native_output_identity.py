"""Read native Float32 coordinate metadata via exact, injective input-node identity.

The solver reads double input coordinates but its vtkPoints output defaults to
float. This adapter never interpolates fields or accepts a geometric tolerance.
Native files remain untouched. A rounded-coordinate collision is a hard error.
"""
from pathlib import Path
import hashlib
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
from vascular_validation.postprocess import SolutionMeasurements
from vascular_validation.validation import parse_result_vtu


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identify_native_grid(points, tetra, grid):
    points=np.asarray(points);tetra=np.asarray(tetra);native=np.asarray(grid.points)
    assert points.dtype==np.float64 and np.isfinite(points).all()
    assert grid.n_points==len(points) and grid.n_cells==len(tetra)
    assert np.all(grid.celltypes==10), 'Only the reviewed linear tetrahedral output is supported'
    assert native.dtype in [np.dtype('float32'),np.dtype('float64')]
    encoded=points.astype(native.dtype)
    assert len(np.unique(encoded,axis=0))==len(points), 'Input coordinates collide at native output precision'
    distance,native_to_canonical=cKDTree(encoded.astype(np.float64)).query(native.astype(np.float64))
    assert np.all(distance==0), 'Native coordinates are not the exact declared input encoding'
    assert len(np.unique(native_to_canonical))==len(points), 'Native node correspondence is not bijective'
    canonical_to_native=np.argsort(native_to_canonical)
    assert np.array_equal(native[canonical_to_native],encoded)
    cells=np.asarray(grid.cells).reshape(-1,5);assert np.all(cells[:,0]==4)
    mapped=native_to_canonical[cells[:,1:]]
    assert np.array_equal(np.sort(mapped,axis=1),np.sort(tetra,axis=1)), 'Native tetra rows or physical node membership changed'
    record=dict(status='EXACT_ENCODING_AND_TOPOLOGY_IDENTITY_PASS',native_coordinate_dtype=str(native.dtype),
        canonical_coordinate_dtype=str(points.dtype),point_count=len(points),tetra_count=len(tetra),
        coordinate_encoding='Exact dtype cast of declared input coordinates; no nearest-distance acceptance tolerance',
        point_permutation_identity=bool(np.array_equal(canonical_to_native,np.arange(len(points)))),
        max_coordinate_serialization_error_m=float(np.linalg.norm(native[canonical_to_native].astype(float)-points,axis=1).max()),
        rounded_coordinate_collisions=0,tetra_rows_identical_ignoring_local_vertex_order=True,
        local_order_changed_cells=int(np.any(mapped!=tetra,axis=1).sum()),
        node_map_sha256=hashlib.sha256(canonical_to_native.astype('<i8').tobytes()).hexdigest(),
        gradients_use_canonical_connectivity=True,field_interpolation=False,scientific_values_modified=False)
    return canonical_to_native,record


class CanonicalSolutionMeasurements(SolutionMeasurements):
    def __init__(self,arrays_path,Q_target,Umean):
        super().__init__(arrays_path,Q_target,Umean)
        self.arrays_path=Path(arrays_path)
        self.native_output_identities={}

    def read(self,path):
        grid,velocity,pressure=parse_result_vtu(path)
        mapping,record=identify_native_grid(self.points,self.tetra,grid)
        record.update(native_file_sha256=sha(path),mesh_arrays_sha256=sha(self.arrays_path),reader_source_sha256=sha(__file__))
        self.native_output_identities[str(path)]=record
        return np.asarray(velocity,dtype=float)[mapping],np.asarray(pressure,dtype=float).reshape(-1)[mapping]


def canonical_flow_grid(mesh_arrays_path,native_path):
    """Create a separately identified field artifact on the actual solver input mesh."""
    with np.load(mesh_arrays_path) as z:points=z['points_m'];tetra=z['tetra']
    native,velocity,pressure=parse_result_vtu(native_path)
    mapping,record=identify_native_grid(points,tetra,native)
    cells=np.column_stack([np.full(len(tetra),4),tetra])
    grid=pv.UnstructuredGrid(cells,np.full(len(tetra),10,dtype=np.uint8),points.copy())
    grid.point_data['Velocity']=np.asarray(velocity,dtype=float)[mapping]
    grid.point_data['Pressure']=np.asarray(pressure,dtype=float).reshape(-1)[mapping]
    for name in native.cell_data:grid.cell_data[name]=native.cell_data[name].copy()
    for name in native.field_data:grid.field_data[name]=native.field_data[name].copy()
    record.update(native_file_sha256=sha(native_path),mesh_arrays_sha256=sha(mesh_arrays_path),reader_source_sha256=sha(__file__),
        canonical_artifact='Declared input coordinates and connectivity with unmodified native nodal fields; native raw VTK preserved separately')
    assert np.array_equal(grid['Velocity'],np.asarray(velocity)[mapping])
    assert np.array_equal(grid['Pressure'],np.asarray(pressure).reshape(-1)[mapping])
    return grid,record
