"""Read-only velocity representation inventory; no flow or divergence integration."""
from pathlib import Path
import argparse
import hashlib
import json
import struct
import time

import numpy as np
import pyvista as pv
from particle_3d.audit import read_frozen
from particle_3d.field import FrozenFEMField


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def array_sha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def compare(a, b):
    a, b = np.asarray(a), np.asarray(b)
    diff = a.astype(float)-b.astype(float)
    scale = np.linalg.norm(a.astype(float))
    return dict(numeric_equal=bool(np.array_equal(a, b)),
                bitwise_equal=bool(a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()),
                source_dtype=str(a.dtype), destination_dtype=str(b.dtype),
                max_absolute_difference=float(np.max(np.abs(diff))),
                rms_difference=float(np.sqrt(np.mean(diff**2))),
                relative_L2_difference=float(np.linalg.norm(diff)/scale) if scale else 0.)


def checkpoint(path, nodes):
    raw = Path(path).read_bytes()
    h = struct.unpack_from('<8i3d', raw)
    if h[:8] != (1, 1, 1, nodes, 0, 4, 0, 71) or len(raw) != 56+64*nodes:
        raise ValueError('Not the audited single-rank, four-DOF, step-71 checkpoint layout')
    return np.frombuffer(raw, dtype='<f8', count=4*nodes, offset=56).reshape(nodes, 4), h


def run(repo, case, source):
    report = repo/'particle_3d/reports/particle9a3b_flowfield_conservation'
    data = report/'data'
    data.mkdir(parents=True, exist_ok=True)
    froot = repo/'formal_3D_flow_solver/FEM_SimVascular'
    summary, mesh, particle, boundaries = read_frozen(froot)
    field = FrozenFEMField.from_grids(mesh, particle)
    native_path = case/'run/1-procs/result_071.vtu'
    frozen_path = case/'frozen_flow/steady_flow_mean_2p0_mmps.vtu'
    native, exported = pv.read(native_path), pv.read(frozen_path)
    cp_path = case/'run/1-procs/stFile_071.bin'
    cp, header = checkpoint(cp_path, mesh.n_points)
    native_mesh = pv.read(case/'SV_MESH/mesh-complete.mesh.vtu')
    assert np.all(native_mesh.celltypes == 10)
    assert np.array_equal(native_mesh.points, native.points)
    entries = []
    evidence = ['run/solver.xml:Add_equation fluid; no Use_taylor_hood_type_basis',
                'read_files.cpp:1742-1776 (THflag=false => nFs=1)',
                'fs.cpp:300-317 (first space inherits mesh)',
                'nn.cpp:87,174-184 (TET4 => Tetra4)',
                'FE/Basis/LagrangeBasis.cpp:30-62 (fixed linear order)',
                'fluid.cpp:496-499,1534,1727-1728 (VMS weak continuity)']

    def entry(layer, path, points, tetra, u, storage, meaning):
        entries.append(dict(layer=layer, path=str(path), SHA256=sha(path) if Path(path).is_file() else None,
            point_count=len(points), cell_count=len(tetra), velocity_array_name='Velocity' if storage=='point data' else 'Y_n[0:3,node]',
            velocity_storage=storage, component_count=3, dtype=str(u.dtype), coordinate_dtype=str(points.dtype),
            velocity_dtype=str(u.dtype), finite_element_order=1, finite_element_space='continuous nodal TET4 P1 velocity; P1 pressure; VMS stabilized',
            order_evidence=evidence, representation=meaning,
            coordinate_array_SHA256=array_sha(points), connectivity_array_SHA256=array_sha(tetra), velocity_array_SHA256=array_sha(u)))

    nt = native.cells.reshape(-1, 5)[:, 1:]
    entry('1_native_solution', cp_path, native_mesh.points.astype(float), nt, cp[:, :3], 'other: raw original FE DOFs',
          'Original resolved nodal FE solution, recoverable from checkpoint; no independently stored higher-order transport field')
    entry('2_native_checkpoint', cp_path, native_mesh.points.astype(float), nt, cp[:, :3], 'other: binary column-major Array<double>(4,N)',
          '56-byte <8i3d header, Y_n then A_n; coordinates/connectivity supplied by authoritative solver input mesh')
    for layer, path, grid in [('3_exported_native_VTU', native_path, native), ('4_frozen_flow_VTU', frozen_path, exported),
                              ('5_particle_frozen_reference', Path(summary['flow_path']), particle)]:
        entry(layer, path, grid.points, grid.cells.reshape(-1, 5)[:, 1:], grid['Velocity'], 'point data', 'Original nodal DOFs serialized; no nodal projection')
    entry('6_FrozenFEMField_internal_arrays', Path(summary['flow_path']), field.points_m, field.tetra, field.velocity_nodes_m_s,
          'other: in-memory float64 arrays', 'Coordinates cast float32 to float64 exactly; canonical connectivity; velocity copied without resampling')
    entry('7_Particle_tetra_local_reconstruction', repo/'particle_3d/src/particle_3d/field.py', field.points_m, field.tetra,
          field.velocity_nodes_m_s, 'other: barycentric evaluation of nodal arrays', 'u(x)=sum_i N_i(x)u_i; affine P1 inside each tetra; no added degrees of freedom')
    entries[0]['coordinate_note'] = 'Solver reads float32 mesh coordinates into double. This does not add significant bits.'
    write(data/'velocity_representation_inventory.json', dict(completed_before_divergence_and_flux_computation=True,
          utc_unix_s=time.time(), distinct_physical_velocity_fields=1, layers=entries,
          note='Seven pipeline layers, one resolved P1 field. Native VMS fine scales are residual-based quadrature quantities, not another saved higher-order nodal velocity.'))
    write(data/'native_to_export_audit.json', dict(checkpoint_header=header, checkpoint_sha256=sha(cp_path),
          checkpoint_layout_evidence=['output.cpp:200-204,296-299,401-416','Array.h:346-359'],
          checkpoint_to_native_VTU_velocity=compare(cp[:, :3], native['Velocity']),
          checkpoint_to_native_VTU_pressure=compare(cp[:, 3], native['Pressure']),
          native_input_to_output_coordinates=compare(native_mesh.points, native.points),
          native_to_frozen_file_bitwise_equal=sha(native_path)==sha(frozen_path),
          conversion_chain=['output.cpp writes current Y_n as double', 'vtk_xml.cpp:925-927,1100-1107 copies current Y by node index',
            'VtkData.cpp:166-190 writes vtkDoubleArray; coordinates vtkPoints at 240-261',
            'scripts/flow_2mmps/validate.py:97 uses shutil.copyfile(final,frozen)'],
          interpolation=False, projection=False, averaging=False, velocity_float_conversion=False,
          coordinate_float_conversion='double -> float32 in VTK writer; input mesh was already float32; observed numeric change=0',
          mapping='single rank; checkpoint Y_n node rows match all VTU velocity and pressure values exactly'))
    write(data/'export_to_particle_audit.json', dict(export_file_sha256=sha(frozen_path), particle_file_sha256=sha(summary['flow_path']),
          exported_to_reference_file_bitwise_equal=sha(frozen_path)==sha(summary['flow_path']),
          coordinates_export_to_reference=compare(exported.points, particle.points),
          connectivity_export_to_reference=compare(exported.cells, particle.cells),
          velocity_export_to_reference=compare(exported['Velocity'], particle['Velocity']),
          coordinates_reference_to_internal=compare(particle.points, field.points_m),
          velocity_reference_to_internal=compare(particle['Velocity'], field.velocity_nodes_m_s),
          flow_to_canonical_local_vertex_permutation=summary['local_vertex_permutation'],
          connectivity_same_node_sets_per_cell=bool(np.array_equal(np.sort(nt,axis=1),np.sort(field.tetra,axis=1))),
          source=['audit.py:39-50,85-89', 'field.py:42-57,76-80,115-122', 'geometry.py:7-10,26-59'],
          resampling=False, coordinate_transformation=False))
    refs = ['read_files.cpp','fs.cpp','nn.cpp','fluid.cpp','output.cpp','vtk_xml.cpp','VtkData.cpp','Array.h',
            'FE/Basis/LagrangeBasis.cpp','FE/Basis/BasisTraits.h','distribute.cpp','vtk_xml_parser.cpp']
    write(data/'native_source_hashes.json', {name:dict(path=str(source/name), sha256=sha(source/name)) for name in refs})
    with np.load(case/'SV_MESH/mesh_arrays.npz') as arrays:
        coordinate_audit = compare(native_mesh.points, arrays['points_m'])
    write(data/'precision_audit.json', dict(coordinate_storage='float32', velocity_storage='float64',
          solver_input_equals_export_coordinates=True, mesh_npz_vs_input=coordinate_audit,
          authoritative_higher_precision_solver_coordinates_available=False,
          authoritative_precision_loss_test='cannot test authoritative precision loss: solver input, NPZ and all exported coordinates have exactly the same quantized values',
          velocity_checkpoint_to_export=compare(cp[:, :3], native['Velocity']),
          float32_to_float64='exact promotion of stored values; no resampling and no recovery of lost coordinate precision'))
    print(json.dumps(dict(inventory_layers=len(entries), points=mesh.n_points, tetra=mesh.n_cells,
          checkpoint_velocity_exact=np.array_equal(cp[:, :3], native['Velocity']), frozen_copy_exact=sha(native_path)==sha(frozen_path))))


if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('--repo',type=Path,required=True);p.add_argument('--case',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True);a=p.parse_args();run(a.repo,a.case,a.source)
