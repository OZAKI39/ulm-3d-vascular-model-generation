"""Frozen P1 field loading: original node/cell identities, canonical gradients."""
from dataclasses import dataclass
import numpy as np
import pyvista as pv
from .paths import ROOT, local_path, read_json, sha256
from flow_solver_support.wss_case import coordinate_identity
from vascular_validation.postprocess import SolutionMeasurements
from particle_3d.field import FrozenFEMField

@dataclass
class FrozenFlow:
    mesh: object
    grid: object
    arrays: dict
    velocity: np.ndarray
    pressure: np.ndarray
    policy: dict
    record: dict
    correspondence: dict

    def sampler(self):
        return FrozenFEMField.from_grids(self.mesh, self.grid)

def load(case='balanced_pressure_final', root=ROOT):
    record=read_json(root/'configs/cfd/cases.json')[case]
    for name, expected in record['sha256'].items():
        assert sha256(local_path(record[name],root))==expected, name
    mesh=pv.read(local_path(record['mesh'],root));flow=pv.read(local_path(record['flow'],root))
    correspondence=coordinate_identity(flow.points,mesh.points)
    canonical=mesh.cells_dict[pv.CellType.TETRA];native=flow.cells_dict[pv.CellType.TETRA]
    assert canonical.shape==native.shape and np.array_equal(np.sort(canonical,axis=1),np.sort(native,axis=1))
    correspondence.update(tetra_rows_same_members=True,local_order_changed_cells=int(np.any(canonical!=native,axis=1).sum()))
    with np.load(local_path(record['mesh_arrays'],root)) as a: arrays={k:a[k] for k in a.files}
    assert np.array_equal(arrays['tetra'],canonical) and np.array_equal(arrays['points_m'],mesh.points)
    policy=read_json(local_path(record['policy'],root))
    measure=SolutionMeasurements(local_path(record['mesh_arrays'],root),policy['Q_target_m3_s'],policy['Umean_m_s'])
    velocity,pressure=measure.read(local_path(record['flow'],root))
    assert np.isfinite(velocity).all() and np.isfinite(pressure).all()
    return FrozenFlow(mesh,flow,arrays,velocity,pressure,policy,record,correspondence)

def measurements(field, root=ROOT):
    measure=SolutionMeasurements(local_path(field.record['mesh_arrays'],root),field.policy['Q_target_m3_s'],field.policy['Umean_m_s'])
    return measure.measure(field.velocity,field.pressure)

