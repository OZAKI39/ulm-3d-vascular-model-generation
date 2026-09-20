import numpy as np
import pyvista as pv
from fem_sync_support import ROOT, read, digest


def test_volume_is_stage_q_input_not_old_fem_mesh():
    m=read('frozen_reference/mesh_manifest.json')
    expected=read('reports/sv1_3q/remote/final_input_integrity.json')['files'][0]
    assert digest(ROOT/m['path'])==m['sha256']==expected['sha256']
    g=pv.read(ROOT/m['path'])
    assert g.n_points==m['nodes'] and g.n_cells==m['tetra']
    assert np.all(g.celltypes==10) and m['units']=='m'
    assert list(g.bounds)==m['bounding_box_m']


def test_mesh_index_and_flow_order_contract_explicit():
    m=read('frozen_reference/mesh_manifest.json')
    assert m['flow_points_identical'] and m['flow_tetra_sets_identical']
    assert not m['flow_cell_connectivity_order_identical']
    assert m['flow_cell_rows_identical_ignoring_local_vertex_order']
    assert sorted(m['flow_local_vertex_permutation_from_volume'])==[0,1,2,3]
    assert 'zero-based' in m['tetra_id_convention'] and 'GlobalElementID' in m['tetra_id_convention']
