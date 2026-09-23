from stage017_helpers import *
import pytest
from fem3d.adaptive_qc import measure_volume
from fem3d.audit import sha256

def test_baseline_recomputed_from_raw_artifacts_matches_history():
 d=measure_volume(dict(np.load(ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz')),source(),source(),CONTRACT,POLICY)
 h=read(ROOT/'outputs/stage01/medium/qc/geometry_qc.json')
 assert d['proxy']==BASE['proxy']
 assert d['quality']['min_sicn']==h['quality']['gmsh_min_sicn']
 assert d['quality']['low_quality_nearest_boundary_counts']==h['quality']['low_quality_nearest_boundary_counts']
 assert d['quality']['total_below_0_1']==h['quality']['advisory_count']
 assert sha256(ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz')==BASE['mesh_sha256']
 assert all(BASE['history_checks'].values())
