import json
from fem3d.audit import sha256


def test_fresh_process_save_and_mpi_reload_keep_mesh_and_tags(stage01_mesh):
    base=stage01_mesh["base"]
    save=json.loads((base/"metadata/dolfinx_save.json").read_text())
    qc=stage01_mesh["qc"]
    run_ids=[save["run_id"]]
    for ranks in (1,2):
        loaded=json.loads((base/f"qc/reload_r{ranks}.json").read_text())
        assert loaded["status"]=="PASS" and loaded["mpi_ranks"]==ranks
        assert loaded["gdim"]==loaded["tdim"]==3
        assert loaded["cell_type"]=="tetrahedron"
        assert loaded["cell_count"]==qc["tetrahedron_count"]
        assert loaded["vertex_count"]==qc["vertex_count"]
        assert loaded["cell_tags"]==qc["cell_tags"]
        assert loaded["topology"]["boundary_counts"]==qc["topology"]["boundary_counts"]
        assert loaded["tetra_geometry_exactly_equal"] and loaded["file_closed"]
        assert loaded["boundary_fidelity"]["maximum_boundary_displacement_m"]==0
        assert sum(row["owned_cells"] for row in loaded["ranks"])==qc["tetrahedron_count"]
        for ext,key in (("xdmf","xdmf_sha256"),("h5","hdf5_sha256")):
            assert loaded[key]==save[key]==sha256(base/f"mesh/fluid.{ext}")
        run_ids.append(loaded["run_id"])
    assert len(set(run_ids))==3
    for run_id in run_ids:
        log=stage01_mesh["root"]/"outputs/stage01/remote_return/logs/stage01"/run_id
        metadata=json.loads((log/"metadata.json").read_text())
        assert metadata["returncode"]==0 and metadata["finished_utc"]
        assert (log/"stdout.log").is_file() and (log/"stderr.log").is_file()
