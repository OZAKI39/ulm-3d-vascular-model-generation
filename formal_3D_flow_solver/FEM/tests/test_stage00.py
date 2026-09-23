"""Stage 0 acceptance checks; tests remain in the repository after success."""
import copy
import json
import math
import subprocess
from pathlib import Path
import pytest

from fem3d.audit import reference_snapshot,sha256,write_json
from fem3d.source import length_to_m,resolve_source_path,validate_contract

ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def contract():
    return json.loads((ROOT/"reports/stage00/source_contract.json").read_text())


@pytest.mark.parametrize("name",["configs","src/fem3d","scripts","tests","reports/stage00","outputs","logs","remote"])
def test_project_structure(name):
    assert (ROOT/name).is_dir()
    assert ROOT==Path("/home/lzy/projects/formal_3D_flow_solver/FEM")


def test_package_import_is_local():
    import fem3d
    assert Path(fem3d.__file__).resolve().is_relative_to(ROOT/"src")


def test_reference_trees_unchanged(tmp_path):
    baseline=json.loads((ROOT/"reports/stage00/reference_baseline.json").read_text())
    current=reference_snapshot(baseline["roots"])
    # Tests must not overwrite an accepted historical stage's evidence.
    write_json(tmp_path/"reference_final.json",current)
    old,new=baseline["files"],current["files"]
    changes=[k for k in old.keys()|new.keys() if old.get(k)!=new.get(k)]
    write_json(tmp_path/"reference_integrity.json",{
        "status":"FAIL" if changes else "PASS","entries":len(new),"changes":changes,
        "baseline_timestamp":baseline["timestamp"],"final_timestamp":current["timestamp"],
        "method":"Full SHA256, size, mode and mtime recomputation, with directory and symlink inventories"})
    assert changes==[]


def test_source_files_exist_and_hashes_recompute(contract):
    assert validate_contract(contract)
    for entry in contract["source_files"]:
        p=Path(entry["path"])
        assert p.is_file()
        assert sha256(p)==entry["sha256"]
    for name in ("source_metadata_files","source_configs"):
        assert all(Path(p).is_file() for p in contract[name])
    assert sha256(contract["geometry_path"])==contract["geometry_sha256"]
    assert sha256(contract["meter_geometry_path"])==contract["meter_geometry_sha256"]


def test_units_and_si_adapter(contract):
    assert contract["geometry_units"]=="um"
    assert contract["meter_geometry_units"]=="m"
    assert contract["solver_units"]=="SI"
    assert length_to_m(100.,"um")==pytest.approx(1e-4)
    with pytest.raises(ValueError):
        length_to_m(1,"unknown")


def test_ports_are_unique_finite_positive_and_explicit(contract):
    assert len(contract["inlets"])>=1 and len(contract["outlets"])>=1
    ports=contract["inlets"]+contract["outlets"]
    assert len({p["port_id"] for p in ports})==len(ports)
    assert len({p["surface_entity_id"] for p in ports})==len(ports)
    for p in ports:
        assert p["area_m2"]>0
        assert all(math.isfinite(x) for x in p["outward_normal"])
        assert sum(x*x for x in p["outward_normal"])>0
        assert p["upstream_normal_dot"]>.999
        assert p["role_source"] in contract["source_metadata_files"]


def test_tags_partition_source_surface_exactly(contract):
    import numpy as np
    import pyvista as pv
    surface=pv.read(contract["geometry_path"])
    tags=np.asarray(surface.cell_data[contract["tag_array"]])
    ports=contract["inlets"]+contract["outlets"]
    assert set(tags)=={contract["wall_entity_id"]}|{p["surface_entity_id"] for p in ports}
    assert sum(p["triangle_count"] for p in ports)+contract["wall_triangle_count"]==surface.n_cells
    for p in ports:
        assert (tags==p["surface_entity_id"]).sum()==p["triangle_count"]
    assert contract["geometry_checks"]["meter_stl_matches_scaled_tagged_triangles"]


@pytest.mark.parametrize("problem",["no_inlet","no_outlet","duplicate_port","zero_area","zero_normal","nan_normal","ambiguous_units"])
def test_bad_contracts_fail_closed(contract,problem):
    c=copy.deepcopy(contract)
    if problem=="no_inlet": c["inlets"]=[]
    if problem=="no_outlet": c["outlets"]=[]
    if problem=="duplicate_port": c["outlets"][0]["port_id"]=c["inlets"][0]["port_id"]
    if problem=="zero_area": c["inlets"][0]["area_m2"]=0
    if problem=="zero_normal": c["inlets"][0]["outward_normal"]=[0,0,0]
    if problem=="nan_normal": c["inlets"][0]["outward_normal"]=[float("nan"),0,1]
    if problem=="ambiguous_units": c["geometry_units"]="unknown"
    with pytest.raises(ValueError): validate_contract(c)


def test_historical_path_resolution_rejects_guessing(tmp_path):
    (tmp_path/"x").write_text("input")
    assert resolve_source_path("E:\\old\\x",tmp_path,"E:/old/")==tmp_path/"x"
    with pytest.raises(ValueError):resolve_source_path("D:\\unknown\\x",tmp_path,"E:/old/")
    with pytest.raises(ValueError):resolve_source_path("../outside",tmp_path,"E:/old/")
    with pytest.raises(FileNotFoundError):resolve_source_path("absent",tmp_path,"E:/old/")


def test_geometry_review_artifact(contract):
    from PIL import Image
    path=ROOT/"reports/stage00/source_geometry_overview.png"
    with Image.open(path) as im:
        assert im.width>=1600 and im.height>=900
        assert len(im.convert("RGB").getcolors(im.width*im.height))>100
    metadata=json.loads((ROOT/"reports/stage00/geometry_figure_metadata.json").read_text())
    assert metadata["geometry_sha256"]==contract["geometry_sha256"]
    assert metadata["image_sha256"]==sha256(path)
    assert not metadata["volume_mesh_generated"] and not metadata["flow_solved"]


@pytest.mark.parametrize("script",["remote_probe.sh","remote_sync.sh","remote_run.sh","remote_fetch.sh"])
def test_remote_shell_syntax_and_dry_run_interface(script):
    subprocess.run(["bash","-n",str(ROOT/"scripts"/script)],check=True)
    r=subprocess.run(["bash",str(ROOT/"scripts"/script),"--help"],capture_output=True,text=True,check=True)
    assert "--dry-run" in r.stdout


def test_remote_roundtrip_receipt():
    receipt=json.loads((ROOT/"outputs/stage00/remote_return/outputs/stage00/roundtrip_result.json").read_text())
    state=json.loads((ROOT/"remote/connection.local.json").read_text())
    assert receipt["status"]=="PASS"
    assert receipt["remote_work_dir"]==state["remote_work_dir"]
    assert receipt["source_of_truth"]==str(ROOT)
    for row in receipt["files"]:
        assert sha256(ROOT/row["project_path"])==row["sha256"]


@pytest.mark.parametrize("ranks",[1,2,4])
def test_remote_mpi_smoke_results(ranks):
    base=ROOT/"outputs/stage00/remote_return/outputs/stage00"
    result=json.loads((base/f"smoke_r{ranks}.json").read_text())
    launch=json.loads((base/f"smoke_r{ranks}_launch.json").read_text())
    assert result["status"]=="PASS" and result["mpi_ranks"]==ranks
    assert result["all_owned_solution_values_finite"] and result["ksp_converged_reason"]>0
    assert result["relative_algebraic_residual"]<1e-10
    assert result["dolfinx_version"].startswith("0.11.")
    assert launch["returncode"]==0 and launch["wall_time_s"]>0
