#!/usr/bin/env python3
"""Fill a frozen discrete boundary; run serially on the compute host."""
import argparse
import json
import os
import platform
import resource
import signal
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import gmsh
from fem3d.audit import environment, sha256, timestamp, write_json
from fem3d.mesh_input import FACET_NAMES, surface_topology
from fem3d.sparse_volume import audit_volume


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--candidate',required=True,choices=['sparse_A','sparse_B','sparse_C'])
    a=p.parse_args()
    policy_path=ROOT/'inputs/stage01_6/acceptance_policy.json'; policy=json.loads(policy_path.read_text())
    lock=json.loads((ROOT/'inputs/stage01_6/freeze_lock.json').read_text())
    assert sha256(policy_path)==lock['policy_sha256']
    contract_path=ROOT/'inputs/stage01_6/planar_port_contract_v2.json'
    assert sha256(contract_path)==lock['contract_sha256']
    contract=json.loads(contract_path.read_text())
    out=ROOT/'outputs/stage01_6'/a.candidate
    gate=json.loads((out/'qc/surface_invariants.json').read_text())
    assert gate['status']=='PASS' and gate['volume_meshing_permitted'] and all(gate['cap_gate']['checks'].values())
    for name in ('mesh','metadata','qc'): (out/name).mkdir(parents=True,exist_ok=True)
    if (out/'mesh/fluid.msh').exists(): raise RuntimeError('No mesh overwrite or automatic parameter search')
    source=out/'surface/tagged_surface_si.npz'
    sm=json.loads((out/'metadata/cap_remesh.json').read_text())
    assert sha256(source)==sm['derived_surface_sha256']
    assert sha256(ROOT/'inputs/stage01/tagged_surface_si.npz')==lock['source_surface_sha256']
    data=np.load(source);points,triangles,tags=data['points_m'],data['triangles'],data['facet_tags']
    topology=surface_topology(points,triangles)
    config={'limits':{'maximum_rss_gib':8,'maximum_wall_time_s':300,'maximum_tetrahedra':1500000},'gmsh':{'optimize_interior':True}}
    control={'bulk_target_size_m':5e-7}
    env = environment()
    metadata = {"run_id": os.environ.get("FEM3D_RUN_ID"), "timestamp": timestamp(),
                "hostname": platform.node(), "environment": env, "gpu_used": False,
                "source_surface_sha256": sha256(source),
                "source_contract_sha256": lock["source_contract_sha256"],
                "configuration_sha256": sha256(policy_path),
                "configuration": policy, "requested_mesh_controls": control,
                "gmsh_version": gmsh.__version__, "python_version": sys.version,
                "units": {"coordinates": "m", "length": "m", "volume": "m^3"},
                "mesh_profile": a.candidate, "facet_tag_meanings": FACET_NAMES,
                "cell_tag_meanings": {100: "FLUID"}, "threads": 1, "mpi_ranks": 1,
                "stdout_path": os.environ.get("FEM3D_STDOUT_LOG"),
                "stderr_path": os.environ.get("FEM3D_STDERR_LOG"),
                "run_metadata_path": os.environ.get("FEM3D_RUN_METADATA"),
                "input_topology": topology, "exit_status": None}
    start = time.perf_counter()
    # Bound this development experiment. No automatic increase in resource limits.
    resource.setrlimit(resource.RLIMIT_AS, (int(config["limits"]["maximum_rss_gib"] * 1024**3),) * 2)
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("Meshing wall-time limit")))
    signal.alarm(config["limits"]["maximum_wall_time_s"])
    try:
        gmsh.initialize()
        gmsh.model.add(a.candidate)
        options = {"General.NumThreads": 1, "Mesh.MaxNumThreads1D": 1,
                   "Mesh.MaxNumThreads2D": 1, "Mesh.MaxNumThreads3D": 1,
                   "Mesh.Algorithm3D": 1, "Mesh.MeshOnlyEmpty": 1,
                   "Mesh.MeshSizeFromPoints": 0, "Mesh.MeshSizeFromCurvature": 0,
                   "Mesh.MeshSizeExtendFromBoundary": 0,
                   "Mesh.MeshSizeMin": control["bulk_target_size_m"],
                   "Mesh.MeshSizeMax": control["bulk_target_size_m"],
                   "Mesh.Optimize": int(config["gmsh"]["optimize_interior"]),
                   "Mesh.OptimizeNetgen": 0, "Mesh.ElementOrder": 1,
                   "Mesh.MshFileVersion": 4.1, "Mesh.Binary": 1}
        for key, value in options.items():
            gmsh.option.setNumber(key, value)
        assert options==policy["volume_meshing"]["effective_gmsh_options"]
        metadata["effective_gmsh_options"] = options
        for tag in FACET_NAMES:
            gmsh.model.addDiscreteEntity(2, tag)
        gmsh.model.mesh.addNodes(2, 1, np.arange(1, len(points)+1), points.ravel())
        for tag in FACET_NAMES:
            indices = np.flatnonzero(tags == tag)
            gmsh.model.mesh.addElementsByType(tag, 2, indices + 1, (triangles[indices]+1).ravel())
        gmsh.model.mesh.reclassifyNodes()
        loop = gmsh.model.geo.addSurfaceLoop(list(FACET_NAMES))
        volume = gmsh.model.geo.addVolume([loop])
        gmsh.model.geo.synchronize()
        for tag, name in FACET_NAMES.items():
            gmsh.model.addPhysicalGroup(2, [tag], tag, name)
        gmsh.model.addPhysicalGroup(3, [volume], 100, "FLUID")
        metadata["discrete_topology"] = {"surface_entities": gmsh.model.getEntities(2),
            "volume_entities": gmsh.model.getEntities(3),
            "volume_boundary": gmsh.model.getBoundary([(3, volume)], oriented=True),
            "method": "5 premeshed discrete surfaces; one built-in surface loop and volume; MeshOnlyEmpty=1; no surface parametrization/remeshing"}
        write_json(out / "metadata/meshing.json", metadata)
        gmsh.model.mesh.generate(3)
        types, _, _ = gmsh.model.mesh.getElements(3)
        if list(types) != [4]:
            raise RuntimeError(f"Expected only linear tetrahedra, got {types}")
        element_ids, connectivity = gmsh.model.mesh.getElementsByType(4)
        if len(element_ids) > config["limits"]["maximum_tetrahedra"]:
            raise RuntimeError(f"Tetra count {len(element_ids)} exceeds development limit; no escalation")
        node_ids, coords, _ = gmsh.model.mesh.getNodes()
        order = np.argsort(node_ids)
        node_ids, coords = node_ids[order], coords.reshape(-1, 3)[order]
        tetra = np.searchsorted(node_ids, connectivity).reshape(-1, 4)
        btri, btags = [], []
        for tag in FACET_NAMES:
            _, indices = gmsh.model.mesh.getElementsByType(2, tag)
            btri.append(np.searchsorted(node_ids, indices).reshape(-1, 3))
            btags.extend([tag] * (len(indices)//3))
        boundary = np.concatenate(btri)
        quality = gmsh.model.mesh.getElementQualities(element_ids, "minSICN")
        gmsh.write(str(out / "mesh/fluid.msh"))
        np.savez_compressed(out / "mesh/volume_mesh.npz", points_m=coords,
            tetra=tetra, boundary_triangles=boundary, facet_tags=np.array(btags, dtype=np.int32),
            cell_tags=np.full(len(tetra), 100, dtype=np.int32),
            gmsh_element_ids=element_ids, min_sicn=quality)
        metadata.update(actual_tetra_count=len(tetra), actual_vertex_count=len(coords),
            boundary_facet_count=len(boundary),
            boundary_counts={name: int(np.count_nonzero(np.array(btags)==tag)) for tag, name in FACET_NAMES.items()},
            physical_groups=[{"dim": d, "tag": t, "name": gmsh.model.getPhysicalName(d,t)} for d,t in gmsh.model.getPhysicalGroups()],
            exit_status=0)
    except BaseException as exc:
        metadata.update(exit_status=1, error=repr(exc))
        raise
    finally:
        signal.alarm(0)
        metadata.update(wall_time_s=time.perf_counter()-start,
                        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        memory_limit_scope="8 GiB virtual address-space upper bound (stricter than RSS)")
        write_json(out / "metadata/meshing.json", metadata)
        if gmsh.isInitialized():
            gmsh.finalize()
    mesh_path=out/'mesh/volume_mesh.npz'
    result=audit_volume(np.load(mesh_path),data,np.load(ROOT/'inputs/stage01/tagged_surface_si.npz'),contract,policy)
    result.update(timestamp=timestamp(),mesh_sha256=sha256(mesh_path),candidate=a.candidate,policy_sha256=sha256(policy_path))
    write_json(out/'qc/volume_quality.json',result)
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
