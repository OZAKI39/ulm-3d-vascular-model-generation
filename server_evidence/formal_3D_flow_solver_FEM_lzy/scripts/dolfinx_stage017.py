#!/usr/bin/env python3
"""Save or reload mesh/tags in separate processes; no FEM function spaces."""
import argparse
import json
import os
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
import numpy as np
from mpi4py import MPI
import dolfinx
from dolfinx import io, mesh
from dolfinx.io import gmsh as gmsh_io
from scipy.spatial import cKDTree
from fem3d.audit import sha256, timestamp, write_json
from fem3d.mesh_qc import boundary_partition, tetra_volumes
from fem3d.adaptive_qc import measure_volume


def main():
    p=argparse.ArgumentParser()
    p.add_argument("action",choices=("save","reload"))
    a=p.parse_args()
    comm=MPI.COMM_WORLD
    out=ROOT/'outputs/stage01_7/selected'
    policy=json.loads((ROOT/'inputs/stage01_7/acceptance_policy.json').read_text())
    lock=json.loads((ROOT/'inputs/stage01_7/freeze_lock.json').read_text())
    assert out.exists(),'Only a selected feasible iteration may roundtrip'
    assert sha256(out/'planar_port_contract_v2.json')==lock['contract_sha256']
    assert sha256(ROOT/'inputs/stage01_7/acceptance_policy.json')==lock['policy_sha256']
    path=out/"mesh/fluid.xdmf"
    if a.action=="save":
        if comm.size!=1:
            raise RuntimeError("Save once with one rank; reload with 1/2 ranks independently")
        if path.exists():
            raise RuntimeError("Refusing to overwrite saved mesh")
        data=gmsh_io.read_from_msh(str(out/"mesh/fluid.msh"),comm,rank=0,gdim=3)
        msh, ct, ft=data.mesh,data.cell_tags,data.facet_tags
        msh.name="fluid_mesh"
        ct.name="cell_tags"
        ft.name="facet_tags"
        msh.topology.create_connectivity(2,3)
        with io.XDMFFile(comm,str(path),"w") as file:
            file.write_mesh(msh)
            file.write_meshtags(ct,msh.geometry)
            file.write_meshtags(ft,msh.geometry)
        write_json(out/"metadata/dolfinx_save.json", {"status":"PASS","timestamp":timestamp(),
            "pid":os.getpid(),"run_id":os.environ.get("FEM3D_RUN_ID"),"dolfinx_version":dolfinx.__version__,
            "mpi_ranks":comm.size,"mesh_name":msh.name,"cell_tag_name":ct.name,"facet_tag_name":ft.name,
            "physical_groups":{k:list(v) for k,v in data.physical_groups.items()},
            "xdmf_sha256":sha256(path),"hdf5_sha256":sha256(path.with_suffix(".h5")),
            "artifact_format":"XDMF/HDF5", "file_closed":True})
        print("DOLFINx mesh and tags saved; exiting process.")
        return
    with io.XDMFFile(comm,str(path),"r") as file:
        msh=file.read_mesh(name="fluid_mesh",ghost_mode=mesh.GhostMode.shared_facet)
        msh.topology.create_connectivity(2,3)
        msh.topology.create_connectivity(3,2)
        ct=file.read_meshtags(msh,name="cell_tags")
        ft=file.read_meshtags(msh,name="facet_tags")
    assert msh.geometry.dim==3 and msh.topology.dim==3
    assert msh.topology.cell_type==mesh.CellType.tetrahedron
    nc=msh.topology.index_map(3).size_local
    nf=msh.topology.index_map(2).size_local
    owned_cells=ct.indices<nc
    # Only owned cells/facets participate in global counts; never double-count ghosts.
    assert np.array_equal(ct.indices[owned_cells],np.arange(nc))
    assert np.all(ct.values[owned_cells]==100)
    owned_facets=ft.indices<nf
    ext=mesh.exterior_facet_indices(msh.topology)
    ext=ext[ext<nf]
    assert np.array_equal(np.sort(ext),np.sort(ft.indices[owned_facets]))
    assert len(np.unique(ft.indices[owned_facets]))==np.count_nonzero(owned_facets)
    fgeom=mesh.entities_to_geometry(msh,2,ft.indices[owned_facets],permute=False)
    cells=msh.geometry.x[msh.geometry.dofmap[:nc]]
    facets=msh.geometry.x[fgeom]
    gathered=comm.gather((cells,facets,ft.values[owned_facets]),root=0)
    local_counts=comm.gather({"rank":comm.rank,"owned_cells":nc,"owned_tagged_exterior_facets":len(facets),
                              "ghost_cells":msh.topology.index_map(3).num_ghosts,"pid":os.getpid()},root=0)
    if comm.rank==0:
        coordinates=np.concatenate([g[0] for g in gathered])
        points,inverse=np.unique(coordinates.reshape(-1,3),axis=0,return_inverse=True)
        tetra=inverse.reshape(-1,4)
        final_facets=np.concatenate([g[1] for g in gathered])
        distances,indices=cKDTree(points).query(final_facets.reshape(-1,3))
        assert np.all(distances==0)
        triangles=indices.reshape(-1,3)
        tags=np.concatenate([g[2] for g in gathered])
        raw=np.load(out/"mesh/volume_mesh.npz")
        dist,raw_indices=cKDTree(raw["points_m"]).query(points)
        assert np.all(dist==0) and len(points)==len(raw["points_m"])
        def canonical(array):
            s=np.sort(array,axis=1)
            return s[np.lexsort(s[:,::-1].T)]
        assert np.array_equal(canonical(raw_indices[tetra]),canonical(raw["tetra"]))
        exterior, exterior_tags, partition=boundary_partition(points,tetra,triangles,tags)
        source=np.load(out/'surface/tagged_surface_si.npz')
        original=np.load(ROOT/'inputs/stage01/tagged_surface_si.npz')
        contract=json.loads((out/'planar_port_contract_v2.json').read_text())
        # Re-evaluate minSICN on the loaded tetra geometry. No mesh generation.
        positive=tetra.copy();signed=tetra_volumes(points,positive);flip=signed<0
        positive[flip]=positive[flip][:,[1,0,2,3]]
        import gmsh
        gmsh.initialize()
        try:
            gmsh.model.add('reload_quality_only');gmsh.model.addDiscreteEntity(3,100)
            gmsh.model.mesh.addNodes(3,100,np.arange(1,len(points)+1),points.ravel())
            eid=np.arange(1,len(positive)+1);gmsh.model.mesh.addElementsByType(100,4,eid,(positive+1).ravel())
            quality=np.asarray(gmsh.model.mesh.getElementQualities(eid,'minSICN'))
        finally:gmsh.finalize()
        measured=measure_volume({'points_m':points,'tetra':positive,'boundary_triangles':triangles,'facet_tags':tags,'cell_tags':np.full(len(tetra),100),'gmsh_element_ids':eid,'min_sicn':quality},source,original,contract,policy)
        stored=json.loads((out/'qc/volume_quality.json').read_text())
        for key in ('minimum','P1','P5','median','P95','maximum'):
            assert abs(measured['quality']['min_sicn'][key]-stored['quality']['min_sicn'][key])<=1e-12
        assert measured['quality']['low_quality_nearest_boundary_counts']==stored['quality']['low_quality_nearest_boundary_counts']
        assert measured['proxy']==stored['proxy']
        boundary,ports=measured['boundary_fidelity'],measured['ports']
        volumes=np.abs(tetra_volumes(points,tetra));assert np.isfinite(volumes).all() and np.all(volumes>0)
        assert partition['connected_fluid_components']==1
        saved=json.loads((out/'metadata/dolfinx_save.json').read_text())
        assert all(row['pid']!=saved['pid'] for row in local_counts)
        assert sha256(path)==saved['xdmf_sha256'] and sha256(path.with_suffix('.h5'))==saved['hdf5_sha256']
        report={"status":"PASS","timestamp":timestamp(),"run_id":os.environ.get("FEM3D_RUN_ID"),
            "dolfinx_version":dolfinx.__version__,"mpi_ranks":comm.size,"ranks":local_counts,
            "gdim":3,"tdim":3,"cell_type":"tetrahedron","cell_count":len(tetra),"vertex_count":len(points),
            "cell_tags":{"100":len(tetra)},"topology":partition,"boundary_fidelity":boundary,"ports":ports,
            "tetra_geometry_exactly_equal":True,"positive_geometric_volume":True,
            "volume_m3":float(volumes.sum()),
            "orientation_note":"DOLFINx may permute tetra vertex ordering. Positivity of Gmsh's oriented tetrahedra was checked before export; reload verifies identical unordered tetrahedra and positive geometric volume.",
            "xdmf_sha256":sha256(path),"hdf5_sha256":sha256(path.with_suffix(".h5")),
            "quality_summary_recomputed":measured["quality"],"quality_summary_matches":True,"no_mesh_generation_in_reload":True,"fem_space_created":False,"file_closed":True}
        write_json(out/f"qc/reload_r{comm.size}.json",report)
        print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
