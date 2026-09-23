"""Independent synthetic analytic pipe. Never opens a vascular mesh."""
import json
from pathlib import Path
import numpy as np
from .audit import sha256, timestamp, write_json
from .analytic import pipe_frame


def load_benchmark_config(root):
    root=Path(root)
    bundle=json.loads((root/"inputs/stage02/pipe_config.json").read_text())
    if bundle["yaml_sha256"]!=sha256(root/"configs/stage02_pipe_benchmark.yaml"):
        raise ValueError("Benchmark YAML and remote JSON snapshot disagree")
    config=bundle["config"]
    assert config["units"]=="SI" and config["experimental_condition"] is False
    assert set(config["profiles"])=={"pipe_coarse","pipe_medium","pipe_fine"}
    assert config["boundary_tags"]=={"WALL":1,"OUTLET":2,"INLET":4}
    return config,bundle["yaml_sha256"]


def build_pipe(root,profile):
    import os
    import platform
    import resource
    import time
    import gmsh
    import dolfinx
    from dolfinx import fem,io
    from dolfinx.io import gmsh as gmsh_io
    from mpi4py import MPI
    import ufl
    root=Path(root)
    config,config_hash=load_benchmark_config(root)
    assert MPI.COMM_WORLD.size==1, "Gmsh meshing is serial"
    for ranks in (1,2):
        api=json.loads((root/f"outputs/stage02/api/real_r{ranks}.json").read_text())
        assert api["status"]=="PASS" and api["C_nonzero"] and api["spaces"]["real_global_dofs"]==1
    gate=json.loads((root/"inputs/stage02/api_block_gate.json").read_text())
    assert gate["status"]=="PASS"
    geometry=config["geometry"]
    origin,axis,_,_=pipe_frame(geometry["origin_m"],geometry["axis"])
    radius,length=geometry["radius_m"],geometry["length_m"]
    size=config["profiles"][profile]["target_size_m"]
    base=root/"outputs/stage02/meshes"/profile
    base.mkdir(parents=True,exist_ok=True)
    mshpath=base/"pipe.msh"
    if mshpath.exists():
        raise RuntimeError("Refusing to overwrite canonical mesh")
    start=time.perf_counter()
    gmsh.initialize()
    try:
        gmsh.model.add(profile)
        for key,value in {"General.NumThreads":1,"Mesh.MaxNumThreads1D":1,"Mesh.MaxNumThreads2D":1,
            "Mesh.MaxNumThreads3D":1,"Mesh.MeshSizeMin":size,"Mesh.MeshSizeMax":size,
            "Mesh.Algorithm3D":1,"Mesh.MshFileVersion":4.1,"Mesh.Binary":1}.items():
            gmsh.option.setNumber(key,value)
        volume=gmsh.model.occ.addCylinder(*origin,*(axis*length),radius)
        gmsh.model.occ.synchronize()
        surfaces={}
        for dim,tag in gmsh.model.getBoundary([(3,volume)],oriented=False):
            center=np.array(gmsh.model.occ.getCenterOfMass(dim,tag))
            axial=float((center-origin)@axis)
            physical=4 if abs(axial)<length*1e-9 else (2 if abs(axial-length)<length*1e-9 else 1)
            assert physical not in surfaces
            surfaces[physical]=tag
        assert set(surfaces)=={1,2,4}
        for name,physical in config["boundary_tags"].items():
            gmsh.model.addPhysicalGroup(2,[surfaces[physical]],physical,name)
        gmsh.model.addPhysicalGroup(3,[volume],100,"FLUID")
        gmsh.model.mesh.generate(3)
        cell_ids,_=gmsh.model.mesh.getElementsByType(4)
        if len(cell_ids)>config["limits"]["maximum_tetrahedra"]:
            raise RuntimeError("Canonical mesh exceeds frozen development cell limit")
        gmsh.model.mesh.setOrder(geometry["coordinate_degree"])
        types,blocks,_=gmsh.model.mesh.getElements(3)
        assert list(types)==[11], "Quadratic geometry tetrahedra expected"
        determinant=gmsh.model.mesh.getElementQualities(blocks[0],"minDetJac")
        assert np.isfinite(determinant).all() and np.all(determinant>0)
        gmsh.write(str(mshpath))
        imported=gmsh_io.model_to_mesh(gmsh.model,MPI.COMM_WORLD,0,gdim=3)
        domain,ct,ft=imported.mesh,imported.cell_tags,imported.facet_tags
        domain.name="pipe_mesh"; ct.name="cell_tags"; ft.name="facet_tags"
        domain.topology.create_connectivity(2,3)
        with io.XDMFFile(domain.comm,str(base/"pipe.xdmf"),"w") as file:
            file.write_mesh(domain)
            file.write_meshtags(ct,domain.geometry)
            file.write_meshtags(ft,domain.geometry)
        ds=ufl.Measure("ds",domain=domain,subdomain_data=ft,metadata={"quadrature_degree":6})
        dx=ufl.Measure("dx",domain=domain,metadata={"quadrature_degree":6})
        cache={"cache_dir":str(root/"outputs/stage02/jit_cache")}
        def integrate(expr):
            return float(fem.assemble_scalar(fem.form(expr,jit_options=cache)))
        areas={name:integrate(1*ds(tag)) for name,tag in config["boundary_tags"].items()}
        measured_volume=integrate(1*dx)
        normal=ufl.FacetNormal(domain)
        axis_ufl=ufl.as_vector(axis)
        normal_checks={name:integrate(ufl.dot(normal,axis_ufl)*ds(tag))/areas[name] for name,tag in {"INLET":4,"OUTLET":2}.items()}
        points=domain.geometry.x
        axial=(points-origin)@axis
        radial=np.linalg.norm(points-origin-axial[:,None]*axis,axis=1)
        relative_area={name:abs(areas[name]-(2*np.pi*radius*length if name=="WALL" else np.pi*radius**2))/(2*np.pi*radius*length if name=="WALL" else np.pi*radius**2) for name in areas}
        relative_volume=abs(measured_volume-np.pi*radius**2*length)/(np.pi*radius**2*length)
        assert abs(axial.min())<length*1e-10 and abs(axial.max()-length)<length*1e-10
        assert abs(radial.max()-radius)<radius*1e-10
        assert max(relative_area.values())<1e-3 and relative_volume<1e-3
        assert normal_checks["INLET"]<-.999999 and normal_checks["OUTLET"]>.999999
        metadata={"status":"PASS","purpose":config["purpose"],"experimental_condition":False,"profile":profile,
            "timestamp":timestamp(),"hostname":platform.node(),"run_id":os.environ.get("FEM3D_RUN_ID"),
            "geometry":geometry,"target_size_m":size,"units":"SI","mesh_coordinate_degree":2,
            "tetrahedron_count":int(domain.topology.index_map(3).size_global),
            "geometry_node_count":len(points),"boundary_tags":config["boundary_tags"],"cell_tags":config["cell_tags"],
            "tag_contract":config["tag_contract"],"areas_m2":areas,"area_relative_errors":relative_area,
            "volume_m3":measured_volume,"volume_relative_error":relative_volume,
            "axial_bounds_m":[float(axial.min()),float(axial.max())],"maximum_nodal_radius_m":float(radial.max()),
            "normal_dot_axis":normal_checks,"minimum_jacobian_determinant_m3":float(determinant.min()),
            "wall_time_s":time.perf_counter()-start,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "gmsh_version":gmsh.__version__,"dolfinx_version":dolfinx.__version__,"mpi_ranks":1,"threads":1,"gpu_used":False,
            "config_sha256":config_hash,"msh_sha256":sha256(mshpath),"xdmf_sha256":sha256(base/"pipe.xdmf"),"hdf5_sha256":sha256(base/"pipe.h5"),
            "pde_solved":False}
        write_json(base/"metadata.json",metadata)
        print(json.dumps(metadata,indent=2))
        return metadata
    finally:
        gmsh.finalize()
