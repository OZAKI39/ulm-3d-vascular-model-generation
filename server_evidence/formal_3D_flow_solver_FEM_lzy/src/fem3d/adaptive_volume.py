"""Unchanged Stage 1 medium discrete-boundary volume meshing path."""
import os,platform,resource,signal,time
import numpy as np
from .audit import write_json,sha256,timestamp
from .mesh_input import FACET_NAMES


def generate_volume(surface,out,policy):
    import gmsh
    out.mkdir(parents=True,exist_ok=True);meshdir=out/'mesh';meshdir.mkdir(exist_ok=True)
    if (meshdir/'fluid.msh').exists():raise RuntimeError('Volume evidence already exists')
    options=policy['volume_meshing']['effective_gmsh_options'];start=time.perf_counter()
    limits=policy['limits'];bound=int(limits['volume_virtual_memory_gib']*1024**3)
    resource.setrlimit(resource.RLIMIT_AS,(bound,bound))
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('Volume meshing time limit')))
    signal.alarm(limits['volume_wall_time_s'])
    points,tri,tags=surface['points_m'],surface['triangles'],surface['facet_tags']
    meta={'timestamp':timestamp(),'run_id':os.environ.get('FEM3D_RUN_ID'),'hostname':platform.node(),'gpu_used':False,'effective_gmsh_options':options,'gmsh_version':gmsh.__version__}
    try:
        gmsh.initialize();gmsh.model.add('adaptive_port_volume')
        for k,v in options.items():gmsh.option.setNumber(k,v)
        for tag in FACET_NAMES:gmsh.model.addDiscreteEntity(2,tag)
        gmsh.model.mesh.addNodes(2,1,np.arange(1,len(points)+1),points.ravel())
        for tag in FACET_NAMES:
            idx=np.flatnonzero(tags==tag);gmsh.model.mesh.addElementsByType(tag,2,idx+1,(tri[idx]+1).ravel())
        gmsh.model.mesh.reclassifyNodes();loop=gmsh.model.geo.addSurfaceLoop(list(FACET_NAMES));volume=gmsh.model.geo.addVolume([loop]);gmsh.model.geo.synchronize()
        for tag,name in FACET_NAMES.items():gmsh.model.addPhysicalGroup(2,[tag],tag,name)
        gmsh.model.addPhysicalGroup(3,[volume],100,'FLUID');gmsh.model.mesh.generate(3)
        types,_,_=gmsh.model.mesh.getElements(3)
        if list(types)!=[4]:raise ValueError('Non-tetra volume')
        eid,conn=gmsh.model.mesh.getElementsByType(4);nodes,coords,_=gmsh.model.mesh.getNodes();order=np.argsort(nodes);nodes=nodes[order];points=coords.reshape(-1,3)[order]
        tetra=np.searchsorted(nodes,conn).reshape(-1,4);bt=[];values=[]
        for tag in FACET_NAMES:
            _,conn=gmsh.model.mesh.getElementsByType(2,tag);bt.append(np.searchsorted(nodes,conn).reshape(-1,3));values.extend([tag]*(len(conn)//3))
        data={'points_m':points,'tetra':tetra,'boundary_triangles':np.vstack(bt),'facet_tags':np.array(values,dtype=np.int32),
              'cell_tags':np.full(len(tetra),100,np.int32),'gmsh_element_ids':eid,'min_sicn':gmsh.model.mesh.getElementQualities(eid,'minSICN')}
        gmsh.write(str(meshdir/'fluid.msh'));np.savez_compressed(meshdir/'volume_mesh.npz',**data)
        meta.update(status='PASS',mesh_sha256=sha256(meshdir/'volume_mesh.npz'),actual_tetra=len(tetra),actual_vertices=len(points))
        return data
    except BaseException as exc:meta.update(status='FAIL',error=repr(exc));raise
    finally:
        signal.alarm(0)
        if gmsh.isInitialized():gmsh.finalize()
        meta.update(wall_time_s=time.perf_counter()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write_json(out/'metadata/meshing.json',meta)


def dolfinx_loadability(path):
    from mpi4py import MPI
    import dolfinx
    from dolfinx.io import gmsh as gmsh_io
    from dolfinx import mesh
    d=gmsh_io.read_from_msh(str(path),MPI.COMM_SELF,rank=0,gdim=3);m=d.mesh
    assert m.topology.dim==m.geometry.dim==3 and m.topology.cell_type==mesh.CellType.tetrahedron
    nc=m.topology.index_map(3).size_local;m.topology.create_connectivity(2,3);nf=m.topology.index_map(2).size_local
    assert len(d.cell_tags.indices)==nc and np.all(d.cell_tags.values==100)
    ext=mesh.exterior_facet_indices(m.topology);tagged=d.facet_tags.indices[d.facet_tags.indices<nf]
    assert np.array_equal(np.sort(ext),np.sort(tagged)) and set(d.facet_tags.values)==set(FACET_NAMES)
    return {'status':'PASS','dolfinx_version':dolfinx.__version__,'tdim':3,'gdim':3,'cell_type':'tetrahedron','owned_cells':nc,'tagged_exterior':len(tagged),'fem_space_created':False}
