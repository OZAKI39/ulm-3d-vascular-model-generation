"""Exact nodal coefficient checkpoint plus DOLFINx XDMF geometry.

Coordinates identify the CG nodal DOFs across partitionings. This restores the
original P2/P1 coefficients; interpolation to visualization spaces is separate.
"""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .audit import sha256
from .spaces import create_spaces,real_owned_value


def load_primary(root,case,comm):
    from dolfinx import fem,io,mesh
    root=Path(root)
    base=root/"outputs/stage02/cases"/case
    metadata=json.loads((base/"metadata/run.json").read_text())
    meshbase=root/"outputs/stage02/meshes"/metadata["profile"]
    assert sha256(meshbase/"pipe.h5")==metadata["mesh_sha256"]
    assert sha256(base/"solution/primary_checkpoint.npz")==metadata["checkpoint_sha256"]
    with io.XDMFFile(comm,str(meshbase/"pipe.xdmf"),"r") as file:
        domain=file.read_mesh(name="pipe_mesh",ghost_mode=mesh.GhostMode.shared_facet)
        domain.topology.create_connectivity(2,3)
        facet_tags=file.read_meshtags(domain,name="facet_tags")
    V,P,R,_=create_spaces(domain)
    saved=np.load(base/"solution/primary_checkpoint.npz")
    fields={}
    displacements={}
    for name,space in (("velocity",V),("pressure",P)):
        coords=space.tabulate_dof_coordinates()
        count=space.dofmap.index_map.size_local
        blocksize=space.dofmap.index_map_bs
        source=saved[name+"_coordinates_m"]
        assert len(np.unique(source,axis=0))==len(source)
        distance,index=cKDTree(source).query(coords[:count])
        tolerance=512*np.finfo(float).eps*max(np.abs(source).max(),1e-12)
        assert distance.max(initial=0)<=tolerance and len(np.unique(index))==len(index)
        assert len(source)==space.dofmap.index_map.size_global
        field=fem.Function(space,name="velocity" if name=="velocity" else "pressure_gauge_pa")
        field.x.array[:count*blocksize]=saved[name+"_values"][index].ravel()
        field.x.scatter_forward()
        assert np.array_equal(field.x.array[:count*blocksize], saved[name+"_values"][index].ravel())
        fields[name]=field
        displacements[name]=float(distance.max(initial=0))
    multiplier=fem.Function(R,name="inlet_normal_traction_multiplier")
    multiplier.x.array[:R.dofmap.index_map.size_local]=float(saved["lambda_pa"])
    multiplier.x.scatter_forward()
    fields.update(multiplier=multiplier,lambda_pa=real_owned_value(multiplier))
    return domain,facet_tags,fields,metadata,displacements
