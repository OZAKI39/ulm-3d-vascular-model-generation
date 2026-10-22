"""Native SimVascular/TetGen interface and unchanged tetrahedron QA algorithms."""
from pathlib import Path
import subprocess
import numpy as np
from .paths import ROOT, read_json, sha256, fresh_directory, write_json
from vascular_validation.mesh_diagnostics import volume_validity, tetra_sicn, quality_summary

def check_existing(root=ROOT):
    with np.load(root/'inputs/mesh/SV_MESH/mesh_arrays.npz') as data:
        points,tetra=data['points_m'],data['tetra']
        validity=volume_validity(points,tetra)
        quality=quality_summary(tetra_sicn(points,tetra))
        assert quality['q_min']>0
        assert set(np.unique(data['facet_tags']))=={1,2,3,4,5}
        return dict(status='PASS',points=len(points),tetra=len(tetra),validity=validity,quality=quality)

def prepare(output, root=ROOT):
    output=fresh_directory(output,root)
    native=read_json(root/'configs/cfd/native.json')
    executable=Path(native['mesher'])
    if sha256(executable)!=native['mesher_sha256']:
        raise ValueError('Mesher launcher identity changed')
    if sha256(native['mesher_binary'])!=native['mesher_binary_sha256']:
        raise ValueError('SimVascular executable identity changed')
    command=[str(executable),'-python','--',str(root/'native/tetgen_mesh.py'),
             str(root/'inputs/geometry/labelled_core_m.vtp'),str(output/'generated')]
    write_json(output/'PLAN.json',dict(command=command,executed=False,
        source_surface_sha256=sha256(root/'inputs/geometry/labelled_core_m.vtp'),global_edge_size_m=.0002,
        generator='SimVascular 2023-05-31 native TetGen',boundary_layers=False,automatic_CFD=False))
    return command

def run(command, output):
    with (Path(output)/'execution.log').open('x') as stream:
        subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True)

