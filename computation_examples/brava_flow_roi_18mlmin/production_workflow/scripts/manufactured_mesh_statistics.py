"""Actual final tetra sizes for one manufactured case; no mesh mutation."""
from pathlib import Path
import argparse
import itertools
import json
import numpy as np
import pyvista as pv
from pipe_probe_geometry import strict_cells

ROOT=Path(__file__).resolve().parents[1]


def stats(x):
    return dict(zip(['min','P5','median','P95','max'],map(float,np.quantile(x,[0,.05,.5,.95,1]))))


def main(name):
    mesh=ROOT/'manufactured_pipe'/name/'mesh'
    config=json.loads((ROOT/'straight_pipe/manufactured'/name/'INPUT.json').read_text())
    with np.load(mesh/'BOUNDARY_TOPOLOGY.npz') as data:
        points,tetra=data['points_m'],data['tetra']
        mask=data['facet_tags']==1
        wall,owner=data['boundary_triangles'][mask],data['adjacent_tetra'][mask]
    xyz=points[wall];centers=xyz.mean(axis=1)
    n=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
    height=np.abs(np.einsum('nvi,ni->nv',points[tetra[owner]]-centers[:,None,:],n)).max(axis=1)
    length=config['length_m'];radius=config['radius_m']
    central=(centers[:,2]>=.4*length)&(centers[:,2]<=.6*length)
    probes=np.c_[np.zeros((21,2)),np.linspace(.4,.6,21)*length]
    grid=pv.read(mesh/'volume.vtu');cells,_=strict_cells(grid,points,tetra,probes);assert np.all(cells>=0)
    edge=np.array(list(itertools.combinations(range(4),2)))
    h=np.linalg.norm(points[tetra[cells][:,edge[:,0]]]-points[tetra[cells][:,edge[:,1]]],axis=2).mean(axis=1)
    result=dict(nodes=len(points),tetra=len(tetra),final_cell_types={'10':len(tetra)},h_bulk_over_R_requested=config['policy']['h_bulk_over_R'],
                final_wall_adjacent_tetra_normal_height_m=stats(height[central]),
                final_wall_adjacent_tetra_normal_height_over_R=stats(height[central]/radius),
                central_axis_D_over_mean_tet_edge=stats(2*radius/h),
                actual_layer_count=None,actual_prism_cells=0,
                interpretation='Measured final tetra normal altitudes, not a prism-layer claim. Axis D/h is a proxy, not literal counted cells across the lumen.')
    (mesh/'MESH_STATISTICS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');main(p.parse_args().case)
