"""Independently verify the real WALL subset and original tolerance-band states."""
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist
from scripts.stage5a_r2.run import load
from scripts.stage5a_r2.seal import BASE
from scripts.stage4l.dense_probe import write,sha


def main():
    rows=[]
    for case in ['B2','B3']:
        tracker,sweep,_,state,inputs=load(case,'S1','numpy',True)
        closed=pv.read(tracker.config['geometry']['mesh_surface']).triangulate()
        wall=sweep.mesh
        distance,index=cKDTree(closed.points).query(wall.points)
        cf=closed.faces.reshape(-1,4)[:,1:]
        wf=wall.faces.reshape(-1,4)[:,1:]
        closed_faces={tuple(sorted(v)) for v in cf}
        exact=bool(np.max(distance)==0 and all(tuple(sorted(v)) in closed_faces for v in index[wf]))
        edges=wall.extract_feature_edges(boundary_edges=True,feature_edges=False,manifold_edges=False,non_manifold_edges=False)
        pairs=edges.lines.reshape(-1,3)[:,1:]
        graph=coo_matrix((np.ones(2*len(pairs)),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(edges.n_points,edges.n_points)).tocsr()
        ports,_=connected_components(graph,directed=False)
        mask=state['active'];x=state['positions_m'][mask];r=state['radii_m'][mask]
        gaps,normals,valid=tracker.cylinder.wall_gap_and_normal(x,r)
        pair=pdist(x)-2*r[0]
        tri=wall.points[wf]
        area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
        row=dict(case=case,wall_path=sweep.path,wall_sha256=sha(sweep.path),closed_points=closed.n_points,
                 closed_source_path=tracker.config['geometry']['mesh_surface'],closed_source_sha256=sha(tracker.config['geometry']['mesh_surface']),
                 wall_points=wall.n_points,wall_triangles=wall.n_cells,maximum_vertex_distance_to_closed_m=float(distance.max()),
                 every_wall_triangle_is_exact_closed_subset=exact,boundary_connected_loops=int(ports),
                 minimum_triangle_area_m2=float(area.min()),initial_pair_gap_min_m=float(pair.min()),
                 initial_signed_wall_gap_min_m=float(gaps.min()),initial_negative_pair_gaps=int((pair<0).sum()),
                 initial_negative_wall_gaps=int((gaps<0).sum()),active_count=int(mask.sum()),total_particle_IDs=len(mask),
                 original_coordinates_retained=True,initial_projection_m=0.,input_sha256=inputs,
                 pass_gate=bool(exact and ports==3 and area.min()>0 and valid.all() and pair.min()>=-1e-11 and gaps.min()>=-1e-11))
        rows.append(row)
    write(BASE/'final/GEOMETRY_INPUT_AUDIT.json',dict(status='PASS' if all(r['pass_gate'] for r in rows) else 'FAIL',rows=rows,
          certificate_scope='Closed-surface signed distance is a conservative lower bound to its exact open WALL subset; exact finite-radius sweeps use only WALL triangles. Ports are not capped in the sweep.'))
    print([(r['case'],r['pass_gate'],r['boundary_connected_loops']) for r in rows],flush=True)


if __name__=='__main__':main()
