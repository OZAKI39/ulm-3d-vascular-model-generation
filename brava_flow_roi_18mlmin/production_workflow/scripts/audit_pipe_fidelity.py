"""Finite-facet cylinder and cap trace QA; independent of any CFD solution."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pyvista as pv
from pipe_surface_qa import surface_qa

ROOT = Path(__file__).resolve().parents[1]


def summary(values):
    return dict(zip(['min', 'median', 'P95', 'P99', 'max'], map(float, np.quantile(values, [0, .5, .95, .99, 1]))))


def boundary_profile(poly, radius, mean=.1):
    """Exact degree-four quadrature of the squared P1-trace error on cap triangles."""
    xyz = np.asarray(poly.points, float)
    tri = poly.faces.reshape(-1, 4)[:, 1:]
    vertices = xyz[tri]
    areas = np.linalg.norm(np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]), axis=1) / 2
    # Perimeter IDs from topology, never a radial tolerance that could erase valid nodes.
    edges, counts = np.unique(np.sort(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1), axis=0, return_counts=True)
    rim = np.unique(edges[counts == 1])
    values = 2 * mean * (1 - (xyz[:, 0]**2 + xyz[:, 1]**2) / radius**2)
    rim_raw = float(np.max(abs(values[rim])))
    # The analytical rim value is zero; native Float32 coordinates can differ by one ULP.
    values[rim] = 0.
    bary, weights = [], []
    for a, b, w in [(.445948490915965, .108103018168070, .223381589678011),
                    (.091576213509771, .816847572980459, .109951743655322)]:
        bary.extend([[a, a, b], [a, b, a], [b, a, a]]); weights.extend([w] * 3)
    bary, weights = np.array(bary), np.array(weights)
    query = np.einsum('qi,tij->tqj', bary, vertices)
    exact = 2 * mean * (1 - np.sum(query[..., :2]**2, axis=-1) / radius**2)
    actual = np.einsum('qi,ti->tq', bary, values[tri])
    integrate = lambda x: float(np.einsum('t,q,tq->', areas, weights, x))
    q = float(np.dot(areas, values[tri].mean(axis=1)))
    result = dict(profile_relative_L2=float(np.sqrt(integrate((actual-exact)**2) / integrate(exact**2))),
                  Q_discrete_m3_s=q, Q_exact_m3_s=np.pi * radius**2 * mean,
                  Q_relative_error=abs(q / (np.pi * radius**2 * mean) - 1),
                  maximum_nodal_velocity_m_s=float(values.max()), rim_nodes=len(rim),
                  rim_velocity_m_s=0., rim_analytic_coordinate_roundoff_m_s=rim_raw,
                  profile_scaling_factor=1., node_count=poly.n_points, triangles=len(tri))
    result['PASS'] = bool(result['profile_relative_L2'] <= .01 and result['Q_relative_error'] <= .005 and rim_raw < 1e-7)
    return result, values


def main(name, directory, input_only=False):
    source = ROOT / 'straight_pipe/manufactured' / name
    config = json.loads((source / 'INPUT.json').read_text())
    radius, length = config['radius_m'], config['length_m']
    path = source / 'model.vtp' if input_only else directory / 'audited_exterior.vtp'
    mesh = pv.read(path); tags = mesh['ModelFaceID']
    wall = mesh.extract_cells(tags == 1).extract_surface(algorithm='dataset_surface')
    central = wall.clip(normal=[0, 0, 1], origin=[0, 0, .4 * length], invert=False)
    central = central.clip(normal=[0, 0, 1], origin=[0, 0, .6 * length], invert=True).triangulate()
    tri = np.asarray(central.points, float)[central.faces.reshape(-1, 4)[:, 1:]]
    xy = tri[..., :2]
    minima = np.full(len(xy), np.inf)
    for i, j in [(0, 1), (1, 2), (2, 0)]:
        a, delta = xy[:, i], xy[:, j] - xy[:, i]
        denom = np.sum(delta**2, axis=1)
        t = np.divide(-np.sum(a * delta, axis=1), denom, out=np.zeros(len(a)), where=denom > 0)
        closest = a + np.clip(t, 0, 1)[:, None] * delta
        minima = np.minimum(minima, np.linalg.norm(closest, axis=1))
    cross = lambda a, b: a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    signs = np.stack([cross(xy[:, j]-xy[:, i], -xy[:, i]) for i,j in [(0,1),(1,2),(2,0)]], axis=1)
    area2 = cross(xy[:,1]-xy[:,0],xy[:,2]-xy[:,0])
    origin_inside = ((signs > 0).all(axis=1) | (signs < 0).all(axis=1)) & (abs(area2) > 0)
    minima[origin_inside] = 0.
    inward_bound = radius - float(minima.min())
    outward_bound = float(np.linalg.norm(xy, axis=2).max()) - radius
    azimuth = 2 * np.pi * (np.arange(2048) + .193) / 2048
    directions = np.c_[np.cos(azimuth), np.sin(azimuth)]
    samples = []
    planes = np.linspace(.4, .6, 21)
    for z in planes:
        cut = wall.slice(normal=[0,0,1], origin=[0,0,z * length])
        lines = cut.lines.reshape(-1,3); assert np.all(lines[:,0] == 2)
        a = np.asarray(cut.points[lines[:,1],:2],float)
        delta = np.asarray(cut.points[lines[:,2],:2],float) - a
        distance = np.empty(len(directions))
        for start in range(0,len(directions),128):
            d = directions[start:start+128]
            denominator = d[:,None,0]*delta[None,:,1]-d[:,None,1]*delta[None,:,0]
            with np.errstate(divide='ignore',invalid='ignore'):
                rr = cross(a,delta)[None,:]/denominator
                f = (a[None,:,0]*d[:,None,1]-a[None,:,1]*d[:,None,0])/denominator
            good=(rr>0)&(f>=-1e-10)&(f<=1+1e-10)&np.isfinite(rr)
            distance[start:start+len(d)] = np.where(good,rr,np.inf).min(axis=1)
        assert np.isfinite(distance).all()
        samples.append(radius-distance)
    qa = surface_qa(path,radius,length)
    p = mesh.points[mesh.faces.reshape(-1,4)[:,1:]]
    areas = np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1)/2
    caps={}
    for tag,label in [(2,'INLET'),(3,'OUTLET')]:
        cap=mesh.extract_cells(tags==tag).extract_surface(algorithm='dataset_surface')
        caps[label],_=boundary_profile(cap,radius)
    geometry = dict(central_z_over_L=[.4,.6], radial_samples=len(planes)*len(azimuth),
                    inward_deviation_m=summary(np.asarray(samples).ravel()),
                    finite_triangle_max_inward_bound_m=inward_bound,
                    finite_triangle_max_outward_bound_m=outward_bound,
                    maximum_absolute_wall_deviation_m=max(inward_bound,outward_bound),
                    finite_triangle_method='Central clipped finite triangles; exact minimum radial norm over projected triangle edges/interior, maximum at vertices. Dense rays provide quantiles, not the acceptance bound.',
                    cap_area_m2={k:v['cap_area_m2'] for k,v in qa['ports'].items()},
                    cap_area_relative_error={k:v['cap_area_relative_error'] for k,v in qa['ports'].items()},
                    wall_area_m2=float(areas[tags==1].sum()),
                    wall_area_relative_error=abs(float(areas[tags==1].sum())/(2*np.pi*radius*length)-1),
                    volume_m3=qa['signed_volume_m3'],
                    volume_relative_error=abs(qa['signed_volume_m3']/(np.pi*radius**2*length)-1))
    geometry['PASS'] = bool(qa['status']=='PASS' and geometry['maximum_absolute_wall_deviation_m'] <= .05e-6 and max(geometry['cap_area_relative_error'].values()) <= .0005)
    result=dict(status='PASS' if geometry['PASS'] and all(c['PASS'] for c in caps.values()) else 'FAIL',
                source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), input_only=input_only,
                surface=geometry,boundary=caps,wall_stress_probe_m=1e-7,wall_stress_is_production_hard_gate=False,
                thresholds=dict(surface_max_m=.05e-6,cap_area_relative=.0005,profile_L2=.01,Q_relative=.005),
                note='Cap velocity trace is evaluated on both ports. Outlet uses separate exact Cauchy traction; its solution velocity is not prescribed.')
    output = source if input_only else directory
    (output/'FIDELITY_QA.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(output/'FIDELITY_SAMPLES.npz', z_over_L=planes,azimuth_rad=azimuth,inward_deviation_m=samples)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('case');parser.add_argument('directory',type=Path);parser.add_argument('--input-only',action='store_true')
    args=parser.parse_args();main(args.case,args.directory,args.input_only)
