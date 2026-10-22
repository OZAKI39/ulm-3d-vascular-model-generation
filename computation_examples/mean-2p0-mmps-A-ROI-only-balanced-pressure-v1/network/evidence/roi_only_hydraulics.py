"""Fixed ROI geometry, inlet flow and three outlet pressures; no exterior solve."""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve
from .audit import sha256
from .hydraulic_resistance import linear_radius_resistance, MU_PA_S, ROI_TARGET_Q_M3_S

MODEL_NAME = 'ROI_ONLY_PRESSURE_BOUNDARY_0D_V1'
PORTS = ('INLET', 'O1', 'O2', 'O3')
OUTLETS = PORTS[1:]


def readonly(values):
    a = np.ascontiguousarray(values)
    return np.frombuffer(a.tobytes(), dtype=a.dtype).reshape(a.shape)


@dataclass(frozen=True)
class ROIHydraulicCache:
    node_ids: np.ndarray
    xyz_m: np.ndarray
    radius_m: np.ndarray
    source_node_indices: np.ndarray
    roi_edge_indices: np.ndarray
    edge_nodes: np.ndarray
    edge_length_m: np.ndarray
    edge_resistance_pa_s_m3: np.ndarray
    port_node_indices: np.ndarray
    port_edge_indices: np.ndarray
    port_coefficients: np.ndarray
    geometry_sha256: str
    geometry_file_sha256: str
    port_mapping_sha256: str
    mu_pa_s: float

    def summary(self):
        return dict(model_name=MODEL_NAME, node_count=len(self.node_ids), edge_count=len(self.edge_nodes),
                    connected=True, full_A_network_used_for_boundary_generation=False,
                    geometry_sha256=self.geometry_sha256, geometry_file_sha256=self.geometry_file_sha256,
                    port_mapping_sha256=self.port_mapping_sha256, dynamic_viscosity_pa_s=self.mu_pa_s,
                    ports=[dict(port=p, node_id=int(self.node_ids[n]), xyz_m=self.xyz_m[n].tolist(),
                                local_node_index=int(n), local_edge_index=int(e),
                                source_edge_index=int(self.roi_edge_indices[e]), coefficient=int(c),
                                positive='INTO_ROI' if p == 'INLET' else 'OUT_OF_ROI')
                           for p,n,e,c in zip(PORTS,self.port_node_indices,self.port_edge_indices,self.port_coefficients)])


def load_roi_cache(geometry_npz, ports_json, *, mu_pa_s=MU_PA_S):
    mapping = json.loads(Path(ports_json).read_text())
    if mapping.get('status') != 'PASS' or mapping.get('extra_saved_roi_ports') != 0 or mapping.get('expected_fem_ports') != 4:
        raise ValueError('Verified four-port mapping required')
    ports = {p['name']:p for p in mapping['ports']}
    if len(mapping['ports']) != 4 or set(ports) != set(PORTS):
        raise ValueError('Exactly INLET/O1/O2/O3 required')
    with np.load(geometry_npz, allow_pickle=False) as g:
        mask = g['roi_internal_edge_mask']
        if mask.dtype != bool or mask.shape != (len(g['edges']),):
            raise ValueError('Expected boolean internal ROI edge mask')
        ei = np.flatnonzero(mask)
        retained = g['edges'][ei]
        ni = np.unique(retained)
        lookup = np.full(len(g['ids']), -1, dtype=int); lookup[ni] = np.arange(len(ni))
        ids, xyz, radii = g['ids'][ni], g['xyz_m'][ni], g['radius_m'][ni]
        edges = lookup[retained]
    if len(np.unique(ids)) != len(ids) or not np.isfinite(xyz).all() or not np.isfinite(radii).all() or np.any(radii <= 0):
        raise ValueError('Invalid ROI geometry')
    u,v = edges.T
    graph = coo_matrix((np.ones(2*len(u)), (np.r_[u,v],np.r_[v,u])), shape=(len(ids),len(ids))).tocsr()
    if connected_components(graph,directed=False)[0] != 1 or len(edges) != len(ids)-1:
        raise ValueError('Expected connected ROI tree')
    pn,pe,pc = [],[],[]
    for name in PORTS:
        port = ports[name]; matches = np.flatnonzero(ids == port['network_node_id'])
        if len(matches) != 1: raise ValueError('Missing/duplicate port node')
        n = int(matches[0]); incident = np.flatnonzero(np.any(edges == n,axis=1))
        if len(incident) != 1: raise ValueError('Port needs exactly one ROI edge')
        e = int(incident[0]); a,b = edges[e]; other = b if a == n else a
        normal = np.asarray(port['outward_normal']); inward = xyz[other]-xyz[n]
        if not np.isclose(np.linalg.norm(normal),1.,rtol=1e-9) or np.dot(inward/np.linalg.norm(inward),normal) > -.99:
            raise ValueError('Port outward normal disagrees with retained ROI geometry')
        # Algebraic storage u->v is not the blood-flow direction.
        into = 1 if a == n else -1
        pn.append(n); pe.append(e); pc.append(into if name == 'INLET' else -into)
    if set(np.flatnonzero(np.diff(graph.indptr) == 1)) != set(pn):
        raise ValueError('Unmapped ROI boundary leaf')
    length = np.linalg.norm(xyz[v]-xyz[u],axis=1)
    resistance = linear_radius_resistance(length,radii[u],radii[v],mu=mu_pa_s)
    if np.any(length <= 0) or not np.isfinite(resistance).all() or np.any(resistance <= 0):
        raise ValueError('Invalid positive resistance')
    arrays = [ids,xyz,radii,ni,ei,edges,length,resistance,np.array(pn),np.array(pe),np.array(pc)]
    digest = hashlib.sha256()
    for a in arrays:
        digest.update(str((a.shape,a.dtype.str)).encode()); digest.update(np.ascontiguousarray(a).tobytes())
    return ROIHydraulicCache(*map(readonly,arrays),digest.hexdigest(),sha256(geometry_npz),sha256(ports_json),float(mu_pa_s))


def solve_roi_mixed_bc(cache, outlet_pressure_pa, Qin_target_m3_s=ROI_TARGET_Q_M3_S):
    pb = np.asarray(outlet_pressure_pa,dtype=float)
    if pb.shape != (3,) or not np.isfinite(pb).all() or not np.isfinite(Qin_target_m3_s) or Qin_target_m3_s <= 0:
        raise ValueError('Three finite pressures and positive inlet flow required')
    n = len(cache.node_ids); u,v = cache.edge_nodes.T; g = 1/cache.edge_resistance_pa_s_m3
    G = coo_matrix((np.r_[g,g,-g,-g],(np.r_[u,v,u,v],np.r_[u,v,v,u])),shape=(n,n)).tocsr()
    fixed = cache.port_node_indices[1:]; free = np.setdiff1d(np.arange(n),fixed)
    rhs = np.zeros(n); rhs[cache.port_node_indices[0]] = Qin_target_m3_s
    p = np.zeros(n); p[fixed] = pb
    # Conditioning rescale is numerical only; no operating-point flow scaling.
    scale = float(g.max())
    p[free] = spsolve(G[free][:,free]/scale,(rhs[free]-G[free][:,fixed]@pb)/scale)
    q = (p[u]-p[v])*g
    ports = cache.port_coefficients*q[cache.port_edge_indices]
    residual = G@p-rhs
    mass = dict(max_free_node_relative_residual=float(np.max(abs(residual[free]))/Qin_target_m3_s),
                inlet_relative_error=float(abs(ports[0]-Qin_target_m3_s)/Qin_target_m3_s),
                ROI_relative_residual=float(abs(ports[0]-sum(ports[1:]))/Qin_target_m3_s))
    mass['status'] = 'PASS' if max(mass.values()) < 1e-9 else 'FAIL'
    if not np.isfinite(p).all() or mass['status'] != 'PASS': raise ValueError('Mixed-BC solve/mass failure')
    return dict(pressure_pa=p, edge_flow_m3_s=q, port_pressure_pa=p[cache.port_node_indices],
                port_flow_m3_s=ports,outlet_fraction=ports[1:]/Qin_target_m3_s,mass_audit=mass)
