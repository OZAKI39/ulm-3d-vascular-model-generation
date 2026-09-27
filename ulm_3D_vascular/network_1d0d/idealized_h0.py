"""H0: explicitly assumed structural source and equal-pressure structural leaves.

Reuse the v1 SI graph, exact cut nodes, resistance integral and sparse solver.
No physiological identity is required or inferred by this idealized model.
"""
from dataclasses import dataclass
from pathlib import Path
import csv
import json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .audit import sha256
from .hydraulic_resistance import linear_radius_resistance, MU_PA_S, ROI_TARGET_Q_M3_S
from .network_solver import solve_network, operating_point_scale
from .boundary_conditions import require_source

MODEL_NAME = 'A-Network Idealized Hydraulic Baseline Model v2'
SOURCE_ID = 2410
TERMINAL_MODEL = 'ALL_STRUCTURAL_LEAVES_REFERENCE_PRESSURE'


@dataclass
class H0Domain:
    ids: np.ndarray
    xyz_m: np.ndarray
    radius_m: np.ndarray
    edges: np.ndarray
    raw_edge: np.ndarray
    fractions: np.ndarray
    internal_roi: np.ndarray
    ports: list
    port_nodes: np.ndarray
    signs: list
    source: int
    terminals: np.ndarray
    exterior_labels: np.ndarray

    def resistances(self, radii=None):
        r = self.radius_m if radii is None else np.asarray(radii)
        e = self.edges
        return linear_radius_resistance(np.linalg.norm(self.xyz_m[e[:, 1]]-self.xyz_m[e[:, 0]], axis=1), r[e[:, 0]], r[e[:, 1]], MU_PA_S)


def h0_boundary(ids, edges, source_id=SOURCE_ID):
    ids = np.asarray(ids); edges = np.asarray(edges)
    source = np.flatnonzero(ids == source_id)
    if len(source) != 1: raise ValueError('H0 source must exist exactly once')
    degree = np.bincount(edges.ravel(), minlength=len(ids))
    if degree[source[0]] != 1: raise ValueError('Expected the verified structural root endpoint')
    terminals = np.flatnonzero((degree == 1) & (ids != source_id))
    if not len(terminals): raise ValueError('H0 requires at least one reference-pressure terminal')
    require_source(dict(hydraulic_source_id=int(source_id), identification_kind='EXPLICIT_MODEL_ASSUMPTION',
                        evidence_reference='User-defined H0 model, sections 1-3: structural root 2410 and all other leaves'))
    return int(source[0]), terminals


def port_signs(ids, xyz, edges, internal, ports):
    """Orient by the actual retained ROI side, independently of edge storage order."""
    result = []
    for port in ports:
        node = int(np.flatnonzero(ids == port['network_node_id'])[0])
        incident = np.flatnonzero(np.any(edges == node, axis=1))
        inside_edges = [int(e) for e in incident if internal[e]]
        if len(inside_edges) != 1: raise ValueError('Port must have one retained ROI-side segment')
        eid = inside_edges[0]; u, v = edges[eid]
        neighbor = int(v if u == node else u)
        to_inside = xyz[neighbor]-xyz[node]; to_inside /= np.linalg.norm(to_inside)
        normal = np.asarray(port['outward_normal'])
        dot = float(np.dot(to_inside, normal))
        if dot > -.99: raise ValueError('Saved normal disagrees with actual ROI-side geometry')
        # Stored edge q is u->v, merely an algebraic convention, never assumed blood flow.
        into_roi_coefficient = 1 if u == node else -1
        coefficient = into_roi_coefficient if port['name'] == 'INLET' else -into_roi_coefficient
        result.append(dict(port=port['name'], node_index=node, node_id=int(ids[node]), roi_side_neighbor_id=int(ids[neighbor]),
                           roi_segment=eid, algebraic_edge=[int(ids[u]), int(ids[v])], coefficient=coefficient,
                           positive_definition='INTO_ROI' if port['name'] == 'INLET' else 'OUT_OF_ROI',
                           to_roi_unit_vector=to_inside.tolist(), outward_normal=normal.tolist(), to_roi_dot_outward=dot,
                           flow_direction_source='pressure solution + ROI-side topology, not SWC orientation'))
    return result


def load_h0(v1):
    v1 = Path(v1); provenance = json.loads((v1/'data/a_network_provenance.json').read_text())
    for path, digest in json.loads((v1/'data/protected_input_hashes.json').read_text()).items():
        if sha256(path) != digest: raise ValueError('Protected input changed: '+path)
    mapping = json.loads((v1/'data/roi_ports_in_a.json').read_text())
    if mapping['status'] != 'PASS' or mapping['extra_saved_roi_ports'] or mapping['expected_fem_ports'] != 4:
        raise ValueError('ROI port mapping is not verified')
    with np.load(Path(provenance['lineage']['analysis_swc']).with_suffix('.npz')) as analysis:
        wanted = set(map(int, analysis['node_ids']))
    ports = mapping['ports']; wanted.update(p['network_node_id'] for p in ports)
    with np.load(v1/'data/a_graph_with_exact_cuts_si.npz') as raw:
        keep = np.array([int(i) in wanted for i in raw['ids']]); nodes = np.flatnonzero(keep)
        lookup = np.full(len(keep), -1, dtype=int); lookup[nodes] = np.arange(len(nodes))
        edge_keep = np.all(keep[raw['edges']], axis=1)
        ids, xyz, radius = raw['ids'][keep], raw['xyz_m'][keep], raw['radius_m'][keep]
        edges = lookup[raw['edges'][edge_keep]]; raw_edge = raw['original_edge'][edge_keep]; fractions = raw['fractions'][edge_keep]
    with open(v1/'data/roi_edge_provenance.csv') as f:
        retained = {int(r['raw_edge_id']):(float(r['start_fraction']), float(r['end_fraction'])) for r in csv.DictReader(f)}
    inside = np.array([eid in retained and a >= retained[eid][0]-1e-12 and b <= retained[eid][1]+1e-12 for eid, (a,b) in zip(raw_edge, fractions)])
    source, terminals = h0_boundary(ids, edges)
    u, v = edges.T; adj = coo_matrix((np.ones(2*len(edges)), (np.r_[u,v], np.r_[v,u])), shape=(len(ids),len(ids)))
    components, _ = connected_components(adj, directed=False)
    if components != 1 or len(edges)-len(ids)+1 != 0: raise ValueError('Analysis A topology differs from verified v1')
    u,v = edges[~inside].T
    _, exterior = connected_components(coo_matrix((np.ones(2*len(u)),(np.r_[u,v],np.r_[v,u])),shape=adj.shape), directed=False)
    signs = port_signs(ids, xyz, edges, inside, ports)
    return H0Domain(ids, xyz, radius, edges, raw_edge, fractions, inside, ports, np.array([s['node_index'] for s in signs]), signs, source, terminals, exterior)


def solve_operating_point(domain, radii=None, removed_terminal=None):
    terminals = domain.terminals[domain.terminals != removed_terminal]
    if not len(terminals): raise ValueError('INVALID_VARIANT: no sink')
    bc = {int(i):0. for i in terminals}; bc[domain.source] = 1.
    unit = solve_network(len(domain.ids), domain.edges, domain.resistances(radii), bc)
    q_unit = np.array([unit.flow[x['roi_segment']]*x['coefficient'] for x in domain.signs])
    if q_unit[0] <= 0:
        return dict(status='IDEALIZED_SOURCE_ORIENTATION_CONFLICT', unit=unit, unit_roi_q=q_unit, scaling_lambda=None)
    scale = operating_point_scale(q_unit[0], ROI_TARGET_Q_M3_S)
    p, q, nodal, port_q = unit.pressure*scale, unit.flow*scale, unit.node_outflow*scale, q_unit*scale
    interior = np.setdiff1d(np.arange(len(domain.ids)), np.r_[domain.source, terminals])
    qsource = float(nodal[domain.source]); qterm = float(-nodal[terminals].sum())
    imbalance = nodal[interior]
    mass = dict(max_absolute_nodal_imbalance_m3s=float(np.max(abs(imbalance))), rms_nodal_imbalance_m3s=float(np.sqrt(np.mean(imbalance**2))),
                max_relative_nodal_imbalance_to_A_source=float(np.max(abs(imbalance))/abs(qsource)),
                A_total_source_inflow=qsource, A_total_terminal_outflow=qterm,
                global_absolute_residual_m3s=abs(qsource-qterm), global_relative_residual=abs(qsource-qterm)/abs(qsource),
                ROI_absolute_residual_m3s=float(abs(port_q[1:].sum()-port_q[0])),
                ROI_relative_residual=float(abs(port_q[1:].sum()-port_q[0])/port_q[0]),
                ROI_target_relative_error=float(abs(port_q[0]-ROI_TARGET_Q_M3_S)/ROI_TARGET_Q_M3_S),
                no_numerical_values_forced_to_zero=True)
    mass['status'] = 'PASS' if max(mass[k] for k in ['max_relative_nodal_imbalance_to_A_source','global_relative_residual','ROI_relative_residual']) < 1e-9 else 'FAIL'
    positive = np.maximum(port_q[1:], 0); positive /= positive.sum()
    r = domain.radius_m if radii is None else radii
    rows = []
    for j, port in enumerate(domain.ports):
        rows.append(dict(port=port['name'], node_id=int(domain.ids[domain.port_nodes[j]]), pressure_realcut_Pa=float(p[domain.port_nodes[j]]),
                         signed_Q_m3s=float(port_q[j]), mean_velocity_estimate_m_s=float(port_q[j]/(np.pi*r[domain.port_nodes[j]]**2)),
                         signed_fraction=float(port_q[j]/port_q[0]), positive_outflow_normalized_fraction=float(positive[j-1]) if j else None,
                         backflow=bool(port_q[j] < 0), pressure_definition='MODEL_DEFINED_REFERENCE_TERMINAL' if port['name']=='O3' and removed_terminal!=domain.port_nodes[j] else 'NETWORK_SOLUTION'))
    return dict(status='PASS' if mass['status']=='PASS' else 'A_NETWORK_H0_NUMERICAL_FAIL', unit=unit, unit_roi_q=q_unit,
                scaling_lambda=scale, pressure=p, flow=q, node_outflow=nodal, internal_nodes=interior, terminals=terminals,
                ports=rows, fractions=port_q[1:]/port_q[0], mass=mass, solver_audit=unit.audit)


def external_radius_variant(domain, port_name, factor):
    j = next(i for i,p in enumerate(domain.ports) if p['name']==port_name)
    if port_name == 'O3': raise ValueError('O3 has no saved downstream geometry')
    node = domain.port_nodes[j]
    mask = domain.exterior_labels == domain.exterior_labels[node]
    mask[node] = False  # Keep real-cut radius/ROI geometry fixed; first exterior segment tapers continuously.
    radii = domain.radius_m.copy(); radii[mask] *= factor
    return radii, mask


def downstream_resistance(domain, baseline, name):
    if name == 'O3': return dict(status='NOT_AVAILABLE', closure='DIRICHLET_REFERENCE_TERMINAL', pressure_Pa=0.)
    j = next(i for i,p in enumerate(domain.ports) if p['name']==name); node=domain.port_nodes[j]
    nodes=np.flatnonzero(domain.exterior_labels==domain.exterior_labels[node]); index=np.full(len(domain.ids),-1,int);index[nodes]=np.arange(len(nodes))
    edge_mask=~domain.internal_roi & np.all(index[domain.edges]>=0,axis=1); e=index[domain.edges[edge_mask]]
    refs=[int(index[t]) for t in domain.terminals if index[t]>=0]
    bc={t:0. for t in refs};bc[int(index[node])]=1.
    sol=solve_network(len(nodes),e,domain.resistances()[edge_mask],bc);q=sol.node_outflow[index[node]];R=1/q
    op=baseline['ports'][j]; ratio=op['pressure_realcut_Pa']/op['signed_Q_m3s']
    return dict(status='PASS', R_Pa_s_m3=float(R), full_A_p_over_Q=float(ratio), relative_equivalence_error=float(abs(R-ratio)/R),
                external_node_count=len(nodes), external_edge_count=len(e), reference_terminal_count=len(refs),
                method='independent exterior sparse unit-port-pressure solve')
