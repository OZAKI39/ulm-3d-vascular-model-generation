"""Geometry-only feasibility certificate for the requested balanced design.

This is a necessary-condition diagnostic, not a replacement forward solver.
All fitted/evaluated states still come from the full A network. No CFD input.
"""
import numpy as np
from .hydraulic_resistance import linear_radius_resistance


def balance_metrics(fractions):
    f = np.asarray(fractions, dtype=float)
    if f.shape != (3,) or not np.isfinite(f).all():
        raise ValueError('Three finite signed outlet fractions required')
    return dict(mean_fraction=float(f.mean()), target_mean_fraction=1/3,
                max_abs_deviation=float(np.max(abs(f-1/3))),
                flow_fraction_range=float(np.ptp(f)),
                flow_fraction_std=float(np.std(f)))


def require_balanced_fit(summary, tolerance=1e-5):
    """Fail closed BEFORE freeze/handoff/case creation/solver launch."""
    if summary.get('status') != 'FIT_CONVERGED':
        raise ValueError('BALANCED_0D_TARGET_NOT_REACHED: fit not accepted')
    audit = summary['identifiability_audit']
    if (audit['free_parameters'] != ['s_O1', 's_O2']
            or audit['final_data_jacobian']['rank'] != 2):
        raise ValueError('BALANCED_0D_TARGET_NOT_REACHED: parameter/rank gate')
    p = summary['final_parameters']
    if (p['mu_pa_s'] != .00345312 or p['distal_reference_pa'] != 0
            or p['o3_mode'] != 'legacy_reference'
            or p['terminal_resistance_O3_pa_s_m3'] is not None):
        raise ValueError('The prescribed balanced model must remain unchanged')
    prediction = summary['prediction']
    if (prediction['mass_audit']['status'] != 'PASS'
            or any(prediction['reverse_flow_audit']['port_backflow'].values())):
        raise ValueError('BALANCED_0D_TARGET_NOT_REACHED: mass/backflow gate')
    metrics = balance_metrics([prediction['outlet_flow_fraction'][p] for p in ('O1','O2','O3')])
    if metrics['max_abs_deviation'] >= tolerance:
        raise ValueError('BALANCED_0D_TARGET_NOT_REACHED: balance tolerance')
    return metrics


def passive_o1_o3_certificate(cache, *, mu_pa_s, target_roi_flow_m3_s):
    """Prove a necessary O1/O3 relation from fixed paths in the full graph.

    Requires the two outlets to be unbranched children of the same junction,
    an unmodified ROI tree, and O1's downstream component to touch only O1.
    With a passive positive-R tree to common Pd and Q1>0, P1>=Pd.
    O3 legacy P3=Pd gives Q1/Q3 <= R(J2,O3)/R(J2,O1).
    """
    if not np.isfinite(target_roi_flow_m3_s) or target_roi_flow_m3_s <= 0:
        raise ValueError('Positive SI target flow required')
    resistance = linear_radius_resistance(cache.edge_length_m, cache.baseline_r0_m,
                                         cache.baseline_r1_m, mu=mu_pa_s)
    adjacency = {}
    for edge in np.flatnonzero(cache.internal_roi_edge_mask):
        u, v = map(int, cache.edge_nodes[edge])
        adjacency.setdefault(u, []).append((v, int(edge)))
        adjacency.setdefault(v, []).append((u, int(edge)))
    root = int(cache.port_indices[0])
    parent, parent_edge, stack = {root: None}, {}, [root]
    while stack:
        u = stack.pop()
        for v, edge in adjacency[u]:
            if v not in parent:
                parent[v], parent_edge[v] = u, edge
                stack.append(v)
    leaves = {u for u, a in adjacency.items() if len(a) == 1}
    if (len(parent) != len(adjacency)
            or cache.internal_roi_edge_mask.sum() != len(parent)-1
            or leaves != set(map(int, cache.port_indices))):
        raise ValueError('Certificate requires a connected four-port ROI tree')

    def path(node):
        nodes = [int(node)]
        while parent[nodes[-1]] is not None:
            nodes.append(parent[nodes[-1]])
        return nodes[::-1]

    p1, p3 = path(cache.port_indices[1]), path(cache.port_indices[3])
    common = [a for a, b in zip(p1, p3) if a == b]
    junction = common[-1]
    paths = {}
    for name, nodes in [('O1', p1), ('O3', p3)]:
        nodes = nodes[nodes.index(junction):]
        if any(len(adjacency[n]) != 2 for n in nodes[1:-1]):
            raise ValueError('Additional branch invalidates the series-path proof')
        edges = [parent_edge[n] for n in nodes[1:]]
        paths[name] = dict(node_indices=nodes, original_node_ids=cache.node_ids[nodes],
                           edge_indices=edges, resistance_pa_s_m3=float(resistance[edges].sum()))
    if len(adjacency[junction]) != 3:
        raise ValueError('Expected one parent and two outlet arms')
    mask = np.asarray(cache.downstream_node_masks['O1'])
    crossing = cache.edge_nodes[np.logical_xor(mask[cache.edge_nodes[:,0]], mask[cache.edge_nodes[:,1]])]
    outside = [int(v if mask[u] else u) for u,v in crossing]
    if not outside or set(outside) != {int(cache.port_indices[1])}:
        raise ValueError('O1 downstream component has an additional connection')
    if not mask[cache.terminal_indices].any():
        raise ValueError('No passive distal reference terminal')
    r1, r3 = [paths[p]['resistance_pa_s_m3'] for p in ('O1','O3')]
    k = r3/r1
    # Closest equal-sigma LS point in the larger necessary feasible half-plane.
    # Its attainability for finite radius scales is NOT asserted.
    t = (1+k)/(2*(k*k+k+1))
    limiting = np.array([k*t, 1-(k+1)*t, t]) if k < 1 else np.full(3,1/3)
    return dict(status='EQUAL_SPLIT_IMPOSSIBLE_UNDER_PRESCRIBED_MODEL' if k < 1 else 'NOT_EXCLUDED_BY_THIS_CERTIFICATE',
                scope='necessary condition from fixed full-A geometry; no substitute ROI forward model',
                required_o3_mode='legacy_reference', required_common_distal_reference=True,
                source_geometry_sha256=cache.geometry_sha256,
                node_count=len(cache.node_ids), edge_count=len(cache.edge_nodes),
                roi_node_count=len(adjacency), roi_edge_count=int(cache.internal_roi_edge_mask.sum()),
                junction_node_index=junction, junction_original_id=int(cache.node_ids[junction]),
                junction_xyz_m=cache.xyz_m[junction], paths=paths,
                ratio_bound_Q_O1_over_Q_O3=k,
                identity='P_O1-P_O3 = R_J2_O3*Q_O3 - R_J2_O1*Q_O1',
                passive_constraint='P_O1 >= Pd = P_O3; hence f_O1 <= ratio_bound*f_O3',
                required_P_O1_minus_Pd_for_equal_pa=float((r3-r1)*target_roi_flow_m3_s/3),
                necessary_max_abs_deviation_lower_bound=float(max(0.,(1-k)/(3*(1+k)))),
                equal_sigma_least_squares_boundary_projection=limiting,
                projection_is_not_a_finite_valid_parameter_fit=True)
