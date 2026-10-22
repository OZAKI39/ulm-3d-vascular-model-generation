"""Pure ROI pressure design. No CFD, full-A forward solve or radius parameters."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.optimize import least_squares
import yaml
from .roi_only_hydraulics import MODEL_NAME, PORTS, OUTLETS, solve_roi_mixed_bc
from .hydraulic_resistance import ROI_TARGET_Q_M3_S


def metrics(f):
    f = np.asarray(f)
    return dict(max_abs_deviation_from_one_third=float(max(abs(f-1/3))),range=float(np.ptp(f)),
                std=float(np.std(f)),J_balance=float(np.sum((f-1/3)**2)))


def regress_h0(cache, operating_point, Qin=ROI_TARGET_Q_M3_S):
    rows = {r['port']:r for r in operating_point['ports']}
    p = [rows[k]['pressure_realcut_Pa'] for k in OUTLETS]
    result = solve_roi_mixed_bc(cache,p,Qin)
    expected = np.array([rows[k]['signed_Q_m3s'] for k in PORTS])
    error = float(max(abs(result['port_flow_m3_s']-expected))/Qin)
    if error >= 1e-9: raise ValueError('STOP: full-A extracted port state regression failed')
    return dict(status='PASS',maximum_port_flow_error_relative_Qin=error,
                reference_port_flow_m3_s=expected.tolist(),reference_port_pressure_pa=[rows[k]['pressure_realcut_Pa'] for k in PORTS],
                reconstructed={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in result.items()},
                tolerance_relative_Qin=1e-9)


def pressure_basis(cache, Qin=ROI_TARGET_Q_M3_S):
    q0 = solve_roi_mixed_bc(cache,[0,0,0],Qin)['port_flow_m3_s'][1:]
    B = np.column_stack([solve_roi_mixed_bc(cache,p,Qin)['port_flow_m3_s'][1:]-q0 for p in ([1,0,0],[0,1,0])])
    rank = int(np.linalg.matrix_rank(B)); condition = float(np.linalg.cond(B))
    if rank != 2: raise ValueError('STOP: rank-deficient pressure basis')
    return dict(q0_m3_s=q0.tolist(),B_m3_s_per_pa=B.tolist(),rank=rank,condition_number=condition,
                singular_values_m3_s_per_pa=np.linalg.svd(B,compute_uv=False).tolist(),basis_step_pa=1.)


def fixed_path_identity(cache,Qin=ROI_TARGET_Q_M3_S):
    adj = [[] for _ in cache.node_ids]
    for e,(u,v) in enumerate(cache.edge_nodes): adj[u].append((v,e)); adj[v].append((u,e))
    root = int(cache.port_node_indices[0]); parent = {root:None}; pedge = {}; stack = [root]
    while stack:
        u = stack.pop()
        for v,e in adj[u]:
            if v not in parent: parent[v] = u; pedge[v] = e; stack.append(v)
    def path(n):
        nodes = [int(n)]
        while parent[nodes[-1]] is not None: nodes.append(parent[nodes[-1]])
        return nodes[::-1]
    p1,p3 = path(cache.port_node_indices[1]),path(cache.port_node_indices[3])
    j = [a for a,b in zip(p1,p3) if a == b][-1]; paths = {}
    for name,p in [('O1',p1),('O3',p3)]:
        nodes = p[p.index(j):]
        if any(len(adj[n]) != 2 for n in nodes[1:-1]): raise ValueError('Not a series arm')
        edges = [pedge[n] for n in nodes[1:]]
        paths[name] = dict(node_ids=cache.node_ids[nodes].tolist(),local_edge_indices=edges,
                           resistance_pa_s_m3=float(cache.edge_resistance_pa_s_m3[edges].sum()))
    return dict(junction_node_id=int(cache.node_ids[j]),junction_xyz_m=cache.xyz_m[j].tolist(),paths=paths,
                equal_split_dP_O1_O3_pa=(paths['O3']['resistance_pa_s_m3']-paths['O1']['resistance_pa_s_m3'])*Qin/3)


def design(cache, regression, certificate, Qin=ROI_TARGET_Q_M3_S):
    if regression['status'] != 'PASS': raise ValueError('Regression must pass before design')
    basis = pressure_basis(cache,Qin); B = np.asarray(basis['B_m3_s_per_pa'])
    x = np.linalg.lstsq(B,np.full(3,Qin/3)-basis['q0_m3_s'],rcond=None)[0]
    trace = []
    def residual(dp):
        r = solve_roi_mixed_bc(cache,[*dp,0.],Qin); f = r['outlet_fraction']
        trace.append(dict(evaluation=len(trace)+1,dP_O1_O3_pa=float(dp[0]),dP_O2_O3_pa=float(dp[1]),
                          **dict(zip(OUTLETS,f.tolist())),**metrics(f)))
        return f-1/3
    initial = np.asarray(regression['reference_port_pressure_pa'])[1:3]-regression['reference_port_pressure_pa'][3]
    opt = least_squares(residual,initial,jac='3-point',x_scale=1000.,ftol=1e-13,xtol=1e-13,gtol=1e-14)
    if not opt.success or max(abs(opt.x-x)) > 1e-6: raise ValueError('STOP: basis/iterative cross-check failed')
    identity = fixed_path_identity(cache,Qin)
    for p in ('O1','O3'):
        np.testing.assert_allclose(identity['paths'][p]['resistance_pa_s_m3'],certificate['paths'][p]['resistance_pa_s_m3'],rtol=1e-12)
    np.testing.assert_allclose(x[0],identity['equal_split_dP_O1_O3_pa'],rtol=1e-10,atol=1e-7)
    relative = np.r_[x,0.]; C = float(-min(relative))
    solved = solve_roi_mixed_bc(cache,relative,Qin); shifted = solve_roi_mixed_bc(cache,relative+C,Qin)
    gauge_error = float(max(abs(shifted['port_flow_m3_s']-solved['port_flow_m3_s']))/Qin)
    np.testing.assert_allclose(shifted['pressure_pa'],solved['pressure_pa']+C,rtol=1e-11,atol=1e-7)
    np.testing.assert_allclose((relative+C)[:,None]-(relative+C),relative[:,None]-relative,rtol=1e-12,atol=1e-12)
    if gauge_error > 1e-9: raise ValueError('STOP: gauge flow invariance failed')
    met = metrics(solved['outlet_fraction'])
    solution = dict(model_name=MODEL_NAME,workflow_kind='ROI_BOUNDARY_DESIGN',
        status='ROI_EQUAL_SPLIT_REACHED' if met['max_abs_deviation_from_one_third']<1e-8 else 'ROI_BEST_BALANCE_FOUND',
        full_A_network_used_for_boundary_generation=False,ROI_geometry_source='existing verified A ROI internal edges only',
        pressure_interpretation='Numerical boundary design; not measured physiology or downstream-network prediction',
        Qin_target_m3_s=Qin,dynamic_viscosity_pa_s=cache.mu_pa_s,
        relative_realcut_pressure_pa=dict(zip(OUTLETS,relative.tolist())),gauge_shift_pa=C,
        gauge_realcut_pressure_pa=dict(zip(OUTLETS,(relative+C).tolist())),
        inlet_pressure_relative_pa=float(solved['port_pressure_pa'][0]),inlet_pressure_pa=float(shifted['port_pressure_pa'][0]),
        outlet_flow_m3_s=dict(zip(OUTLETS,solved['port_flow_m3_s'][1:].tolist())),
        outlet_fraction=dict(zip(OUTLETS,solved['outlet_fraction'].tolist())),**met,mass_audit=solved['mass_audit'],
        geometry_sha256=cache.geometry_sha256,geometry_file_sha256=cache.geometry_file_sha256,
        port_mapping_sha256=cache.port_mapping_sha256,gauge_flow_error_relative_Qin=gauge_error,
        direct_relative_pressure_pa=x.tolist(),iterative_relative_pressure_pa=opt.x.tolist(),
        direct_minus_iterative_pa=(x-opt.x).tolist(),iterative_nfev=opt.nfev,iterative_evaluations=len(trace),
        iterative_message=opt.message,initial_relative_pressure_pa=initial.tolist(),
        fixed_path_identity=identity,regression_status=regression['status'],CFD_feedback_used=False,particle_RBC_calls=0)
    return solution,basis,trace


def content_hash(payload):
    return hashlib.sha256(json.dumps(payload,sort_keys=True,allow_nan=False,separators=(',',':')).encode()).hexdigest()


def freeze_design(path,solution):
    if solution['regression_status'] != 'PASS' or solution['mass_audit']['status'] != 'PASS': raise ValueError('Failed design gate')
    payload = dict(solution); payload['artifact_sha256'] = content_hash(solution)
    with Path(path).open('x') as f: yaml.safe_dump(payload,f,sort_keys=False,allow_unicode=True)
    return load_frozen_design(path)


def load_frozen_design(path):
    p = yaml.safe_load(Path(path).read_text()); digest = p.pop('artifact_sha256')
    if digest != content_hash(p): raise ValueError('Frozen ROI design modified')
    if p['model_name'] != MODEL_NAME or p['workflow_kind'] != 'ROI_BOUNDARY_DESIGN' or p['full_A_network_used_for_boundary_generation'] is not False:
        raise ValueError('Not an independent ROI design')
    p['artifact_sha256'] = digest
    return p
