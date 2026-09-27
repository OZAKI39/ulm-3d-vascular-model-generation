"""Read-only provenance audit of the current A/ROI case; no legacy flow reuse."""
import ast
import csv
import importlib.util
import json
from pathlib import Path
import time
import resource
import platform
import numpy as np
import yaml
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .audit import sha256, write_csv, write_json, terminal_audit
from .swc_graph import load_swc, audit_graph
from .roi_mapping import match_edge_position, insert_cuts, clip_interval
from .hydraulic_resistance import linear_radius_resistance, um_to_m, MU_PA_S, RHO_KG_M3, ROI_TARGET_Q_M3_S
from .boundary_conditions import require_source, SourceAmbiguous
from .network_solver import conductance_matrix


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def read_csv(path):
    with open(path, encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def remap_historical(path, root):
    """Use the literal project-relative suffix saved in Windows-era metadata."""
    normalized = str(path).replace('\\', '/').replace('//', '/')
    marker = 'ulm_3D_vascular/'
    if marker in normalized:
        return root/normalized.split(marker, 1)[1]
    return Path(path)


def source_inventory(root):
    """Read the local import closure without executing historical pipelines."""
    queue = [root/'s1-1_swc_roi_generate_mouse.py', root/'s2_swc_stl_model_generate.py']
    seen = set(); records = []
    while queue:
        path = queue.pop()
        if path in seen or not path.is_file(): continue
        seen.add(path); code = path.read_text(encoding='utf-8-sig'); tree = ast.parse(code)
        relative = path.relative_to(root); module = list(relative.with_suffix('').parts)
        package = module[:-1]
        imports = []
        for node in ast.walk(tree):
            candidates = []
            if isinstance(node, ast.Import): candidates = [x.name for x in node.names]
            elif isinstance(node, ast.ImportFrom):
                base = package[:len(package)-node.level+1] if node.level else []
                if node.module: base += node.module.split('.')
                candidates = ['.'.join(base)] + ['.'.join(base+[x.name]) for x in node.names if x.name != '*']
            for name in candidates:
                if not name.startswith('utils'): continue
                target = root/Path(*name.split('.'))
                for f in [target.with_suffix('.py'), target/'__init__.py']:
                    if f.is_file(): imports.append(str(f.relative_to(root))); queue.append(f)
        records.append(dict(path=str(relative), sha256=sha256(path), lines=len(code.splitlines()),
                            local_imports=sorted(set(imports)), functions=[dict(name=n.name, line=n.lineno) for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]))
    return sorted(records, key=lambda x: x['path'])


def graph_rows(g):
    nodes = [dict(original_id=int(g.ids[i]), x_um=g.xyz_um[i, 0], y_um=g.xyz_um[i, 1], z_um=g.xyz_um[i, 2],
                  radius_um=g.radius_um[i], parent=int(g.parent[i]), degree=int(g.degree[i]), component_id=int(g.component[i])) for i in range(len(g.ids))]
    edges = [dict(edge_id=i, u=int(g.ids[u]), v=int(g.ids[v]), length_um=g.length_um[i], r_u_um=g.radius_um[u], r_v_um=g.radius_um[v], component_id=int(g.component[u])) for i, (u, v) in enumerate(g.edges)]
    return nodes, edges


def run_audit(root, output, flow_root, fem_root):
    start = time.perf_counter()
    output = Path(output); data = output/'data'; data.mkdir(parents=True, exist_ok=True)
    (output/'logs').mkdir(exist_ok=True); (output/'figures').mkdir(exist_ok=True)
    flow = Path(flow_root)/'formal_3D_flow_solver/FEM_SimVascular'
    case = flow/'flow_cases/mean-2p0-mmps'
    port_contract_path = flow/'inputs/fem_reference/port_contract.json'
    pc = read_json(port_contract_path)
    source_contract_path = Path(fem_root)/pc['source_contract_path']
    contract = read_json(source_contract_path); lineage = contract['lineage']
    if sha256(source_contract_path) != pc['source_contract_sha256']:
        raise ValueError('FEM source contract SHA mismatch')
    rodent = Path(lineage['rodent_run']); sampling = Path(lineage['sampling_run'])
    recon = Path(lineage['reconstruction_run']); preprocess = Path(lineage['preprocess_run'])
    source_surface_run = Path(contract['source_run'])
    sample = Path(lineage['analysis_swc']).parent.parent
    saved_config_path = rodent/'source_swc_roi_generate.yaml'
    saved = yaml.safe_load(saved_config_path.read_text()); current = yaml.safe_load((root/'configs/swc_roi_generate.yaml').read_text())
    manifest = read_json(sample/'preprocess_manifest.json')
    original = remap_historical(manifest['record']['swc_path'], root)
    if original != Path(lineage['source_swc']): raise ValueError('Two provenance chains disagree about original A')
    input_dir = root/saved['paths']['input_dir']
    if not original.is_relative_to(input_dir) or original.stem != saved['pipeline']['sample_id']:
        raise ValueError('Saved s1 configuration does not select the original file')
    scale = np.array(manifest['spacing_xyz_um'])
    g = load_swc(original, scale); ga = audit_graph(g)
    normalized_path = sample/'normalized/analysis_swc_single_component.npz'
    norm = np.load(normalized_path, allow_pickle=False)
    ix = g.index; indices = np.array([ix[int(n)] for n in norm['node_ids']])
    normalized_checks = dict(ids_parent_exact=np.array_equal(g.parent[indices], norm['parent_ids']),
                             xyz_exact=np.array_equal(g.xyz_um[indices], norm['points_um']),
                             radius_exact=np.array_equal(g.radius_um[indices], norm['radius_raw_um']))
    if not all(normalized_checks.values()): raise ValueError('Normalized A differs from raw source')
    component = int(g.component[indices[0]])
    if set(indices) != set(np.flatnonzero(g.component == component)):
        raise ValueError('Analysis SWC is not exactly one original component')
    nodes, edges = graph_rows(g)
    write_csv(data/'a_nodes.csv', nodes); write_csv(data/'a_edges.csv', edges)
    write_json(data/'a_graph_audit.json', ga)
    comp_rows = []
    for c in np.unique(g.component):
        ni = np.flatnonzero(g.component == c); ei = np.flatnonzero(g.component[g.edges[:, 0]] == c)
        comp_rows.append(dict(component_id=int(c), nodes=len(ni), edges=len(ei), cycle_rank=len(ei)-len(ni)+1,
                             total_length_um=float(g.length_um[ei].sum()), selected=bool(c == component),
                             structural_roots=';'.join(map(str, g.ids[ni][g.parent[ni] == -1]))))
    write_csv(data/'a_components.csv', comp_rows)
    terminals = terminal_audit(g, manifest['auxiliary_volume_shape_zyx'][::-1], scale, component)
    write_csv(data/'a_terminal_audit.csv', terminals)
    # Preserve every original resistance, not the Ultraliser feed radius multiplier.
    r = linear_radius_resistance(um_to_m(g.length_um), um_to_m(g.radius_um[g.edges[:, 0]]), um_to_m(g.radius_um[g.edges[:, 1]]))
    G = conductance_matrix(len(g.ids), g.edges, r)
    write_json(data/'network_solver_audit.json', dict(matrix_size=list(G.shape), nnz=G.nnz,
               representation='scipy.sparse.csr_matrix', solve_executed=False, solver='scipy.sparse.linalg.splu (tested on analytic networks only)',
               condition_number=None, residual=None, reason='Full A has 43 components; no hydraulic source or certified terminal Dirichlet set has been assigned. Ungrounded Laplacian is singular by construction.'))
    write_csv(data/'a_edge_hydraulic_resistance.csv', [row | dict(length_m=row['length_um']*1e-6, r_u_m=row['r_u_um']*1e-6, r_v_m=row['r_v_um']*1e-6, R_Pa_s_m3=r[i]) for i, row in enumerate(edges)])
    roi = np.load(lineage['roi_archive'], allow_pickle=False)
    metadata = read_json(recon/'input/metadata.json')
    global_edges_path = sampling/'manifests/global_edges.csv'
    global_edges = {int(row['global_edge_id']): (int(row['upstream_global_node_id']), int(row['downstream_global_node_id'])) for row in read_csv(global_edges_path) if row['source_model_id'] == str(roi['source_model_id'])}
    raw_eid = {(int(g.ids[u]), int(g.ids[v])): eid for eid, (u, v) in enumerate(g.edges)}
    # Compare all saved ROI node and edge geometry, not only the four ports.
    canonical = np.loadtxt(lineage['canonical_roi_swc']); canonical_by_id = {int(row[0]): row for row in canonical}
    roi_node_rows = []
    for local, original_id in enumerate(roi['local_node_global_ids']):
        saved_xyz = roi['local_node_positions_um'][local]; saved_radius = roi['local_node_radius_um'][local]
        core_id = metadata['swc_node_id_by_local_node_id'][str(local)]
        core = canonical_by_id[core_id]
        if np.linalg.norm(core[2:5]-saved_xyz) > 1e-8 or abs(core[5]-saved_radius) > 1e-8:
            raise ValueError('Canonical ROI SWC geometry mismatch')
        if original_id >= 0:
            i = ix[int(original_id)]
            if not np.array_equal(g.xyz_um[i], saved_xyz) or g.radius_um[i] != saved_radius:
                raise ValueError('ROI source-node mismatch')
        roi_node_rows.append(dict(local_id=local, canonical_swc_id=core_id, original_id=int(original_id) if original_id >= 0 else '', virtual=bool(original_id < 0), x_um=saved_xyz[0], y_um=saved_xyz[1], z_um=saved_xyz[2], radius_um=saved_radius))
    write_csv(data/'roi_node_provenance.csv', roi_node_rows)
    roi_edge_rows = []; retained = {}
    for local, global_id in enumerate(roi['local_edge_global_ids']):
        u_id, v_id = global_edges[int(global_id)]; u, v = ix[u_id], ix[v_id]; ts = []
        for point, radius in zip(roi['local_edge_points_um'][local], roi['local_edge_radius_um'][local]):
            ts.append(match_edge_position(g.xyz_um[u], g.xyz_um[v], g.radius_um[u], g.radius_um[v], point, radius)['fraction_parent_to_child'])
        if ts[1] <= ts[0]: raise ValueError('Reversed ROI edge')
        retained[raw_eid[(u_id, v_id)]] = ts
        roi_edge_rows.append(dict(local_edge_id=local, analysis_global_edge_id=int(global_id), raw_edge_id=raw_eid[(u_id, v_id)], u=u_id, v=v_id, start_fraction=ts[0], end_fraction=ts[1]))
    write_csv(data/'roi_edge_provenance.csv', roi_edge_rows)
    # Independent scan of every original edge through the box.
    crossing_rows = []
    for eid, (u, v) in enumerate(g.edges):
        interval = clip_interval(g.xyz_um[u], g.xyz_um[v], roi['bbox_min_um'], roi['bbox_max_um'])
        if interval is None: continue
        for end, t in zip(('entry', 'exit'), interval):
            if 1e-12 < t < 1-1e-12:
                xyz = g.xyz_um[u]+t*(g.xyz_um[v]-g.xyz_um[u])
                crossing_rows.append(dict(raw_edge_id=eid, u=int(g.ids[u]), v=int(g.ids[v]), fraction=t, end=end,
                                          x_um=xyz[0], y_um=xyz[1], z_um=xyz[2], component_id=int(g.component[u]),
                                          belongs_to_saved_roi=eid in retained))
    write_csv(data/'all_A_roi_bbox_crossings.csv', crossing_rows)
    expected_crossings = sum(row['belongs_to_saved_roi'] for row in crossing_rows)
    # A real retained node with an incident unsaved edge would be an extra port.
    real_ids = set(int(x) for x in roi['local_node_global_ids'] if x >= 0)
    missing_incident = [eid for eid, (u, v) in enumerate(g.edges) if eid not in retained and (int(g.ids[u]) in real_ids or int(g.ids[v]) in real_ids)]
    local_degree = np.bincount(roi['local_edges'].ravel(), minlength=len(roi['local_node_ids']))
    if expected_crossings != 3 or missing_incident or int(np.sum(local_degree == 1)) != 4:
        raise ValueError('ROI_PORT_SET_INCOMPLETE')
    classification_path = preprocess/'roi/port_classification.csv'
    upstream_ports = read_csv(classification_path)
    extension_geometry_path = source_surface_run/'qc/extension_geometry_qc.json'
    extension_pressure_path = source_surface_run/'qc/extension_pressure_geometry_qc.json'
    eg = {x['port_id']: x for x in read_json(extension_geometry_path)['boundaries']}
    ep = {x['port_id']: x for x in read_json(extension_pressure_path)['rows']}
    names = ['INLET', 'O1', 'O2', 'O3']; fem_names = ['inlet', 'outlet_01', 'outlet_02', 'outlet_03']
    ports = []; cuts = {}
    # Read actual production mesh cap triangles, not just the old surface cap.
    import pyvista as pv
    for name, fem_name, row in zip(names, fem_names, upstream_ports):
        u_id, v_id = global_edges[int(row['global_edge_id'])]; u, v = ix[u_id], ix[v_id]
        xyz = np.array([float(row[a+'_um']) for a in 'xyz']); radius = float(row['radius_um'])
        match = match_edge_position(g.xyz_um[u], g.xyz_um[v], g.radius_um[u], g.radius_um[v], xyz, radius)
        outward = np.array([float(row['outward_normal_'+a]) for a in 'xyz'])
        mesh_name = 'INLET' if name == 'INLET' else 'OUTLET_0'+name[1:]
        mesh = pv.read(case/'SV_MESH/mesh-surfaces'/f'{mesh_name}.vtp')
        tris = mesh.faces.reshape(-1, 4)[:, 1:]; points = np.asarray(mesh.points)
        vec = np.cross(points[tris[:, 1]]-points[tris[:, 0]], points[tris[:, 2]]-points[tris[:, 0]])/2
        areas = np.linalg.norm(vec, axis=1); centroid = np.average(points[tris].mean(1), axis=0, weights=areas)
        axial_um = float(np.dot(centroid*1e6-xyz, outward))
        eq = eg[row['port_id']]; pressure_profile = ep[row['port_id']]
        fractions = np.array(pressure_profile['station_fractions']); radii = um_to_m(pressure_profile['equivalent_radius_um'])
        # Extrapolate only the end half-bins as constants; integrate interior linear-radius profile exactly.
        fractions = np.r_[0., fractions, 1.]; radii = np.r_[radii[0], radii, radii[-1]]
        L = axial_um*1e-6
        estimate = float(linear_radius_resistance(np.diff(fractions)*L, radii[:-1], radii[1:]).sum())
        virtual_id = -10001-len(ports) if row['boundary_origin'] == 'CUT_PORT' else int(row['global_node_id'])
        if row['boundary_origin'] == 'CUT_PORT': cuts[raw_eid[(u_id, v_id)]] = [(match['fraction_parent_to_child'], virtual_id)]
        port = dict(name=name, port_id=row['port_id'], boundary_origin=row['boundary_origin'], original_edge=[u_id, v_id],
                    raw_edge_id=raw_eid[(u_id, v_id)], analysis_global_edge_id=int(row['global_edge_id']),
                    original_node_id=int(row['global_node_id']) if row['global_node_id'] else None,
                    network_node_id=virtual_id, real_cut_xyz_um=xyz, raw_coordinate_xyz=xyz/scale, radius_um=radius,
                    outward_normal=outward, parent_direction_dot_outward=float(np.dot(match['tangent_parent_to_child'], outward)),
                    orientation_basis='SWC_PARENT_TO_CHILD_STRUCTURAL_NOT_MEASURED',
                    roi_positive_flow_relative_to_parent=1, **match,
                    fem_cap=dict(centroid_um=centroid*1e6, area_m2=float(areas.sum()), equivalent_radius_um=float(np.sqrt(areas.sum()/np.pi)*1e6),
                                 real_cut_to_cap_axial_um=axial_um, plane_offset_error_m=float(np.max(abs((points-np.array(pc['ports'][fem_name]['plane_origin_m'])) @ pc['ports'][fem_name]['outward_normal'])))),
                    extension=dict(kind='FEM_ARTIFICIAL_EXTENSION', planned_length_um=eq['planned_extension_length_um'],
                                   saved_actual_axial_um=eq['actual_axial_length_um'], proximal_equivalent_radius_um=eq['equivalent_radius_original_um'],
                                   distal_equivalent_radius_um=eq['equivalent_radius_distal_um'],
                                   historical_20_station_R_Pa_s_m3=pressure_profile['extension_resistance_pa_s_m3'],
                                   current_length_piecewise_linear_station_radius_R_Pa_s_m3=estimate,
                                   estimate_limitation='Historical surface equivalent circular cross-sections; 20 saved samples, constant half-bin endpoints; current SV cap areas differ and current extension profile is not remeasured',
                                   historical_integration_fraction_interval=[float(fractions[1]), float(fractions[-2])],
                                   historical_trapezoid_omitted_end_fraction=.05,
                                   recomputed_to_historical_R_ratio=estimate/pressure_profile['extension_resistance_pa_s_m3'],
                                   pressure_transfer='p_cap = p_real - R_extension * Q_outward', pressure_drop_Pa=None))
        ports.append(port)
    mapping = dict(status='PASS', ports=ports, mapping_method='saved original IDs and exact edge interpolation; no nearest-neighbor',
                   roi_node_count=len(roi['local_node_ids']), original_node_count=len(real_ids), roi_edges=len(retained),
                   true_graph_terminal_count=1, graph_cut_count=expected_crossings, expected_fem_ports=4,
                   extra_saved_roi_ports=0, missing_incident_edges=missing_incident,
                   bbox_all_raw_crossings=len(crossing_rows), bbox_analysis_crossings=sum(x['component_id'] == component for x in crossing_rows),
                   bbox_other_crossings_are_not_roi_ports=True)
    write_json(data/'roi_ports_in_a.json', mapping)
    split = insert_cuts(g.ids, um_to_m(g.xyz_um), um_to_m(g.radius_um), g.edges, cuts)
    np.savez_compressed(data/'a_graph_with_exact_cuts_si.npz', **split)
    is_internal = np.array([eid in retained and a >= retained[eid][0]-1e-12 and b <= retained[eid][1]+1e-12 for eid, (a, b) in zip(split['original_edge'], split['fractions'])])
    ext_e = split['edges'][~is_internal]; u, v = ext_e.T; n = len(split['ids'])
    adj = coo_matrix((np.ones(2*len(ext_e)), (np.r_[u, v], np.r_[v, u])), shape=(n, n))
    _, labels = connected_components(adj, directed=False)
    port_indices = [int(np.flatnonzero(split['ids'] == p['network_node_id'])[0]) for p in ports]
    conn = labels[port_indices, None] == labels[port_indices]
    ext_degrees = np.bincount(ext_e.ravel(), minlength=n)[port_indices]
    topo = dict(status='TOPOLOGY_ONLY_NO_HYDRAULIC_BC', port_order=names, external_same_component=conn,
                exterior_degree=ext_degrees, removed_roi_segment_count=int(is_internal.sum()),
                external_segments=int((~is_internal).sum()), O3_external_edge_count=int(ext_degrees[3]),
                interpretation='Distinct exterior port components; O3 has no exterior vessel. No finite downstream R3 can be inferred without terminal closure.',
                multiport_Y=None, multiport_Z=None, coupling_metric=None)
    write_json(data/'external_network_topology.json', topo)
    roots = norm['root_ids'].tolist()
    source = dict(structural_root_ids=roots, structural_root_xyz_um=[g.xyz_um[ix[int(x)]].tolist() for x in roots],
                  hydraulic_source_id=None, identification_kind=None, evidence_reference=None,
                  legacy_root_policy='single_structural_root', legacy_role='ASSUMED_GLOBAL_INLET',
                  reason='The saved pipeline explicitly assumes SWC parent-to-current flow. It identifies serialization root 2410, not a dataset-labelled hydraulic source. Source identity remains unresolved under this task source gate.',
                  legacy_evidence_path=str(preprocess/'config/source_cfd_preprocess.yaml'))
    try: require_source(source)
    except SourceAmbiguous: status = 'A_NETWORK_SOURCE_AMBIGUOUS'
    else: raise RuntimeError('Source evidence changed: this audit-only case requires a new reviewed solve configuration')
    write_json(data/'source_identification_audit.json', source)
    radius_sensitivity = dict(status='HYDRAULIC_SENSITIVITY_BLOCKED_BY_SOURCE_GATE', empirical_uncertainty='NO_EMPIRICAL_RADIUS_UNCERTAINTY_AVAILABLE',
                              global_scaling=[dict(radius_factor=f, exact_R_factor=f**-4, conditional_fixed_target_pressure_factor=f**-4,
                                           predicted_split_change_for_uniform_scaling='exactly zero for unchanged linear homogeneous boundary model', solved=False) for f in [.9, .95, 1., 1.05, 1.1]],
                              branchwise_status='DETERMINISTIC_RADIUS_PATTERNS_SAVED_ONLY; no pressure or flow sensitivity computed')
    # Deterministic branch groups = maximal chains between non-degree-2 vertices.
    incident = [[] for _ in g.ids]
    for eid, (u, v) in enumerate(g.edges): incident[u].append(eid); incident[v].append(eid)
    groups = np.full(len(g.edges), -1, dtype=int); group = 0
    for start_node in np.r_[np.flatnonzero(g.degree != 2), np.arange(len(g.ids))]:
        for first in incident[start_node]:
            if groups[first] >= 0: continue
            node, edge = int(start_node), first
            while groups[edge] < 0:
                groups[edge] = group; u, v = g.edges[edge]; node = int(v if u == node else u)
                if g.degree[node] != 2: break
                nxt = [e for e in incident[node] if groups[e] < 0]
                if not nxt: break
                edge = nxt[0]
            group += 1
    node_group = [min((groups[e] for e in incident[i]), default=-1) for i in range(len(g.ids))]
    factors = np.array([1+.05*(1 if k % 2 == 0 else -1) for k in node_group])
    write_csv(data/'radius_branch_patterns.csv', [dict(original_id=int(g.ids[i]), branch_group=int(node_group[i]), plus_pattern_factor=factors[i], minus_pattern_factor=2-factors[i],
                                                        plus_radius_um=g.radius_um[i]*factors[i], minus_radius_um=g.radius_um[i]*(2-factors[i])) for i in range(len(g.ids))])
    radius_sensitivity['branch_group_count'] = group
    radius_sensitivity['junction_policy'] = 'Each raw node uses the lowest incident branch ID; one radius per junction, both deterministic +/- patterns. Edge profiles remain continuous.'
    radius_sensitivity['branchwise_resistance_only'] = []
    for name, f in [('plus_pattern', factors), ('minus_pattern', 2-factors)]:
        perturbed = g.radius_um*f
        rp = linear_radius_resistance(um_to_m(g.length_um), um_to_m(perturbed[g.edges[:, 0]]), um_to_m(perturbed[g.edges[:, 1]]))
        radius_sensitivity['branchwise_resistance_only'].append(dict(pattern=name, min_R_ratio=float(np.min(rp/r)), max_R_ratio=float(np.max(rp/r)),
                                                                    network_pressure_scale=None, flow_split=None))
    write_json(data/'radius_sensitivity.json', radius_sensitivity)
    terminal_sensitivity = dict(status='NOT_SOLVED_SOURCE_AND_TERMINAL_SETS_UNRESOLVED', MODEL_A=dict(trusted_outer_terminals=[], pressure_Pa=0, valid=False),
                                MODEL_B=dict(strict_outer_terminals=[], uncertain_policy='NO_FLOW', valid=False),
                                reason='No certified outer-terminal list; closing every uncertain endpoint leaves no valid sink. Neither model can be called a validated baseline.',
                                geometric_candidate_counts=dict(radius_ball_reaches_extent=sum(x['radius_ball_reaches_image_extent'] for x in terminals),
                                                                one_voxel_to_extent=sum(x['within_one_voxel_of_image_extent'] for x in terminals),
                                                                deep_interior=sum(x['deep_interior_candidate'] for x in terminals)))
    write_json(data/'terminal_sensitivity.json', terminal_sensitivity)
    imported = source_inventory(root); write_json(data/'source_helper_inventory.json', imported)
    protected = [original, normalized_path, Path(lineage['roi_archive']), Path(lineage['canonical_roi_swc']), saved_config_path,
                 recon/'input/source_swc_stl_model_generate.yaml', source_contract_path, port_contract_path, global_edges_path,
                 classification_path, extension_geometry_path, extension_pressure_path, Path(contract['geometry_path']),
                 recon/'geometry/lumen_surface_um.stl', case/'frozen_flow/steady_flow_mean_2p0_mmps.vtu',
                 case/'FLOW_2MMPS_VALIDATION.json', flow/'inputs/MANIFEST.json', root/'s1-1_swc_roi_generate_mouse.py', root/'s2_swc_stl_model_generate.py',
                 root/'configs/swc_roi_generate.yaml', root/'configs/swc_stl_model_generate.yaml',
                 flow/'frozen_reference/mesh_manifest.json', flow/'scripts/flow_2mmps/prepare.py',
                 root/'utils/cfd_surface_prepare/vmtk_qc.py', flow/'configs/mesh_policy.json',
                 flow/'scripts/sv_generate_mesh.py', flow/'scripts/prepare_surface.py',
                 flow/'reports/sv1/geometry_qc.json', flow/'reports/sv1/mesh_validity.json']
    input_hashes = read_json(case/'input_hashes.json')
    hash_checks = {path: sha256(case/path) == expected for path, expected in input_hashes.items()}
    if not all(hash_checks.values()): raise ValueError('Existing production case differs from its frozen hash manifest')
    protected += [case/path for path in input_hashes]
    hashes = {str(p): sha256(p) for p in protected}
    write_json(data/'protected_input_hashes.json', hashes)
    provenance = dict(original_A_path=str(original), original_A_sha256=sha256(original), original_audit=ga,
                      definition='Original release SWC = all 43 components; historical displayed analysis A = selected component 42. Never equate the two silently.',
                      analysis_component=comp_rows[component], normalized_checks=normalized_checks,
                      saved_config=str(saved_config_path), saved_config_content=saved,
                      current_default_input=current['paths']['input_dir'], current_default_differs_from_historical=current['paths']['input_dir'] != saved['paths']['input_dir'],
                      coordinate_transform=dict(raw_xyz_unit='voxel', spacing_xyz_um=scale, radius_raw_unit='um', physical_xyz='raw_xyz * [1,1,2]', hydraulic_SI='physical_um * 1e-6',
                                                normalized_text_swc_warning='Historical analysis .swc stores voxel xyz; NPZ points_um is physical. Do not read text xyz as um or apply spacing twice.'),
                      lineage=lineage, roi_anchor=int(roi['anchor_id']), roi_center_um=roi['anchor_position_um'], roi_bbox_um=[roi['bbox_min_um'], roi['bbox_max_um']],
                      roi_raw_node_ids=sorted(real_ids), roi_raw_edge_ids=sorted(retained), analysis_edge_ids=roi['local_edge_global_ids'],
                      roi_node_reindex_map=metadata['swc_node_id_by_local_node_id'], ids_preserved_until_canonical_roi_swc=True,
                      raw_swc_geometry_modified=False, clipping='Exact segment/AABB intersection, linear radius interpolation; keep anchor connected component only',
                      reconstruction_metadata=metadata, reconstruction_qc=read_json(recon/'qc/run_summary.json'),
                      resampling='Graph display derived branch geometry is resampled at 1 um; ROI sampling reloads raw normalized SWC, not derived branches',
                      surface_extension='Official VMTK thinplatespline, boundarynormal, local collar remesh, distal cap; full records hashed',
                      fem_source_contract=str(source_contract_path), fem_mesh_provenance=read_json(flow/'inputs/MANIFEST.json'),
                      actual_sv_mesh_manifest=read_json(flow/'frozen_reference/mesh_manifest.json'),
                      sv_mesh_generation_policy=read_json(flow/'configs/mesh_policy.json'),
                      sv_mesh_geometry_changes=read_json(flow/'reports/sv1/geometry_qc.json'),
                      sv_mesh_validity=read_json(flow/'reports/sv1/mesh_validity.json'),
                      current_case_preparation='scripts/flow_2mmps/prepare.py copies frozen_reference/SV_MESH byte-for-byte; changes inlet Q only',
                      current_mesh=read_json(case/'FLOW_2MMPS_VALIDATION.json')['mesh'], input_hash_checks=hash_checks,
                      source_root=source, helper_inventory='source_helper_inventory.json', physical_boundary_metadata_found=False,
                      dataset_sources=['https://www.frontiersin.org/journals/neuroinformatics/articles/10.3389/fninf.2026.1809341/full', 'https://zenodo.org/records/19184697'])
    write_json(data/'a_network_provenance.json', provenance)
    summary = dict(status=status, original_A_path=str(original), original_A_sha256=sha256(original),
                   roi_source_script_sha256=sha256(root/'s1-1_swc_roi_generate_mouse.py'), roi_stl_script_sha256=sha256(root/'s2_swc_stl_model_generate.py'),
                   A_nodes=len(g.ids), A_edges=len(g.edges), A_components=ga['components'], A_cycles=ga['cycle_rank'],
                   analysis_A=comp_rows[component], source_root=source,
                   terminal_counts_by_class=dict(SOURCE_ROOT=0, OUTER_BOUNDARY_TERMINAL=0, INTERIOR_DEAD_END=0, UNKNOWN_TERMINAL=len(terminals)),
                   roi_ports=ports, roi_target_Q=ROI_TARGET_Q_M3_S, mu=MU_PA_S, rho=RHO_KG_M3,
                   full_A_source_pressure_scaled=None, full_A_total_inflow=None, ROI_inlet_pressure=None, ROI_inlet_flow=None,
                   O1_pressure=None, O1_flow=None, O1_fraction=None, O2_pressure=None, O2_flow=None, O2_fraction=None,
                   O3_pressure=None, O3_flow=None, O3_fraction=None, current_ROI_zero_pressure_split=[.04257917, .85205026, .10537057],
                   split_difference=None, independent_R_if_valid=None, multiport_Y=None, multiport_Z=None, coupling_metric=None,
                   radius_sensitivity=radius_sensitivity, terminal_sensitivity=terminal_sensitivity,
                   blocker='Hydraulic source not identified by dataset; structural root known. All endpoint roles remain uncertified.',
                   provenance_status='PASS', roi_mapping_status='PASS', roi_port_set_status='PASS', numerical_tests_status='SEE_TEST_LOG',
                   full_A_solve_status='NOT_RUN_SOURCE_GATE', network_reduction_status='TOPOLOGY_ONLY_SOURCE_AND_TERMINAL_GATE',
                   runtime_seconds=time.perf_counter()-start, peak_memory_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                   workers=1, execution_host=platform.node(), server_used=False)
    summary['3D_validation_status'] = 'NOT_CREATED_NOT_RUN_PRECONDITIONS_FAILED'
    summary['3D_validation_split'] = None
    write_json(data/'final_summary.json', summary)
    write_json(data/'network_reduction_audit.json', topo | dict(status='NOT_COMPUTED_SOURCE_AND_TERMINAL_GATE', symmetry_error=None, condition_number=None, passive_check=None, full_reduced_equivalence=None))
    write_json(data/'roi_fixed_pressure_bc_from_A.json', dict(status='BLOCKED_DO_NOT_USE', reference_pressure_definition='Proposed certified outer terminals 0 Pa gauge; terminal set unresolved',
               raw_pressure_Pa=None, shifted_pressure_Pa=None, pressure_differences=None, source_A_solution_SHA=None,
               roi_port_mapping_SHA=sha256(data/'roi_ports_in_a.json'), pressure_location='REAL_A_NETWORK_CUT; FEM cap needs p_real-R_ext Q_outward',
               extension_pressure_correction_Pa=None, blocked_by=status))
    # Empty files are schema-only; explicit manifest prevents interpreting no rows as zero.
    schemas = {'a_network_solution_raw.csv':['node_id', 'pressure_Pa'], 'a_edge_flow_solution.csv':['edge_id', 'Q_m3s', 'flow_direction', 'mean_velocity_estimate_m_s'],
               'roi_boundary_admittance_matrix.csv':['port', *names], 'roi_boundary_resistance_matrix.csv':['port', *names]}
    for name, fields in schemas.items(): write_csv(data/name, [], fields)
    write_json(data/'uncomputed_outputs.json', dict(reason=status, schema_only=list(schemas),
                    omitted_figures=['02_full_A_pressure_map', '03_full_A_flow_map', '04_current_vs_network_split', '05_outlet_coupling_matrix', '06_three_generation_split'],
                    rule='Missing physical model data are null/absent, never zero-valued invented solutions'))
    print(json.dumps({k: summary[k] for k in ['status', 'A_nodes', 'A_edges', 'A_components', 'runtime_seconds', 'peak_memory_MiB']}, indent=2))
    return summary
