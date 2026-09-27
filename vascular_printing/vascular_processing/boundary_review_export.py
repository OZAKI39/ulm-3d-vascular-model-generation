"""Read-only boundary evidence from frozen BG001 RMCA artifacts.

This module does not fit, vote, relabel, extract ROIs, or build vessel surfaces.
Branch polylines are decoded from the saved labeled VTP, and all classification
statistics are read from the saved CSV/NPZ rather than recomputed.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
import csv
import hashlib
import itertools
import json
import math
import builtins
import io
import os
from pathlib import Path
from unittest.mock import patch

import networkx as nx
import numpy as np

from .topbrain_qc import sha256, write_json

ROOT=Path(__file__).resolve().parents[1]
HUMAN_SHA='3361cea44a97e8bfa907ff1504b0f977e329d35eb59cb9f57cc1b3a2a55cd729'
RAW_SHA='5b69c6c406e4724ffc15a2a4c6b5958cbb169b822a05098d6ebedeecbbb3afbe'
COLOR_SHA='733b255aeeb7795de45799f64e73b19a93ada45ca63518e9aa479c203658ef0e'
TYPES=('PROXIMAL_M1_TO_MEVO','DISTAL_MEVO_TO_UNKNOWN','DISTAL_MEVO_TO_TERMINAL',
       'UNKNOWN_TO_MEVO','MEVO_TO_M1_REVERSE','MEVO_COMPONENT_ROOT','DETACHED_MEVO_COMPONENT')
TRANSITIONS={('M1','MeVO'):TYPES[0],('MeVO','UNKNOWN'):TYPES[1],
             ('UNKNOWN','MeVO'):TYPES[3],('MeVO','M1'):TYPES[4]}


class ReviewSourceMismatch(ValueError):
    def __init__(self,message):super().__init__('REVIEW_SOURCE_MISMATCH: '+message)


def read_csv(path):
    with Path(path).open(newline='',encoding='utf-8') as stream:return list(csv.DictReader(stream))


def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:raise ValueError('Cannot infer columns for an empty review table')
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def number(value):
    if value in (None,''):return None
    value=float(value)
    return value if math.isfinite(value) else None


@contextmanager
def protect_inputs(snapshot):
    """Reject Python file writes/unlinks to protected inputs before they happen."""
    protected={str(Path(p).resolve()) for p in snapshot}
    def check(path):
        if isinstance(path,(str,bytes,os.PathLike)) and str(Path(os.fsdecode(path)).resolve()) in protected:
            raise PermissionError('REVIEW_READ_ONLY_INPUT: '+str(path))
    original_open=builtins.open;original_io=io.open;original_os=os.open
    def wrap_open(original):
        def guarded(path,mode='r',*args,**kwargs):
            if any(letter in str(mode) for letter in 'wax+'):check(path)
            return original(path,mode,*args,**kwargs)
        return guarded
    def guarded_os(path,flags,*args,**kwargs):
        if flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND):check(path)
        return original_os(path,flags,*args,**kwargs)
    original_unlink=os.unlink;original_remove=os.remove;original_rename=os.rename;original_replace=os.replace
    def deletion(original):
        def guarded(path,*args,**kwargs):check(path);return original(path,*args,**kwargs)
        return guarded
    def movement(original):
        def guarded(source,destination,*args,**kwargs):check(source);check(destination);return original(source,destination,*args,**kwargs)
        return guarded
    with ExitStack() as stack:
        for name,value in [('builtins.open',wrap_open(original_open)),('io.open',wrap_open(original_io)),('os.open',guarded_os),
            ('os.unlink',deletion(original_unlink)),('os.remove',deletion(original_remove)),
            ('os.rename',movement(original_rename)),('os.replace',movement(original_replace))]:stack.enter_context(patch(name,value))
        yield


def stable_id(upstream,downstream,kind,subject='BG001',side='RMCA'):
    identity=json.dumps([subject,side,upstream,downstream,kind],separators=(',',':'))
    return f'{subject}_{side}_B'+hashlib.sha256(identity.encode()).hexdigest()[:10]


@contextmanager
def prohibit_recomputation():
    """Real export runs with all forbidden production entry points disabled."""
    names=[
        'vascular_processing.similarity_registration.register',
        'vascular_processing.semantic_ensemble.register',
        'vascular_processing.semantic_ensemble.select_donors',
        'vascular_processing.semantic_ensemble.nearest_baseline',
        'vascular_processing.semantic_ensemble.nearest_support',
        'vascular_processing.brava_branch_labels.aggregate_branches',
        'vascular_processing.brava_branch_labels.select_major_branches',
        'vascular_processing.brava_mevo_roi.strict_components',
        'vascular_processing.brava_mevo_roi.extract_and_model',
        'vascular_processing.brava_mevo_roi.process_file',
        'vascular_processing.pipeline.process_file',
        'vascular_processing.pipeline.save_native_surface',
        'third_party.vascularmd.ArterialTree.ArterialTree.model_network',
        'third_party.vascularmd.ArterialTree.ArterialTree.compute_cross_sections',
        'third_party.vascularmd.ArterialTree.ArterialTree.mesh_surface']
    from . import nn_support_gate
    names += ['vascular_processing.nn_support_gate.'+name for name in ('freeze_pilot','calibrate','freeze_predictions','calibrate_gate') if hasattr(nn_support_gate,name)]
    with ExitStack() as stack:
        mocks=[stack.enter_context(patch(name,side_effect=AssertionError('FORBIDDEN_RECOMPUTATION: '+name))) for name in names]
        yield names
        assert all(mock.call_count==0 for mock in mocks),'A forbidden production function was invoked'


@dataclass
class ReviewCache:
    directory: Path
    manifest: dict
    provenance: dict
    snapshot: dict
    branches: dict
    branch_graph: nx.DiGraph
    node_graph: nx.DiGraph
    coordinates: dict
    components: dict
    branch_component: dict
    registrations: list
    support: dict | None
    point_index: dict
    points: dict
    warnings: list

    @property
    def root(self):return next(n for n in self.node_graph if self.node_graph.in_degree(n)==0)

    def verify_unchanged(self):
        changed=[p for p,h in self.snapshot.items() if not Path(p).is_file() or sha256(p)!=h]
        if changed:raise ReviewSourceMismatch('Read-only invariant failed: '+str(changed))


def load_cache(directory):
    """Check hashes before creating any review output; missing optional votes warn."""
    import pyvista as pv
    directory=Path(directory).resolve();manifest_path=directory/'roi_manifest.json'
    try:
        manifest=json.loads(manifest_path.read_text())
        contract_path=directory.parent/'input_provenance.json'
        contract=json.loads(contract_path.read_text())
        source=manifest['source'];raw=Path(source['raw_source']);color=Path(source['color_source'])
        if manifest['side']!='RMCA' or contract['subject']!='BG001' or manifest['major_artery']['type_code']!=4:
            raise ReviewSourceMismatch('This exporter is locked to BG001 RMCA TYPE 4')
        for key,path,expected in [('raw',raw,RAW_SHA),('color',color,COLOR_SHA)]:
            if source[key+'_sha256']!=expected or contract[key+'_sha256']!=expected or sha256(path)!=expected:
                raise ReviewSourceMismatch(key+' source / manifest / input provenance disagree')
        if contract['major_arteries']['RMCA']!=4:raise ReviewSourceMismatch('Major artery mapping changed')
        snapshot={str(p):sha256(p) for p in [manifest_path,contract_path,raw,color]}
        warnings=[]
        for relative,digest in manifest['artifact_hashes'].items():
            path=(directory/relative).resolve()
            if not path.is_relative_to(directory):raise ReviewSourceMismatch('Artifact path escapes input directory')
            if relative=='point_support.npz' and not path.exists():
                warnings.append('DONOR_DETAIL_UNAVAILABLE: point_support.npz is absent');continue
            if not path.is_file() or sha256(path)!=digest:raise ReviewSourceMismatch('Frozen artifact hash: '+relative)
            snapshot[str(path)]=digest
        required=['point_predictions.csv','branch_predictions.csv','labeled_tree.vtp','node_mapping.csv',
                  'branch_mapping.csv','donor_registrations.json','donor_registrations.csv']
        if any(name not in manifest['artifact_hashes'] for name in required):
            raise ReviewSourceMismatch('Required prediction artifact is not bound by the ROI manifest')
        for component in manifest['roi_components']:
            for kind in ['strict','modelable']:
                info=component[kind];p=Path(info['path']).resolve()
                if sha256(p)!=info['sha256'] or str(p) not in snapshot:raise ReviewSourceMismatch('ROI hash / artifact binding: '+str(p))
            qc=Path(component['VascularMD_report']).resolve()
            if str(qc) not in snapshot:raise ReviewSourceMismatch('VascularMD QC not bound to manifest')
        for p in [ROOT/'s1-2_swc_roi_generate_human.py',ROOT/'s1-3_swc_roi_generate_MeVO.py']:
            snapshot[str(p)]=sha256(p)
        if snapshot[str(ROOT/'s1-2_swc_roi_generate_human.py')]!=HUMAN_SHA:raise ReviewSourceMismatch('s1-2 hash changed')
    except (KeyError,OSError,ValueError) as exc:
        if isinstance(exc,ReviewSourceMismatch):raise
        raise ReviewSourceMismatch(str(exc)) from exc
    mesh=pv.read(directory/'labeled_tree.vtp')
    if mesh.n_verts or mesh.n_faces or not np.all(mesh.lines.reshape(-1,3)[:,0]==2):
        raise ReviewSourceMismatch('Expected cached directed two-node SWC line cells')
    ids=mesh.point_data['original_swc_id'].astype(int)
    coords={int(n):np.r_[mesh.points[i],mesh.point_data['radius'][i]] for i,n in enumerate(ids)}
    raw_rows=np.loadtxt(raw);raw_lookup={int(r[0]):r for r in raw_rows}
    for n,xyzr in coords.items():
        if not np.array_equal(xyzr,raw_lookup[n][2:6]):raise ReviewSourceMismatch('Labeled VTP geometry differs from raw SWC')
    graph=nx.DiGraph();graph.add_nodes_from(coords)
    branch_edges=defaultdict(list)
    for cell,(a,b) in enumerate(mesh.lines.reshape(-1,3)[:,1:]):
        u,v=int(ids[a]),int(ids[b]);key=str(mesh.cell_data['branch_id'][cell])
        if int(raw_lookup[v][6])!=u:raise ReviewSourceMismatch('VTP parent direction differs from raw SWC')
        graph.add_edge(u,v,weight=float(np.linalg.norm(coords[v][:3]-coords[u][:3])))
        branch_edges[key].append((u,v))
    branches={}
    for row in read_csv(directory/'branch_predictions.csv'):
        key=row['branch_id'];edges=branch_edges[key];next_node=dict(edges)
        first=set(next_node)-{b for a,b in edges}
        if len(first)!=1:raise ReviewSourceMismatch('Cached branch is not an ordered polyline: '+key)
        path=[first.pop()]
        while path[-1] in next_node:path.append(next_node[path[-1]])
        if len(path)!=len(edges)+1:raise ReviewSourceMismatch('Cached branch path is incomplete')
        branch=dict(branch_id=key,parent_branch=row.get('parent_branch') or None,label=row['label'],
            p_M1=number(row.get('mean_p_M1')),p_MeVO=number(row.get('mean_p_MeVO')),
            known_fraction=number(row.get('known_fraction')),confidence=number(row.get('confidence')),
            length_mm=number(row.get('length_mm')),valid_donor_count=number(row.get('valid_donor_count_min')),
            support_distance=number(row.get('median_support_distance')),node_ids=path,
            radius_median=float(np.median([coords[n][3] for n in path])))
        for field in ['p_M1','p_MeVO','known_fraction','confidence','length_mm','valid_donor_count','support_distance']:
            if branch[field] is None:warnings.append(f'FIELD_NOT_AVAILABLE: {key}.{field}')
        branches[key]=branch
    topology=nx.DiGraph();topology.add_nodes_from(branches)
    for key,branch in branches.items():
        parent=branch['parent_branch']
        if parent in branches:
            topology.add_edge(parent,key)
            if branches[parent]['node_ids'][-1]!=branch['node_ids'][0]:raise ReviewSourceMismatch('Cached parent branch does not share junction')
    if not nx.is_arborescence(graph) or not nx.is_arborescence(topology):raise ReviewSourceMismatch('Cached RMCA is not one directed tree')
    components={c['component']:c for c in manifest['roi_components']};membership={}
    for row in read_csv(directory/'branch_mapping.csv'):
        if row['export_kind']=='strict':
            key=row['original_branch_id'];component=int(row['component'])
            if key in membership and membership[key]!=component:raise ReviewSourceMismatch('Branch belongs to two strict components')
            membership[key]=component
    if set(membership)!={k for k,b in branches.items() if b['label']=='MeVO'}:raise ReviewSourceMismatch('Frozen strict ROI membership differs from saved MeVO labels')
    support=None;point_index={}
    if (directory/'point_support.npz').exists():
        with np.load(directory/'point_support.npz',allow_pickle=False) as data:support={k:data[k] for k in data.files}
        point_index={int(n):i for i,n in enumerate(support['original_swc_id'])}
    points={int(r['original_swc_id']):r for r in read_csv(directory/'point_predictions.csv')}
    provenance=dict(subject='BG001',side='RMCA',major_tree_TYPE=4,mapping_status='DERIVED_SOURCE_SUPPORTED',
        scope='Only the locked BG001 standardized and ColorCoded source hashes',source=source,
        production_mapping_provenance=manifest['major_artery'],branch_predictions_sha256=sha256(directory/'branch_predictions.csv'),
        roi_manifest_sha256=sha256(manifest_path),branch_source_binding='ROI manifest artifact_hashes + source hashes + input_provenance.json',
        coordinates='Original SWC world XYZ in mm; no flips',direction='Native parent-to-child, not measured physiological flow')
    return ReviewCache(directory,manifest,provenance,snapshot,branches,topology,graph,coords,components,membership,
        json.loads((directory/'donor_registrations.json').read_text()),support,point_index,points,warnings)


def component_roots(cache,component):
    keys={k for k,c in cache.branch_component.items() if c==component}
    return sorted(k for k in keys if cache.branches[k]['parent_branch'] not in keys)


def detect_events(branches,components,branch_component,terminals):
    """Detect transitions without changing a single cached label or ROI member."""
    events=[]
    def add(kind,up,down,component,node):
        priority='CRITICAL' if kind=='MEVO_TO_M1_REVERSE' else 'MEDIUM' if kind=='DISTAL_MEVO_TO_TERMINAL' else 'HIGH'
        events.append(dict(boundary_id=stable_id(up,down,kind),boundary_type=kind,priority=priority,
            roi_component=component,boundary_node_id=node,upstream_branch_id=up,downstream_branch_id=down,review_status='UNREVIEWED'))
    for key in sorted(branches):
        b=branches[key];parent=b['parent_branch'];up=branches.get(parent)
        kind=TRANSITIONS.get((up['label'],b['label'])) if up else None
        if kind:add(kind,parent,key,branch_component.get(key,branch_component.get(parent)),b['node_ids'][0])
        if b['label']=='MeVO' and b['node_ids'][-1] in terminals:
            add('DISTAL_MEVO_TO_TERMINAL',key,None,branch_component[key],b['node_ids'][-1])
        component=branch_component.get(key)
        if component is not None and branch_component.get(parent)!=component:
            add('MEVO_COMPONENT_ROOT',parent,key,component,b['node_ids'][0])
            # A detached event means no immediate upstream M1 branch. It does
            # not assert a broken anatomical vessel or an incorrect candidate.
            if up is None or up['label']!='M1':add('DETACHED_MEVO_COMPONENT',parent,key,component,b['node_ids'][0])
    return sorted(events,key=lambda e:({'CRITICAL':0,'HIGH':1,'MEDIUM':2}[e['priority']],e['roi_component'] or 0,e['boundary_id']))


def local_keys(cache,event):
    """Boundary to two upstream branches and two downstream branch levels."""
    keys=set();up=event['upstream_branch_id']
    for _ in range(2):
        if up not in cache.branches:break
        keys.add(up);up=cache.branches[up]['parent_branch']
    node=event['boundary_node_id']
    first=[k for k,b in cache.branches.items() if b['node_ids'][0]==node]
    for key in first:
        keys.add(key);keys.update(cache.branch_graph.successors(key))
    return sorted(keys)


def donor_side(cache,branch_id,side):
    if branch_id not in cache.branches:return dict(node_id=None,reason='NO_BRANCH_ON_THIS_SIDE',donors=[])
    path=cache.branches[branch_id]['node_ids'];node=path[-2] if side=='upstream' else path[1]
    result=dict(node_id=node,sampling='nearest existing distinct SWC sample before/after shared junction; no new query',
        radius_mm=float(cache.coordinates[node][3]),agreement=number(cache.points[node].get('agreement')),donors=[])
    fields=['donor_ids','per_donor_binary_vote','per_donor_native_vote','per_donor_distance_mm','per_donor_normalized_distance']
    if cache.support is None or any(k not in cache.support for k in fields) or node not in cache.point_index:
        result['warning']='DONOR_DETAIL_UNAVAILABLE';return result
    registrations={r['source_case']:r for r in cache.registrations if 'rank' in r};i=cache.point_index[node]
    for j,case in enumerate(cache.support['donor_ids']):
        r=registrations.get(str(case),{})
        result['donors'].append(dict(donor_case=str(case),vote={1:'M1',2:'MeVO'}.get(int(cache.support['per_donor_binary_vote'][j,i]),'UNKNOWN'),
            native_vote=int(cache.support['per_donor_native_vote'][j,i]),nearest_distance=number(cache.support['per_donor_distance_mm'][j,i]),
            normalized_distance=number(cache.support['per_donor_normalized_distance'][j,i]),
            registration_fitness=number(r.get('fitness')),registration_RMSE=number(r.get('inlier_rmse'))))
    return result


def build_evidence(cache):
    terminals={n for n in cache.node_graph if cache.node_graph.out_degree(n)==0}
    events=detect_events(cache.branches,cache.components,cache.branch_component,terminals)
    root_distance=nx.single_source_dijkstra_path_length(cache.node_graph,cache.root,weight='weight')
    bif=[n for n in cache.node_graph if cache.node_graph.out_degree(n)>1]
    bif_distance=nx.multi_source_dijkstra_path_length(cache.node_graph.to_undirected(),bif,weight='weight') if bif else {}
    details={};local_topologies={}
    for event in events:
        node=event['boundary_node_id'];event.update(zip(['world_x','world_y','world_z'],map(float,cache.coordinates[node][:3])))
        directions={};sides={}
        for side in ['upstream','downstream']:
            key=event[side+'_branch_id'];branch=cache.branches.get(key)
            for field,source in [('label','label'),('p_M1','p_M1'),('p_MeVO','p_MeVO'),('known_fraction','known_fraction'),
                ('confidence','confidence'),('median_support_distance','support_distance'),('valid_donor_count','valid_donor_count'),
                ('length_mm','length_mm'),('radius_median','radius_median')]:event[side+'_'+field]=branch.get(source) if branch else None
            sides[side]=donor_side(cache,key,side)
            if branch:
                nodes=branch['node_ids'];a,b=nodes[-2:] if side=='upstream' else nodes[:2]
                vector=cache.coordinates[b][:3]-cache.coordinates[a][:3];norm=float(np.linalg.norm(vector))
                directions[side]=(vector/norm).tolist() if norm>0 else None
            else:directions[side]=None
        r0,r1=event['upstream_radius_median'],event['downstream_radius_median']
        event.update(radius_ratio=r1/r0 if r0 and r1 is not None else None,
            distance_from_RMCA_root_mm=float(root_distance[node]),distance_to_nearest_bifurcation_mm=number(bif_distance.get(node)),
            vascularmd_component_status=cache.components[event['roi_component']]['VascularMD_status'])
        keys=local_keys(cache,event);local_nodes=sorted({n for k in keys for n in cache.branches[k]['node_ids']})
        angle=None
        if all(directions.values()):angle=float(np.degrees(np.arccos(np.clip(np.dot(*directions.values()),-1,1))))
        up=cache.branches.get(event['upstream_branch_id']);down=cache.branches.get(event['downstream_branch_id'])
        parent=up['parent_branch'] if up else None
        daughters=list(cache.branch_graph.successors(down['branch_id'])) if down else []
        grandchildren=[c for d in daughters for c in cache.branch_graph.successors(d)]
        detail=dict(**event,upstream_parent_branch=parent,grandparent_branch=cache.branches.get(parent,{}).get('parent_branch'),
            downstream_daughter_branches=daughters,grandchildren_branches=grandchildren,
            branch_direction_vectors=directions,boundary_angle_degrees=angle,
            radius_before_mm=sides['upstream'].get('radius_mm'),radius_after_mm=sides['downstream'].get('radius_mm'),
            branch_length_before_mm=event['upstream_length_mm'],branch_length_after_mm=event['downstream_length_mm'],
            donor_vote_breakdown=sides,top_donors=[{k:r.get(k) for k in ['source_case','rank','fitness','inlier_rmse']} for r in cache.registrations if 'rank' in r],
            local_topology_edge_list=[list(e) for e in cache.branch_graph.edges if e[0] in keys and e[1] in keys],
            local_node_list=[dict(original_swc_id=n,xyz_mm=cache.coordinates[n][:3].tolist(),radius_mm=float(cache.coordinates[n][3])) for n in local_nodes],
            distance_definition='Arc length along native RMCA tree; nearest bifurcation is undirected geodesic, includes boundary itself',
            radius_ratio_definition='downstream original sample-median radius / upstream original sample-median radius',
            flow_direction_caveat='Parent-to-child direction only; physiological flow was not measured')
        details[event['boundary_id']]=detail
        local_topologies[event['boundary_id']]=dict(boundary_id=event['boundary_id'],review_status='UNREVIEWED',
            upstream_depth_max=2,downstream_depth_max=2,
            branches=[dict(cache.branches[k],is_unknown=cache.branches[k]['label']=='UNKNOWN') for k in keys],
            edges=detail['local_topology_edge_list'],nodes=detail['local_node_list'])
    return events,details,local_topologies


def separation_pairs(cache):
    """Explain saved components using their existing branch-parent paths."""
    result=[];undirected=cache.branch_graph.to_undirected()
    for a,b in itertools.combinations(sorted(cache.components),2):
        aa=[k for k,c in cache.branch_component.items() if c==a];bb=[k for k,c in cache.branch_component.items() if c==b]
        paths=[nx.shortest_path(undirected,u,v) for u in aa for v in bb if nx.has_path(undirected,u,v)]
        if not paths:reason='NO_DIRECT_CONNECTION';path=[]
        else:
            path=min(paths,key=lambda p:(len(p),p))
            gap=[k for k in path if k not in aa and k not in bb]
            if any(cache.branches[k]['label']=='UNKNOWN' for k in gap):reason='UNKNOWN_GAP'
            elif {cache.branches[k]['parent_branch'] for k in component_roots(cache,a)}!={cache.branches[k]['parent_branch'] for k in component_roots(cache,b)}:
                reason='DIFFERENT_PROXIMAL_PARENT'
            elif gap:reason='TOPOLOGICALLY_SEPARATE_DAUGHTERS'
            else:reason='OTHER'
        result.append(dict(component_a=a,component_b=b,reason=reason,connecting_branch_path=path,
            intervening_labels=[cache.branches[k]['label'] for k in path],modification_performed=False))
    return result


def component_summary(cache,events,pairs):
    rows=[]
    for c,component in sorted(cache.components.items()):
        keys=sorted(k for k,v in cache.branch_component.items() if v==c);roots=component_roots(cache,c)
        lengths=[cache.branches[k]['length_mm'] for k in keys]
        def weighted(field):
            if any(x is None for x in lengths) or any(cache.branches[k][field] is None for k in keys):return None
            return float(np.average([cache.branches[k][field] for k in keys],weights=lengths))
        proximal=[e['boundary_id'] for e in events if e['roi_component']==c and e['boundary_type']=='MEVO_COMPONENT_ROOT']
        reasons=sorted({p['reason'] for p in pairs if c in (p['component_a'],p['component_b'])})
        rows.append(dict(component_id=f'roi_part{c:02d}',node_count=component['strict']['node_count'],
            branch_count=component['strict']['original_branch_count'],total_length_mm=sum(lengths) if all(x is not None for x in lengths) else None,
            proximal_branch_id=';'.join(roots),proximal_boundary_id=';'.join(proximal),
            proximal_parent_label=';'.join(sorted({cache.branches.get(cache.branches[k]['parent_branch'],{}).get('label','OUTSIDE_RMCA') for k in roots})),
            distal_terminal_count=sum(e['roi_component']==c and e['boundary_type']=='DISTAL_MEVO_TO_TERMINAL' for e in events),
            distal_unknown_boundary_count=sum(e['roi_component']==c and e['boundary_type']=='DISTAL_MEVO_TO_UNKNOWN' for e in events),
            MeVO_mean_probability=weighted('p_MeVO'),known_fraction=weighted('known_fraction'),
            VascularMD_status=component['VascularMD_status'],VascularMD_context_branch=';'.join(component['context_branch_ids']),
            strict_swc_path=f'evidence/ROI/roi_part{c:02d}.swc',modelable_swc_path=f'evidence/ROI/roi_part{c:02d}_modelable.swc',
            surface_path=component.get('model_outputs',{}).get('surface_vtk'),STL_path=component.get('model_outputs',{}).get('surface_stl'),
            surface_and_STL_included=False,separation_reason=';'.join(reasons),semantic_candidate_status='ANATOMICAL_TRANSFER_CANDIDATE',review_status='UNREVIEWED'))
    return rows
