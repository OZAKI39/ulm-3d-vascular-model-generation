"""Existing BraVa topo_graph branches, exact version mapping and length voting."""
from pathlib import Path
import numpy as np
import networkx as nx
from scipy.spatial import cKDTree
from .swc_export import read_source
from .mevo_graph import load_exact_graph
from .topbrain_qc import sha256

LABEL_NAMES={0:'UNKNOWN',1:'M1',2:'MeVO'}


def color_source_mapping(exact,colored_path):
    """Map the re-rooted/re-numbered four-significant-digit companion to raw.

    Only accept a bijection within decimal serialization precision and identical
    undirected edges. This is version correspondence, not anatomical inference.
    Original coordinates, radii and parent directions are never changed.
    """
    colored=read_source(Path(colored_path),allow_forest=True)
    raw_ids=list(exact.graph);color_ids=list(colored.graph)
    a=np.array([exact.graph.nodes[n]['coords'] for n in raw_ids])
    b=np.array([colored.graph.nodes[n]['coords'] for n in color_ids])
    distance,index=cKDTree(b[:,:3]).query(a[:,:3])
    tolerance=.5*10.**(np.floor(np.log10(np.maximum(abs(a),1e-12)))-3)+1e-10
    if len(a)!=len(b) or len(set(index))!=len(a) or np.any(abs(a-b[index])>tolerance):
        raise ValueError('SOURCE_MAPPING_FAILED: not a unique correspondence within ColorCoded decimal precision')
    mapping={int(n):int(color_ids[j]) for n,j in zip(raw_ids,index)}
    if {frozenset((mapping[a],mapping[b])) for a,b in exact.graph.edges}!={frozenset(e) for e in colored.graph.edges}:
        raise ValueError('SOURCE_MAPPING_FAILED: original undirected topology differs')
    labels={n:colored.graph.nodes[c]['swc_type'] for n,c in mapping.items()}
    return labels,mapping,dict(raw_source=str(exact.source.path),raw_sha256=sha256(exact.source.path),
        color_source=str(colored.path),color_sha256=sha256(colored.path),mapped_nodes=len(mapping),
        coordinate_max_difference_mm=float(np.max(abs(a[:,:3]-b[index,:3]))),
        radius_max_difference_mm=float(np.max(abs(a[:,3]-b[index,3]))),undirected_topology_identical=True,
        policy='Bijective nearest coordinates within four-significant-digit rounding; radius and all undirected edges checked. Raw directed tree retained.')


def topology_branches(exact):
    branches={}
    for a,b,data in exact.topo.edges(data=True):
        nodes=[int(exact.topo.nodes[a]['full_id']),*map(int,data['full_id']),int(exact.topo.nodes[b]['full_id'])]
        key=f'{a}_{b}'
        branches[key]=dict(branch_id=key,topo_start=int(a),topo_end=int(b),node_ids=nodes)
    end_index={r['topo_end']:key for key,r in branches.items()}
    for row in branches.values():row['parent_branch']=end_index.get(row['topo_start'])
    return branches


def select_major_branches(exact,major_labels,type_code):
    branches=topology_branches(exact);selected={}
    for key,row in branches.items():
        codes=[major_labels[n] for n in row['node_ids'][1:]]
        # All non-junction samples must belong to the provided major label;
        # shared endpoints can carry the neighboring artery's annotation.
        internal=[major_labels[n] for n in row['node_ids'][1:-1]]
        if (internal and all(c==type_code for c in internal)) or (not internal and codes[-1]==type_code):
            selected[key]=row
        elif type_code in internal:
            raise ValueError(f'MAJOR_LABEL_CROSSES_TOPO_BRANCH: {key}; cannot split the original branch')
    if not selected:raise ValueError('MAJOR_ARTERY_UNAVAILABLE')
    graph=nx.DiGraph()
    for row in selected.values():
        nodes=row['node_ids'];graph.add_nodes_from((n,dict(exact.graph.nodes[n])) for n in nodes)
        graph.add_edges_from(zip(nodes[:-1],nodes[1:]))
    if not nx.is_arborescence(graph):raise ValueError('MAJOR_ARTERY_NOT_ONE_ORIGINAL_TREE')
    return graph,selected


def arc_weights(points):
    length=np.linalg.norm(np.diff(points,axis=0),axis=1)
    weights=np.r_[length,0]+np.r_[0,length]
    return weights/2


def aggregate_branches(exact,branches,node_ids,features,labels):
    index={int(n):i for i,n in enumerate(node_ids)};rows=[]
    for key,branch in branches.items():
        ids=branch['node_ids'];ii=np.array([index[n] for n in ids])
        coords=np.array([exact.graph.nodes[n]['coords'] for n in ids]);weight=arc_weights(coords[:,:3])
        total=float(weight.sum())
        if total<=0:raise ValueError(f'Degenerate branch length: {key}')
        w=weight/total
        p1=float(w@features['p_M1'][ii]);p2=float(w@features['p_MeVO'][ii])
        known=float(w@(labels[ii]!=0))
        label=0 if known<.6 else 2 if p2>=.6 else 1 if p1>=.6 else 0
        rows.append(dict(branch_id=key,parent_branch=branch['parent_branch'],length_mm=total,
            mean_p_M1=p1,mean_p_MeVO=p2,known_fraction=known,label=LABEL_NAMES[label],label_code=label,
            confidence=max(p1,p2),median_support_distance=float(np.median(features['median_normalized_distance'][ii])),
            valid_donor_count_min=int(features['valid_donor_count'][ii].min()),
            valid_donor_count_median=float(np.median(features['valid_donor_count'][ii])),node_count=len(ids)))
    return rows


def topology_warnings(branches,rows):
    by_id={r['branch_id']:r for r in rows};warnings=[]
    for row in rows:
        ancestor=row['parent_branch'];gap=False
        while ancestor in by_id:
            up=by_id[ancestor]
            if up['label']=='MeVO':
                if row['label']=='M1':
                    warnings.append(dict(status='HIGH_CONFIDENCE_TOPOLOGY_REVERSE' if row['confidence']>=.8 else 'TOPOLOGY_REVERSE',
                        upstream_branch=ancestor,downstream_branch=row['branch_id'],unknown_gap=gap))
                elif row['label']=='MeVO' and gap:
                    warnings.append(dict(status='ROI_INTERRUPTED_BY_UNKNOWN',upstream_branch=ancestor,downstream_branch=row['branch_id']))
                break
            if up['label']=='UNKNOWN':gap=True
            ancestor=up['parent_branch']
    return warnings
