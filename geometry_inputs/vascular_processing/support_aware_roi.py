"""Frozen-data subtree selection for the independent s1-7 development stage."""
from __future__ import annotations
import copy
import csv
import itertools
import json
from pathlib import Path
import shutil
import networkx as nx
import numpy as np
from .topbrain_qc import sha256
from .swc_export import read_source
from .manufacturing_roi import compensate_radius,export_manufacturing
from .sacrificial_fixture_review import save_json
from . import roi_support_metrics as metrics

ROOT=Path(__file__).resolve().parents[1]
RMCA=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA'
DEFAULT_OUTPUT=RMCA/'support_aware_roi_60x50'


def table(path,rows,fields=None):
    """Replace only a table in this new stage; union columns for partial results."""
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(dict.fromkeys(k for r in rows for k in r)) or ['status']
    temporary=path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    temporary.replace(path)


def read_csv(path):
    with Path(path).open() as f:return list(csv.DictReader(f))


def verify_snapshot(snapshot):
    changed=[p for p,h in snapshot.items() if not Path(p).is_file() or sha256(p)!=h]
    return dict(all_unchanged=not changed,changed=changed,checked_files=len(snapshot))


def snapshot_sources(rmca,out):
    files=list(rmca.rglob('*'))+list(ROOT.glob('s1-[2-6]*.py'))+list((ROOT/'vascular_processing').glob('*.py'))+list((ROOT/'config').glob('*.yaml'))+list((ROOT/'third_party/vascularmd').rglob('*.py'))+list((ROOT/'utils').rglob('*.py'))
    # All pre-existing output files are protected, including previous production
    # STLs, model companions, native fitted surfaces, reports and slicer archives.
    excluded={'support_aware_roi.py','roi_support_metrics.py','support_aware_port_design.py','support_aware_evaluation.py','support_aware_review.py','support_aware_roi_60x50.yaml'}
    paths=sorted({p.resolve() for p in files if p.is_file() and not p.resolve().is_relative_to(out.resolve()) and p.name not in excluded})
    return {str(p):sha256(p) for p in paths}


def discover(rmca=RMCA):
    rmca=Path(rmca).resolve();manifest=json.loads((rmca/'roi_manifest.json').read_text())
    if manifest['side']!='RMCA':raise ValueError('RMCA_REQUIRED')
    provenance=json.loads((rmca.parent/'input_provenance.json').read_text())
    if provenance['subject']!='BG001' or provenance['major_arteries']['RMCA']!=4:raise ValueError('BG001_RMCA_REQUIRED')
    raw=read_source(Path(manifest['source']['raw_source']))
    for kind in ('raw','color'):
        if sha256(Path(manifest['source'][kind+'_source']))!=manifest['source'][kind+'_sha256']:raise ValueError('FROZEN_SOURCE_HASH_MISMATCH')
    for rel,digest in manifest['artifact_hashes'].items():
        p=(rmca/rel).resolve()
        if not p.is_relative_to(rmca) or sha256(p)!=digest:raise ValueError('FROZEN_ARTIFACT_MISMATCH: '+rel)
    refined=rmca/'refined_roi';refinement=json.loads((refined/'refinement_summary.json').read_text())
    for rel,digest in json.loads((refined/'artifact_manifest.json').read_text())['files'].items():
        p=(refined/rel).resolve()
        if not p.is_relative_to(refined) or sha256(p)!=digest:raise ValueError('REFINEMENT_HASH_MISMATCH: '+rel)
    branches={r['branch_id']:r for r in read_csv(rmca/'branch_predictions.csv')}
    sources=[];excluded=[];edge_branch={}
    semantic={}
    for e in refinement['exports']:
        p=Path(e['path']);rows=read_csv(refined/'mappings'/(p.stem+'_nodes.csv'));ids={int(r['original_swc_id']) for r in rows}
        edges=read_csv(refined/'mappings'/(p.stem+'_edges.csv'))
        edge_branch.update({(int(r['original_parent_id']),int(r['original_child_id'])):r['original_branch_id'] for r in edges})
        if e['level']=='semantic_refined':semantic[e['derived_component'] if 'derived_component' in e else e['component']]=ids
        if e['level']!='diameter_075':continue
        if sha256(p)!=e['sha256']:raise ValueError('FROZEN_D075_HASH_MISMATCH')
        graph=raw.graph.subgraph(ids).copy()
        # No new classification: intersect the saved diameter result with the
        # frozen strict MeVO branch labels. UNKNOWN historical connector removed.
        rejected=[(a,b) for a,b in graph.edges if branches[edge_branch[a,b]]['label']!='MeVO']
        for a,b in rejected:excluded.append(dict(source=p.name,parent=a,child=b,branch=edge_branch[a,b],label=branches[edge_branch[a,b]]['label']))
        graph.remove_edges_from(rejected);graph.remove_nodes_from(list(nx.isolates(graph)))
        for j,group in enumerate(sorted(nx.weakly_connected_components(graph),key=min)):
            g=graph.subgraph(group).copy()
            if sum(g.out_degree(n)==0 for n in g)<2:continue
            sources.append(dict(id=p.stem+f'_{j+1}',path=str(p),sha256=e['sha256'],component=e['component'],graph=g))
    compact=json.loads((rmca/'compact_manufacturing_roi/compact_manifest.json').read_text())
    t=metrics.rigid_transform(json.loads((rmca/'compact_manufacturing_roi/print_transform.json').read_text())['transform_4x4'])
    return dict(raw=raw,sources=sources,branches=branches,edge_branch=edge_branch,excluded=excluded,
        transform=t,profile=compact['profile'],compact=compact,
        source=dict(subject='BG001',side='RMCA',refinement=str(refined),raw=str(raw.path),
            diameter_policy='EXISTING_DIAMETER_075_INTERSECT_FROZEN_MEVO_BRANCHES',semantic_recomputed=False))


def mevo_gate(graph,edge_branch,branches):
    return bool(len(graph.edges) and all(branches[edge_branch[e]]['label']=='MeVO' for e in graph.edges))


def generate_candidates(data,minimum=20):
    """Native-node cuts of 2/3 root-to-tip paths, deterministic and deduplicated.

    Enumerate proximal and distal native cuts; never interpolate or add a node.
    Sampling is fixed by graph order and physical terminal length, not randomness.
    """
    candidates=[];seen=set();t=data['transform']
    floor=data['profile']['effective']['radius_floor_mm'];blend=data['profile']['manufacturing']['radius_blend_length_mm']
    for source in data['sources']:
        g=source['graph'];root=next(n for n in g if g.in_degree(n)==0)
        paths=nx.single_source_shortest_path(g,root);leaves=sorted(n for n in g if g.out_degree(n)==0)
        for count in (3,2):
            for tips in itertools.combinations(leaves,count):
                stem=set.intersection(*(set(paths[n]) for n in tips));stem=sorted(stem,key=lambda n:len(paths[n]))
                for inlet in stem[:-1]:
                    tip_options=[]
                    for tip in tips:
                        p=paths[tip];last=max(j for j,n in enumerate(p) if n in stem or g.out_degree(n)>1)
                        terminal=p[last:];arc=np.r_[0.,np.cumsum([np.linalg.norm(g.nodes[b]['coords'][:3]-g.nodes[a]['coords'][:3]) for a,b in zip(terminal,terminal[1:])])]
                        valid=[n for n,s in zip(terminal,arc) if s>=3]
                        if not valid:break
                        # At most three deterministic native cut positions/leaf.
                        tip_options.append(sorted({valid[0],valid[len(valid)//2],valid[-1]}))
                    if len(tip_options)!=count:continue
                    for ends in itertools.product(*tip_options):
                        nodes=set.union(*(set(paths[n][paths[n].index(inlet):]) for n in ends));key=tuple(sorted(nodes))
                        if key in seen:continue
                        seen.add(key);sub=g.subgraph(nodes).copy();stats=metrics.topology_metrics(sub)
                        if not stats['topology_valid'] or not mevo_gate(sub,data['edge_branch'],data['branches']):continue
                        comp,_=compensate_radius(sub,floor,blend)
                        xyzr=np.array([comp.nodes[n]['coords'] for n in sorted(comp)])
                        xyz=xyzr[:,:3]@t[:3,:3].T+t[:3,3]
                        # Bounding envelope of compensated centerline spheres,
                        # deliberately labelled proxy; final gate uses STL.
                        dim=metrics.dimensions(np.vstack([xyz-xyzr[:,3,None],xyz+xyzr[:,3,None]]))
                        proxy=metrics.centerline_support(comp,t)
                        candidates.append(dict(source_id=source['id'],source_component=source['component'],native_node_ids=list(key),
                            root_branch=data['edge_branch'][next(iter(sub.out_edges(stats['root_node'])))],
                            mevo_only=True,scale_factor=1.,minimum_radius_mm=float(xyzr[:,3].min()),maximum_radius_mm=float(xyzr[:,3].max()),
                            **stats,**dim,**{k:v for k,v in proxy.items() if k!='segments'}))
    # If no size-valid candidate exists, honour the nearest-candidate fallback
    # before ranking support. All independent terms remain in the JSON/CSV.
    def rank(c):
        return (not c['topology_valid'],not c['mevo_only'],not c['size_target_met'],
            c['size_violation_mm'] if not c['size_target_met'] else 0,
            c['estimated_support_contact_mm2'],c['support_critical_branches'],c['support_near_bifurcation_count'],
            c['branch_count'],c['Z_span_mm'],c['source_component']!=4,c['root_node'],c['native_node_ids'])
    candidates.sort(key=rank)
    for j,c in enumerate(candidates,1):c.update(candidate_id=f'candidate_{j:04}',rank=j)
    if len(candidates)<minimum:raise ValueError(f'INSUFFICIENT_VALID_SUBTREES: {len(candidates)} < {minimum}')
    return candidates


def prepare(out=DEFAULT_OUTPUT):
    out=Path(out).resolve()
    if (out/'stage_status.json').exists():raise FileExistsError('Prepared output exists; use --stage model/evaluate/report: '+str(out))
    data=discover();out.mkdir(parents=True,exist_ok=True)
    for name in ('baseline','candidates','selected','casting_mold','tables','QC','bambu'):(out/name).mkdir(exist_ok=True)
    snapshot=json.loads((out/'protected_before.json').read_text()) if (out/'protected_before.json').exists() else snapshot_sources(RMCA,out)
    if not verify_snapshot(snapshot)['all_unchanged']:raise ValueError('PROTECTED_SOURCE_CHANGED')
    for p in (data['raw'].path,Path(json.loads((RMCA/'roi_manifest.json').read_text())['source']['color_source'])):snapshot[str(p.resolve())]=sha256(p)
    save_json(out/'protected_before.json',snapshot)
    save_json(out/'roi_casting_transform.json',dict(transform_4x4=data['transform'],scale_factor=1.,source=str(RMCA/'compact_manufacturing_roi/print_transform.json'),policy='SAME_FIXED_LEGACY_CASTING_FRAME_FOR_ALL_CANDIDATES'))
    candidates=generate_candidates(data)
    save_json(out/'candidates.json',candidates);table(out/'tables/anatomy_exclusions.csv',data['excluded'])
    table(out/'tables/roi_candidate_summary.csv',[{k:v for k,v in c.items() if k not in ('native_node_ids','bounds_mm')} for c in candidates])
    save_json(out/'source_discovery.json',dict(**data['source'],sources=[{k:v for k,v in s.items() if k!='graph'} for s in data['sources']],
        frozen_profile=data['profile'],excluded_non_mevo_edge_count=len(data['excluded']),candidate_count=len(candidates),
        ranking=['topology','MeVO','size','nearest size if no size-valid result','support contact','critical branches','bifurcation contacts','branch count','Z span','preferred part04','native IDs']))
    for c in candidates[:8]:
        directory=out/'candidates'/c['candidate_id'];directory.mkdir()
        graph=data['raw'].graph.subgraph(c['native_node_ids']).copy();comp,windows=compensate_radius(graph,data['profile']['effective']['radius_floor_mm'],data['profile']['manufacturing']['radius_blend_length_mm'])
        exports={}
        for name,g in [('centerline',graph),('compensated',comp)]:exports[name]=export_manufacturing(directory/(name+'.swc'),g,data['raw'],set(graph),data['edge_branch'],name=='compensated')
        c.update(exports=exports,compensation_windows=windows)
        save_json(directory/'candidate.json',c)
        table(directory/'centerline_support_segments.csv',metrics.centerline_support(comp,data['transform'])['segments'])
    save_json(out/'stage_status.json',dict(stage='PREPARED',candidate_count=len(candidates),top8=[c['candidate_id'] for c in candidates[:8]],source_protection=verify_snapshot(snapshot)))
    print('PREPARED',len(candidates),'candidates; Top8',[c['candidate_id'] for c in candidates[:8]],flush=True)


def model(out=DEFAULT_OUTPUT):
    import pyvista as pv
    import trimesh
    import yaml
    from tools.refine_bg001_rmca import model_one
    from .print_orientation import cap_native_ports
    from .sacrificial_fixture import mesh_qc
    out=Path(out).resolve();snapshot=json.loads((out/'protected_before.json').read_text())
    if not verify_snapshot(snapshot)['all_unchanged']:raise ValueError('PROTECTED_SOURCE_CHANGED')
    t=np.array(json.loads((out/'roi_casting_transform.json').read_text())['transform_4x4'])
    cfg=yaml.safe_load((ROOT/'config/abs_casting_mold_BG001_2mm_X1C_02.yaml').read_text());results=[]
    for directory in sorted((out/'candidates').iterdir()):
        c=json.loads((directory/'candidate.json').read_text())
        for export in c['exports'].values():
            if sha256(export['path'])!=export['sha256']:raise ValueError('PREPARED_CANDIDATE_CHANGED')
        if (directory/'model_result.json').exists():results.append(json.loads((directory/'model_result.json').read_text()));continue
        result=dict(candidate_id=c['candidate_id'],failures=[])
        try:
            native=model_one(directory/'compensated.swc',directory/'VascularMD');result['VascularMD']=native
            if native['status']!='success' or not native['topology_preserved']:raise ValueError('VASCULARMD_FAILED_OR_TOPOLOGY_CHANGED')
            closed,capqc=cap_native_ports(pv.read(native['outputs']['surface_vtk']));result['cap_qc']=capqc
            closed.save(directory/'surface_native.stl');closed.transform(t,inplace=True)
            closed.save(directory/'vascular_core.vtp');closed.save(directory/'vascular_core.stl')
            mesh=trimesh.load_mesh(directory/'vascular_core.stl',process=True);result['mesh']=mesh_qc(mesh,cfg)
            graph=read_source(Path(native['outputs']['swc'])).graph
            fitted_radii=np.array([graph.nodes[n]['coords'][3] for n in graph])
            result['radius_qc']=dict(minimum_fitted_radius_mm=float(fitted_radii.min()),approved_floor_mm=.6,
                passed=bool(fitted_radii.min()>=.6-.01),post_fit_clamp=False)
            result['dimensions']=metrics.dimensions(mesh.vertices)
            result['support'],colors=metrics.surface_support(mesh,graph,t)
            np.savez_compressed(directory/'surface_support.npz',risk_class=colors)
            result['pdms']=metrics.vascular_ligament(mesh,graph,t,cfg)
            if not result['mesh']['passed']:result['failures'].append('SURFACE_MESH_QC_FAILED')
            if not result['radius_qc']['passed']:result['failures'].append('FITTED_RADIUS_BELOW_APPROVED_FLOOR')
            if not result['pdms']['passed']:result['failures'].append('PDMS_LIGAMENT_TOO_THIN')
            if not result['dimensions']['size_target_met']:result['failures'].append('ROI_SIZE_TARGET_NOT_MET')
            table(directory/'pdms_ligament.csv',result['pdms']['all_pairs'])
        except Exception as exc:
            import traceback
            result['failures'].append(type(exc).__name__+': '+str(exc));(directory/'error.log').write_text(traceback.format_exc())
        result['status']='MESH_READY' if not result['failures'] else 'NEEDS_ADJUSTMENT'
        save_json(directory/'model_result.json',result);results.append(result)
        print(c['candidate_id'],result['status'],result.get('dimensions'),result['failures'],flush=True)
    save_json(out/'model_results.json',results)
    rows=[]
    for c in json.loads((out/'candidates.json').read_text()):
        r=next((r for r in results if r['candidate_id']==c['candidate_id']),{})
        rows.append({**{k:v for k,v in c.items() if k not in ('native_node_ids','bounds_mm')},**r.get('dimensions',{}),**r.get('support',{}),
            'minimum_pdms_gap_mm':r.get('pdms',{}).get('minimum_ligament_mm'),'status':r.get('status','SURROGATE_ONLY'),'failures':r.get('failures',[])})
    table(out/'tables/roi_candidate_summary.csv',rows)
    save_json(out/'protected_after_model.json',verify_snapshot(snapshot))
