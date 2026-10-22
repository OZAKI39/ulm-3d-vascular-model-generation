"""Exact original-branch subsets, independent components and native modeling."""
from pathlib import Path
import csv
import json
import io
import tempfile
import networkx as nx
import numpy as np
from .swc_export import write_swc,read_source,validate_tree
from .mevo_graph import load_exact_graph
from .topbrain_qc import sha256,write_json
from .pipeline import process_file,Options
from .vascularmd_adapter import COMMIT


def write_csv(path,rows):
    stream=io.StringIO(newline='')
    if rows:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    data=stream.getvalue().encode('utf-8');path=Path(path)
    if path.exists():
        if path.read_bytes()!=data:raise ValueError('Existing CSV belongs to different results: '+str(path))
    else:
        with path.open('xb') as handle:handle.write(data)


def branch_union(exact,branches,keys):
    graph=nx.DiGraph()
    for key in keys:
        nodes=branches[key]['node_ids']
        graph.add_nodes_from((n,dict(exact.graph.nodes[n])) for n in nodes)
        graph.add_edges_from(zip(nodes[:-1],nodes[1:]))
    return graph


def strict_components(exact,branches,rows):
    keys=[r['branch_id'] for r in rows if r['label']=='MeVO']
    graph=branch_union(exact,branches,keys)
    return [graph.subgraph(c).copy() for c in sorted(nx.weakly_connected_components(graph),key=min)] if graph else []


def export_component(path,graph,exact,branches,strict_nodes,context_keys=()):
    headers=[
        'TopBrain-informed transferred candidate labels; NOT manual anatomical ground truth',
        f'Source: {exact.source.path}',f'Source SHA256: {sha256(exact.source.path)}',
        'Original XYZ, radius, TYPE and retained directed edges preserved; roots only reset to -1',
        'ANATOMICAL_TRANSFER_CANDIDATE; MANUFACTURING_NOT_VALIDATED',
        'Context branches: '+','.join(context_keys)]
    # Resume safely without overwriting a previous exact subset.
    with tempfile.TemporaryDirectory() as temporary:
        candidate=Path(temporary)/path.name
        mapping=write_swc(candidate,graph,exact.source.path,COMMIT,'exact-original-branch-subset',header_lines=headers)
        data=candidate.read_bytes()
        if path.exists():
            if path.read_bytes()!=data:raise ValueError('Existing ROI does not match current exact subset')
        else:
            with path.open('xb') as handle:handle.write(data)
    saved=read_source(path);validate_tree(saved.graph)
    if set(saved.graph.edges)!={(mapping[a],mapping[b]) for a,b in graph.edges}:raise ValueError('ROI connectivity changed')
    node_rows=[]
    for original,new in mapping.items():
        np.testing.assert_array_equal(saved.graph.nodes[new]['coords'],exact.graph.nodes[original]['coords'])
        if saved.graph.nodes[new]['swc_type']!=exact.graph.nodes[original]['swc_type']:raise ValueError('ROI TYPE changed')
        node_rows.append(dict(roi_node_id=new,original_swc_id=original,
            original_parent_id=next(exact.graph.predecessors(original),-1),
            roi_parent_id=next(saved.graph.predecessors(new),-1),
            is_context=original not in strict_nodes,is_anatomical_roi=original in strict_nodes))
    inverse={v:k for k,v in mapping.items()};roi_exact=load_exact_graph(path);branch_rows=[]
    for a,b,data in roi_exact.topo.edges(data=True):
        ids=[int(roi_exact.topo.nodes[a]['full_id']),*map(int,data['full_id']),int(roi_exact.topo.nodes[b]['full_id'])]
        old=[inverse[n] for n in ids]
        originals=sorted({f'{u}_{v}' for s,t in zip(old[:-1],old[1:]) for u,v in [exact.edge_locations[(s,t)]]})
        for key in originals:
            branch_rows.append(dict(roi_branch_id=f'{a}_{b}',original_branch_id=key,
                is_context=key in context_keys,is_anatomical_roi=key not in context_keys))
    write_csv(path.with_name(path.stem+'_node_mapping.csv'),node_rows)
    write_csv(path.with_name(path.stem+'_branch_mapping.csv'),branch_rows)
    root=validate_tree(graph)
    return dict(path=str(path.resolve()),sha256=sha256(path),node_count=len(graph),branch_count=len(roi_exact.topo.edges),
        original_branch_count=len({r['original_branch_id'] for r in branch_rows}),root_original_id=root,
        root_out_degree=graph.out_degree(root),geometry_exact_preserved=True,types_preserved=True,
        directed_edges_preserved=True,node_mapping=node_rows,branch_mapping=branch_rows)


def extract_and_model(exact,branches,rows,output,*,model=True):
    output=Path(output);folder=output/'ROI';folder.mkdir(exist_ok=True)
    components=strict_components(exact,branches,rows);by_key={r['branch_id']:r for r in rows}
    all_nodes=[];all_branches=[];reports=[]
    for i,strict in enumerate(components,1):
        tag=f'roi_part{i:02d}'
        exported=export_component(folder/(tag+'.swc'),strict,exact,branches,set(strict))
        root=validate_tree(strict)
        starting=[k for k,b in branches.items() if b['node_ids'][0]==root and strict.has_edge(*b['node_ids'][:2])]
        parents={branches[k]['parent_branch'] for k in starting}
        context=[k for k in parents if k in by_key and by_key[k]['label']=='M1']
        modelable=strict.copy()
        if context:
            addition=branch_union(exact,branches,context);modelable.update(addition)
        modeled_input=export_component(folder/(tag+'_modelable.swc'),modelable,exact,branches,set(strict),context)
        item=dict(component=i,strict=exported,modelable=modeled_input,context_branch_ids=context,
            status='BRAVA_MEVO_ROI_EXTRACTED',VascularMD_status='NOT_RUN',manufacturing_status='MANUFACTURING_NOT_VALIDATED')
        for kind,info in [('strict',exported),('modelable',modeled_input)]:
            all_nodes.extend(dict(component=i,export_kind=kind,**r) for r in info.pop('node_mapping'))
            all_branches.extend(dict(component=i,export_kind=kind,**r) for r in info.pop('branch_mapping'))
        if model:
            run=output/'VascularMD'/tag
            qc_path=run/(Path(modeled_input['path']).stem+'_vmd_qc.json')
            if qc_path.exists():
                report=json.loads(qc_path.read_text())
                if report['source_sha256']!=modeled_input['sha256']:raise ValueError('VascularMD source hash mismatch')
            else:
                # Preserve an interrupted native run; never mix its partial files
                # into a new attempt. Completed failures are reported, not retried.
                if run.exists() and any(run.iterdir()):
                    import time
                    run.rename(run.with_name(tag+'_interrupted_'+str(time.time_ns())))
                report=process_file(Path(modeled_input['path']),run,Options(surface=True,accept_native_merges=False))
            item['VascularMD_report']=str(next(run.glob('*_qc.json')).resolve())
            item['VascularMD_status']='VASCULARMD_MODELED' if report['status']=='success' else 'VASCULARMD_MODEL_FAILED'
            item['model_network_success']=report['status']=='success'
            item['model_outputs']=report['outputs'];item['model_failures']=report['failed_components']
            if report['status']=='success':
                import pyvista as pv
                surface=pv.read(report['outputs']['surface_vtk'])
                stl=run/(tag+'_surface.stl')
                if not stl.exists():surface.save(stl)
                check=pv.read(stl)
                if not check.n_cells or not np.isfinite(check.points).all():raise ValueError('Invalid STL round-trip')
                item['model_outputs']['surface_stl']=str(stl.resolve())
                item['surface_cells']=surface.n_cells
                item['status']='VASCULARMD_MODELED'
        reports.append(item)
        print(f'{output.name} {tag}: {item["VascularMD_status"]}',flush=True)
    write_csv(output/'node_mapping.csv',all_nodes);write_csv(output/'branch_mapping.csv',all_branches)
    return reports
