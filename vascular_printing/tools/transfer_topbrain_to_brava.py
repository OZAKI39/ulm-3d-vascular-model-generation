#!/usr/bin/env python3
"""BG001 proof of concept using frozen Open3D registration and NN support rules."""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import argparse
import json
import numpy as np
from vascular_processing.topbrain_qc import sha256,write_json
from vascular_processing.topbrain_semantic_points import Geometry,downsample,trunk_direction,load
from vascular_processing.semantic_ensemble import select_donors,nearest_support
from vascular_processing.nn_support_gate import apply_gate
from vascular_processing.brava_branch_labels import (load_exact_graph,color_source_mapping,select_major_branches,
    aggregate_branches,topology_warnings,LABEL_NAMES)
from vascular_processing.brava_mevo_roi import write_csv,extract_and_model


def save_labeled_tree(path,exact,graph,branches,rows,node_ids,features,labels):
    import pyvista as pv
    index={int(n):i for i,n in enumerate(node_ids)};edges=list(graph.edges)
    by_id={r['branch_id']:r for r in rows}
    points=np.array([graph.nodes[n]['coords'] for n in node_ids])
    mesh=pv.PolyData(points[:,:3],lines=np.array([[2,index[a],index[b]] for a,b in edges]).ravel())
    for key,value in dict(original_swc_id=node_ids,radius=points[:,3],label=labels,
        swc_type=np.array([graph.nodes[n]['swc_type'] for n in node_ids]),
        p_M1=features['p_M1'],p_MeVO=features['p_MeVO']).items():mesh.point_data[key]=value
    branch_keys=[f'{a}_{b}' for s,t in edges for a,b in [exact.edge_locations[(s,t)]]]
    mesh.cell_data['branch_id']=branch_keys
    mesh.cell_data['label']=[by_id[k]['label_code'] for k in branch_keys]
    mesh.cell_data['is_roi']=[by_id[k]['label']=='MeVO' for k in branch_keys]
    mesh.field_data['label_codes']=['0=UNKNOWN;1=M1;2=MeVO']
    mesh.save(path)


def transfer_side(exact,major_labels,config,side,output,calibration,pool,source_mapping,*,model=True):
    name=side+'MCA';folder=output/name
    manifest_file=folder/'roi_manifest.json'
    if manifest_file.exists():
        saved=json.loads(manifest_file.read_text())
        if saved['source']['raw_sha256']!=sha256(exact.source.path) or saved['source']['color_sha256']!=config['color_sha256'] or saved['UNKNOWN_gate']!=calibration['selected']['gate'] or saved.get('major_artery',{}).get('type_code',config['major_arteries'].get(name))!=config['major_arteries'].get(name):
            raise ValueError('Existing BG001 result belongs to a different source, laterality mapping or gate')
        if saved.get('complete') and (not model or saved.get('model_requested',True)):
            for relative,digest in saved.get('artifact_hashes',{}).items():
                if not (folder/relative).is_file() or sha256(folder/relative)!=digest:
                    raise ValueError('CACHED_ARTIFACT_CHANGED: '+relative)
            print(f'{name}: cached result',flush=True);return saved
    code=config['major_arteries'].get(name)
    if not isinstance(code,int):raise ValueError(f'MAJOR_ARTERY_MAPPING_MISSING: need verified ColorCoded TYPE for {name}')
    graph,branches=select_major_branches(exact,major_labels,code)
    ids=np.array(sorted(graph));points=np.array([graph.nodes[n]['coords'][:3] for n in ids])
    root=next(n for n in graph if graph.in_degree(n)==0);origin=graph.nodes[root]['coords'][:3]
    representatives,assignment,mass=downsample(points,275)
    target=Geometry(points[representatives],mass,origin,trunk_direction(points,origin),'BG001',side)
    folder.mkdir(parents=True,exist_ok=True)
    registration_path=folder/'donor_registrations.json'
    if registration_path.exists():
        registrations=json.loads(registration_path.read_text())
        ranked=sorted([r for r in registrations if 'rank' in r],key=lambda r:r['rank'])
        selected=[(next(d for d in pool if d.geometry.case_id==r['source_case'] and d.geometry.side==side),r) for r in ranked]
    else:
        selected,registrations=select_donors(pool,target)
        write_json(registration_path,registrations)
        write_csv(folder/'donor_registrations.csv',[{k:r.get(k) for k in ['source_case','target_case','side','status','rank','fitness','inlier_rmse','normalized_rmse','scale','root_error']} for r in registrations])
    valid_count=sum(r['status']=='VALID' for r in registrations)
    print(f'BG001 {name}: {valid_count} valid registrations, {len(selected)} selected',flush=True)
    support_path=folder/'point_support.npz'
    if support_path.exists():
        with np.load(support_path,allow_pickle=False) as data:
            np.testing.assert_array_equal(data['original_swc_id'],ids)
            np.testing.assert_array_equal(data['points'],points)
            features={k:data[k] for k in data.files if k not in {'original_swc_id','points','label'}}
            labels=data['label']
        np.testing.assert_array_equal(labels,apply_gate(features,calibration['selected']['gate']))
    else:
        features=nearest_support(selected,points)  # Original samples; empty support is all UNKNOWN.
        labels=apply_gate(features,calibration['selected']['gate'])
        np.savez_compressed(support_path,original_swc_id=ids,points=points,label=labels,**features)
    rows=aggregate_branches(exact,branches,ids,features,labels)
    warnings=topology_warnings(branches,rows)
    write_csv(folder/'branch_predictions.csv',rows)
    point_rows=[]
    memberships={int(n):[] for n in ids}
    for key,b in branches.items():
        for n in b['node_ids']:memberships[n].append(key)
    for i,n in enumerate(ids):
        x,y,z,r=graph.nodes[n]['coords'];parents=list(graph.predecessors(n))
        incoming=f'{exact.edge_locations[(parents[0],n)][0]}_{exact.edge_locations[(parents[0],n)][1]}' if parents else memberships[n][0]
        point_rows.append(dict(original_swc_id=int(n),branch_id=incoming,branch_memberships=';'.join(memberships[n]),
            x=x,y=y,z=z,radius=r,swc_type=graph.nodes[n]['swc_type'],p_M1=features['p_M1'][i],p_MeVO=features['p_MeVO'][i],
            agreement=features['agreement'][i],vote_margin=features['vote_margin'][i],
            normalized_support_distance=features['median_normalized_distance'][i],valid_donor_count=int(features['valid_donor_count'][i]),label=LABEL_NAMES[int(labels[i])]))
    write_csv(folder/'point_predictions.csv',point_rows)
    if not (folder/'labeled_tree.vtp').exists():save_labeled_tree(folder/'labeled_tree.vtp',exact,graph,branches,rows,ids,features,labels)
    def statistics(records):
        return {key:dict(min=float(min(r[key] for r in records)),median=float(np.median([r[key] for r in records])),max=float(max(r[key] for r in records)))
            for key in ['fitness','inlier_rmse','normalized_rmse','scale','root_error'] if records and all(key in r for r in records)}
    registration_summary=statistics([r for _,r in selected])
    report=dict(method='Open3D similarity + multi-donor NN',side=name,source=source_mapping,
        major_artery=dict(type_code=code,evidence=config['evidence'],mapping_selection=config.get('mapping_selection'),
            mapping_basis=config.get('mapping_basis'),evidence_status=config.get('evidence_status'),
            geometry_source='Official standardized non-smoothed SWC; annotation companion only supplies major-tree labels'),
        status='ANATOMICAL_TRANSFER_CANDIDATE',anatomical_status='TopBrain-informed transferred candidate labels; NOT manual anatomical ground truth',
        UNKNOWN_gate=calibration['selected']['gate'],gate_status=calibration['status'],branch_thresholds=dict(known_fraction=.6,probability=.6),
        branch_aggregation='Trapezoidal centerline arc-length weights over original branch samples, including shared junction endpoints',
        registration_summary=registration_summary,attempted_registration_summary=statistics(registrations),
        attempted_registration_donors=len(registrations),valid_registration_donors=valid_count,selected_donor_count=len(selected),donor_cases=[d.geometry.case_id for d,r in selected],
        registration_geometry_samples=len(representatives),point_count=len(ids),branch_count=len(rows),
        point_label_counts={name:int(np.sum(labels==code)) for code,name in LABEL_NAMES.items()},
        branch_label_counts={name:sum(r['label']==name for r in rows) for name in LABEL_NAMES.values()},
        topology_warnings=warnings,original_geometry_preserved=True,manufacturing_status='MANUFACTURING_NOT_VALIDATED')
    # Persist labels before optional modeling; another side is independent.
    write_json(manifest_file,report)
    components=extract_and_model(exact,branches,rows,folder,model=model)
    report.update(roi_components=components,roi_component_count=len(components),
        roi_node_count=sum(r['strict']['node_count'] for r in components),roi_original_branch_count=sum(r['strict']['original_branch_count'] for r in components),
        VascularMD_status='VASCULARMD_MODELED' if any(r['VascularMD_status']=='VASCULARMD_MODELED' for r in components) else 'VASCULARMD_MODEL_FAILED' if components and model else 'NOT_RUN')
    report['status']='BRAVA_MEVO_ROI_EXTRACTED' if components else 'TRANSFER_SUPPORT_INSUFFICIENT'
    if components and report['VascularMD_status']=='VASCULARMD_MODELED':report['status']='VASCULARMD_MODELED'
    if not selected:report['status']='REGISTRATION_UNSUPPORTED'
    report['model_requested']=model
    report['artifact_hashes']={str(p.relative_to(folder)):sha256(p) for p in sorted(folder.rglob('*')) if p.is_file() and p!=manifest_file and '_interrupted_' not in str(p)}
    report['complete']=True
    write_json(manifest_file,report);return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=ROOT/'configs/brava_major_arteries.json')
    p.add_argument('--output-dir',type=Path,default=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001')
    p.add_argument('--side',choices=['L','R','both'],default='both')
    p.add_argument('--no-model',action='store_true')
    a=p.parse_args();config=json.loads(a.config.read_text())
    if config['subject']!='BG001':raise ValueError('Only BG001 is authorized for this proof of concept')
    raw=ROOT/config['raw_source'];color=ROOT/config['color_source']
    if sha256(raw)!=config['raw_sha256'] or sha256(color)!=config['color_sha256']:raise ValueError('Major-tree annotations belong to another SWC')
    exact=load_exact_graph(raw);major_labels,mapping,source_mapping=color_source_mapping(exact,color)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    mapping_path=a.output_dir/'source_node_mapping.csv'
    if not mapping_path.exists():write_csv(mapping_path,[dict(original_swc_id=n,color_swc_id=c,major_type=major_labels[n]) for n,c in mapping.items()])
    base=ROOT/'outputs/topbrain_brava_transfer'
    calibration=json.loads((base/'nn_production/unknown_gate_summary.json').read_text())
    semantic_files=sorted((base/'semantic_cache').glob('*.npz'))
    contract=dict(subject=config['subject'],raw_sha256=config['raw_sha256'],color_sha256=config['color_sha256'],
        major_arteries=config['major_arteries'],mapping_selection=config.get('mapping_selection'),
        UNKNOWN_gate=calibration['selected']['gate'],
        registration_implementation_sha256=sha256(ROOT/'vascular_processing/similarity_registration.py'),
        geometry_implementation_sha256=sha256(ROOT/'vascular_processing/topbrain_semantic_points.py'),
        donor_cache_hashes={str(p.relative_to(ROOT)):sha256(p) for f in semantic_files for p in [f,f.with_suffix('.json')]})
    contract_path=a.output_dir/'input_provenance.json'
    if contract_path.exists() and json.loads(contract_path.read_text())!=contract:
        raise ValueError('INPUT_PROVENANCE_CHANGED: choose a separate output directory for different inputs')
    if not contract_path.exists():write_json(contract_path,contract)
    pool=[load(path) for path in semantic_files]
    results={}
    for side in ['L','R'] if a.side=='both' else [a.side]:
        try:results[side+'MCA']=transfer_side(exact,major_labels,config,side,a.output_dir,calibration,pool,source_mapping,model=not a.no_model)
        except Exception as exc:
            import traceback
            results[side+'MCA']=dict(status='TRANSFER_FAILED',error=str(exc),traceback=traceback.format_exc());print(traceback.format_exc(),flush=True)
    write_json(a.output_dir/'summary.json',results)
    successful=lambda r:r.get('roi_component_count',0)>0 if a.no_model else r.get('VascularMD_status')=='VASCULARMD_MODELED'
    return 0 if any(successful(r) for r in results.values()) else 2

if __name__=='__main__':raise SystemExit(main())
