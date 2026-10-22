"""One-time Pilot feature freezing and four prespecified support gates."""
from pathlib import Path
import json
import numpy as np
from .topbrain_qc import sha256,write_json
from .topbrain_semantic_points import load
from .semantic_ensemble import nearest_support
from .label_transfer_qc import collapse_labels_to_mevo,evaluate_mevo_binary

MINIMUM_DONORS=3


def apply_gate(features,gate):
    known=(features['valid_donor_count']>=gate['minimum_donors']) & (features['agreement']>=gate['agreement_threshold']) & (features['median_normalized_distance']<=gate['distance_threshold'])
    return np.where(known,features['winner'],0).astype(np.int8)


def read_frozen(output):
    folder=Path(output)/'frozen_nn_pilot_predictions'
    manifest=json.loads((folder/'manifest.json').read_text())
    for path,digest in manifest['source_hashes'].items():
        if sha256(Path(path))!=digest:raise ValueError(f'Frozen Pilot input changed: {path}')
    values=[]
    for entry in manifest['artifacts']:
        path=folder/entry['file']
        if sha256(path)!=entry['sha256']:raise ValueError(f'Frozen Pilot artifact changed: {path}')
        with np.load(path,allow_pickle=False) as source:arrays={k:source[k].copy() for k in source.files}
        values.append((entry,arrays))
    return manifest,values


def freeze_pilot(output):
    output=Path(output).resolve();folder=output/'frozen_nn_pilot_predictions'
    if (folder/'manifest.json').is_file():
        manifest,_=read_frozen(output);print('Frozen NN Pilot cache hit: no point queries or registration',flush=True);return manifest
    old_manifest=output/'mevo_binary_recovered_baseline/manifest.json'
    recovered=json.loads(old_manifest.read_text());source_hashes={str(old_manifest):sha256(old_manifest)}
    pilot_path=output/'pilot_summary.json';pilot=json.loads(pilot_path.read_text())
    source_hashes[str(pilot_path)]=sha256(pilot_path)
    prepared=[]
    def remember(path):source_hashes[str(path)]=sha256(path);return path
    for entry in recovered['artifacts']:
        case,side=entry['case'],entry['side'];path=Path(entry['path'])
        if sha256(path)!=entry['sha256']:raise ValueError('Recovered NN hard predictions changed')
        with np.load(remember(path),allow_pickle=False) as saved:
            original={k:saved[k].copy() for k in saved.files}
        def cached(c):
            path=remember(output/'semantic_cache'/f'MRA{c}_{side}.npz')
            remember(path.with_suffix('.json'));return load(path)
        target=cached(case)
        reg_path=remember(output/'pilot'/f'MRA{case}_{side}'/'donor_registrations.json')
        registrations=json.loads(reg_path.read_text())
        ranked=sorted([r for r in registrations if 'rank' in r],key=lambda r:r['rank'])
        if [r['source_case'] for r in ranked]!=entry['donor_cases']:raise ValueError('Frozen donor ranking differs')
        features=nearest_support([(cached(r['source_case']),r) for r in ranked],target.geometry.points)
        # Support features use the exact representative -> original mapping used
        # by the validated baseline. Persist hard predictions by copying them.
        for key,value in features.items():
            if key=='donor_ids':continue
            features[key]=value[:,target.assignment] if key.startswith('per_donor_') else value[target.assignment]
        replay3=np.array([(features['per_donor_native_vote']==k).sum(axis=0) for k in [1,2,3]]).argmax(axis=0)+1
        np.testing.assert_array_equal(replay3,original['predicted_label'])
        np.testing.assert_array_equal(original['points'],target.native_points)
        np.testing.assert_array_equal(original['truth_for_evaluation_only'],target.native_labels)
        features.update(points=original['points'],truth_for_evaluation_only=original['truth_for_evaluation_only'],
            predicted_label=original['predicted_label'],binary_label=collapse_labels_to_mevo(original['predicted_label']),
            registration_ids=np.array([f'{sha256(reg_path)}:rank{r["rank"]}:MRA{r["source_case"]}_{side}' for r in ranked]))
        prepared.append((case,side,features))
    if [(c,s) for c,s,_ in prepared]!=[(c,s) for c in pilot['target_cases'] for s in ['L','R']]:raise ValueError('Pilot sides differ')
    folder.mkdir(exist_ok=True);artifacts=[]
    for case,side,features in prepared:
        path=folder/f'MRA{case}_{side}.npz'
        with path.open('xb') as stream:np.savez_compressed(stream,**features)
        artifacts.append(dict(case=case,side=side,file=path.name,sha256=sha256(path),sample_count=len(features['points'])))
    manifest=dict(status='FROZEN_NN_PILOT',artifacts=artifacts,source_hashes=source_hashes,
        registration_rerun=False,donor_ranking_rerun=False,hard_predictions_copied=True,
        feature_extraction='Once: nearest distances and per-donor votes from saved transforms at original representatives, expanded with saved assignment',
        distance_normalization='distance / max(registration inlier_rmse, 1e-8); epsilon only, no physical distance threshold',
        binary_baseline='Collapse saved three-class argmax (historical validated NN)',
        production_binary_vote='Combine M2+M3 per donor before voting; retained separately as winner',
        code_hashes={str(Path(__file__).name):sha256(Path(__file__))})
    write_json(folder/'manifest.json',manifest);print('Frozen NN artifacts: 10 sides, original hard predictions unchanged',flush=True)
    return manifest


def calibrate(output):
    output=Path(output).resolve();destination=output/'nn_production/unknown_gate_summary.json'
    frozen,items=read_frozen(output)
    if destination.exists():
        report=json.loads(destination.read_text())
        if report['frozen_manifest_sha256']!=sha256(output/'frozen_nn_pilot_predictions/manifest.json'):raise ValueError('Support reference changed')
        print('Support calibration cache hit: no repeated gate evaluation',flush=True);return report
    correct=[];agreements=[]
    for entry,arrays in items:
        truth=collapse_labels_to_mevo(arrays['truth_for_evaluation_only'],truth=True)
        mask=arrays['binary_label']==truth
        correct.extend(arrays['median_normalized_distance'][mask]);agreements.extend(arrays['agreement'][mask])
    reference={str(q):float(np.percentile(correct,q)) for q in [95,99]}
    fields=['MeVO_f1','false_inclusion','MeVO_recall','false_exclusion_total','UNKNOWN_fraction']
    baseline=json.loads((output/'mevo_binary_comparison.json').read_text())['baseline']
    candidates=[]
    for agreement in [.6,.8]:
        for percentile in [95,99]:
            gate=dict(agreement_threshold=agreement,distance_percentile=percentile,
                distance_threshold=reference[str(percentile)],minimum_donors=MINIMUM_DONORS)
            per_side=[]
            for entry,arrays in items:
                result=evaluate_mevo_binary(arrays['truth_for_evaluation_only'],apply_gate(arrays,gate))
                per_side.append(dict(case=entry['case'],side=entry['side'],**{k:result[k] for k in fields}))
            metrics={k:float(np.mean([r[k] for r in per_side])) for k in fields}
            candidates.append(dict(gate=gate,metrics=metrics,per_side=per_side,
                accepted=metrics['MeVO_f1']>=baseline['MeVO_f1']-.01 and metrics['false_inclusion']<baseline['false_inclusion']-1e-12))
    allowed=[r for r in candidates if r['accepted']]
    chosen=min(allowed,key=lambda r:(r['metrics']['UNKNOWN_fraction'],-r['gate']['distance_percentile'],r['gate']['agreement_threshold'])) if allowed else candidates[1]
    report=dict(status='CLASSIFICATION_VALIDATED_ON_PILOT' if allowed else 'OOD_SUPPORT_GATE_ONLY',
        optimization_status='PILOT_CALIBRATED_NOT_INDEPENDENT_VALIDATION' if allowed else 'NOT_CLASSIFICATION_OPTIMIZED',
        frozen_manifest_sha256=sha256(output/'frozen_nn_pilot_predictions/manifest.json'),reference_count=len(correct),
        reference_distance_percentiles=reference,winning_agreement_reference=dict(min=float(np.min(agreements)),median=float(np.median(agreements))),
        candidates=candidates,selected=chosen,baseline=baseline,classification_rejection_enabled=bool(allowed),
        classification_disabled_fallback='Historical NN binary classifier' if not allowed else None,
        registration_rerun=False,donor_ranking_rerun=False,gate_combinations_evaluated=4,
        reference_policy='Correct frozen historical NN binary predictions, all original semantic samples, pooled; same Pilot used for calibration only')
    destination.parent.mkdir(exist_ok=True);write_json(destination,report)
    print(json.dumps(dict(status=report['status'],selected=chosen),indent=2),flush=True)
    return report
