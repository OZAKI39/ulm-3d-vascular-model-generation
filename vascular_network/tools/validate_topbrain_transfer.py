#!/usr/bin/env python3
"""Small real leave-one-case-out pilot; no target internal labels enter fitting."""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import argparse
import json
import resource
import time
import numpy as np
from vascular_processing.topbrain_semantic_points import build_pool
from vascular_processing.semantic_ensemble import select_donors,ensemble,nearest_baseline
from vascular_processing.partial_fgw_transfer import transfer
from vascular_processing.label_transfer_qc import scores,boundary_error,pilot_gate,GATES
from vascular_processing.topbrain_qc import write_json,sha256
from vascular_processing.transfer_runtime import require_versions


def evaluate_existing(output):
    """Re-evaluate saved predictions' metrics; no registration/transport rerun."""
    path=output/'pilot_summary.json';summary=json.loads(path.read_text())
    for run in summary['all_runs']:
        run['passed']=pilot_gate(run['metrics'],summary['baseline']) and not summary['failures']
        detail=output/'pilot'/f"alpha{run['alpha']}_m{run['mass_fraction']}_summary.json"
        saved=json.loads(detail.read_text());saved['passed']=run['passed'];write_json(detail,saved)
    allowed=[r for r in summary['all_runs'] if r['passed']]
    chosen=max(allowed or summary['all_runs'],key=lambda r:r['metrics']['macro_f1'])
    summary.update(pilot_passed=chosen['passed'],production_permission=chosen['passed'],gates=GATES,
        status='TOPBRAIN_TRANSFER_VALIDATED' if chosen['passed'] else 'TOPBRAIN_LABEL_TRANSFER_NOT_VALIDATED',
        final_parameters=dict(alpha=chosen['alpha'],mass_fraction=chosen['mass_fraction']),metrics=chosen['metrics'],boundaries=chosen['boundaries'],
        baseline_difference=chosen['metrics']['macro_f1']-summary['baseline']['macro_f1'],
        decision_basis='User section 100 requires FGW to improve on the NN baseline. No allowable baseline drop. Saved numeric results unchanged.',
        re_evaluation_only=True)
    summary['source_code_hashes']['vascular_processing/label_transfer_qc.py']=sha256(ROOT/'vascular_processing/label_transfer_qc.py')
    summary['source_code_hashes']['tools/validate_topbrain_transfer.py']=sha256(Path(__file__))
    write_json(path,summary)
    # Read-only semantic point visualization, never a TopBrain centerline/tree.
    import pyvista as pv
    pattern=f"alpha{chosen['alpha']}_m{chosen['mass_fraction']}_predictions.npz"
    for prediction in (output/'pilot').glob('MRA*/'+pattern):
        arrays=np.load(prediction);cloud=pv.PolyData(arrays['points'])
        for name in ['predicted_label','truth_for_evaluation_only','confidence','support_mass']:cloud[name]=arrays[name]
        for i,name in enumerate(['p_M1','p_M2','p_M3']):cloud[name]=arrays['probabilities'][:,i]
        cloud.save(prediction.parent/'semantic_validation_points.vtp')
    print(summary['status'],summary['baseline_difference'],flush=True)
    return 0 if summary['pilot_passed'] else 2


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--topbrain-root',type=Path,default=ROOT/'data/TopBrain')
    p.add_argument('--output-dir',type=Path,default=ROOT/'outputs/topbrain_brava_transfer')
    group=p.add_mutually_exclusive_group();group.add_argument('--pilot',action='store_true');group.add_argument('--full-validation',action='store_true')
    p.add_argument('--cache-only',action='store_true');p.add_argument('--maximum-samples',type=int,default=275)
    p.add_argument('--debug-save-matrices',action='store_true')
    p.add_argument('--evaluate-existing',action='store_true',help='Re-evaluate saved metrics and export point VTPs; does not rerun fitting')
    a=p.parse_args();require_versions()
    if a.evaluate_existing:return evaluate_existing(a.output_dir)
    if not 50<=a.maximum_samples<=300:raise ValueError('Expected 50 to 300 samples')
    output=a.output_dir;output.mkdir(parents=True,exist_ok=True)
    pool,excluded=build_pool(a.topbrain_root,output,a.maximum_samples)
    if a.cache_only:return 0
    eligible=sorted({x.geometry.case_id for x in pool}|{x['case_id'] for x in excluded if x['status']!='MISSING_NATIVE_SEGMENT'})
    targets=eligible if a.full_validation else eligible[:5]
    prepared=[];failures=[]
    for case_id in targets:
        for side in ['L','R']:
            target=next((s for s in pool if s.geometry.case_id==case_id and s.geometry.side==side),None)
            if target is None:
                cause=next((e for e in excluded if e['case_id']==case_id and e['side']==side),{})
                if cause.get('status')!='MISSING_NATIVE_SEGMENT':failures.append(dict(case_id=case_id,side=side,status=cause.get('status','UNAVAILABLE')))
                continue
            start=time.perf_counter();folder=output/'pilot'/f'MRA{case_id}_{side}';folder.mkdir(parents=True,exist_ok=True)
            selected,registrations=select_donors(pool,target.geometry)
            write_json(folder/'donor_registrations.json',registrations)
            if not selected:
                failures.append(dict(case_id=case_id,side=side,status='INVALID_REGISTRATION'));continue
            baseline=nearest_baseline(selected,target.geometry)
            base_full=baseline[target.assignment]
            baseline_scores=scores(target.native_labels,base_full)
            prepared.append(dict(target=target,selected=selected,folder=folder,baseline=baseline_scores,
                                 registration_seconds=time.perf_counter()-start))
            print(f'MRA{case_id}/{side}: {sum(r["status"]=="VALID" for r in registrations)} valid registrations; {len(selected)} donors; baseline F1 {baseline_scores["macro_f1"]:.4f}',flush=True)
    if not prepared:raise RuntimeError('No valid target alignment; registration requires review')
    baseline={k:float(np.mean([x['baseline'][k] for x in prepared])) for k in prepared[0]['baseline']}
    grid=[(.5,.8),(.4,.75),(.4,.9),(.7,.75),(.7,.9)];runs=[]
    for alpha,mass in grid:
        records=[]
        for item in prepared:
            start=time.perf_counter();target=item['target'];folder=item['folder'];results=[]
            for donor,registration in item['selected']:
                debug=folder/f'alpha{alpha}_m{mass}_donor{donor.geometry.case_id}.npz' if a.debug_save_matrices else None
                result=transfer(donor,target.geometry,registration,alpha=alpha,mass_fraction=mass,debug_path=debug)
                results.append(result)
            combined=ensemble(results);prediction=combined['predicted_label'][target.assignment]
            metric=scores(target.native_labels,prediction)
            record=dict(case_id=target.geometry.case_id,side=target.geometry.side,metrics=metric,
                        M1_M2_boundary_median_error_mm=boundary_error(target.native_points,target.native_labels,prediction,(1,2)),
                        M2_M3_boundary_median_error_mm=boundary_error(target.native_points,target.native_labels,prediction,(2,3)),
                        valid_donors=len(results),donor_cases=[x['source_case'] for x in results],fgw_losses=[x['loss'] for x in results],
                        runtime_seconds=time.perf_counter()-start+item['registration_seconds'],process_peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                        baseline=item['baseline'])
            records.append(record)
            name=f'alpha{alpha}_m{mass}'
            np.savez_compressed(folder/(name+'_predictions.npz'),points=target.native_points,predicted_label=prediction,
                                probabilities=combined['probabilities'][target.assignment],truth_for_evaluation_only=target.native_labels,
                                confidence=combined['confidence'][target.assignment],support_mass=combined['support_mass'][target.assignment])
            write_json(folder/(name+'_metrics.json'),record)
            print(f'{name} MRA{target.geometry.case_id}/{target.geometry.side}: F1={metric["macro_f1"]:.4f}, unknown={metric["UNKNOWN_fraction"]:.3f}',flush=True)
        aggregate={k:float(np.mean([x['metrics'][k] for x in records])) for k in records[0]['metrics']}
        boundary={k:(float(np.median([x[k] for x in records if x[k] is not None])) if any(x[k] is not None for x in records) else None)
                  for k in ['M1_M2_boundary_median_error_mm','M2_M3_boundary_median_error_mm']}
        run=dict(alpha=alpha,mass_fraction=mass,metrics=aggregate,boundaries=boundary,targets=records,
                 passed=pilot_gate(aggregate,baseline) and not failures)
        runs.append(run);write_json(output/'pilot'/f'alpha{alpha}_m{mass}_summary.json',run)
        print(f'Aggregate alpha={alpha} m={mass}: {aggregate}; gate={run["passed"]}',flush=True)
        if len(runs)==1 and run['passed']:break
    # Gate-passing configurations take priority, then balanced class F1.
    admissible=[r for r in runs if r['passed']]
    chosen=max(admissible or runs,key=lambda r:(r['metrics']['macro_f1'],min(r['metrics']['M2_f1'],r['metrics']['M3_f1'])))
    codefiles=[Path(__file__),*list((ROOT/'vascular_processing').glob('*transfer*.py')),
               ROOT/'vascular_processing/topbrain_semantic_points.py',ROOT/'vascular_processing/similarity_registration.py',ROOT/'vascular_processing/semantic_ensemble.py']
    summary=dict(status='TOPBRAIN_TRANSFER_VALIDATED' if chosen['passed'] else 'TOPBRAIN_LABEL_TRANSFER_NOT_VALIDATED',
        pilot_passed=chosen['passed'],target_cases=targets,target_sides=len(prepared)+len(failures),evaluated_target_sides=len(prepared),failures=failures,
        donor_pool_case_count=len({x.geometry.case_id for x in pool}),donor_pool_side_count=len(pool),top_k=5,
        default_parameters=dict(alpha=.5,mass_fraction=.8),parameter_grid_run=len(runs)>1,parameter_runs=len(runs),
        final_parameters=dict(alpha=chosen['alpha'],mass_fraction=chosen['mass_fraction']),baseline=baseline,metrics=chosen['metrics'],
        boundaries=chosen['boundaries'],gates=GATES,all_runs=[{k:v for k,v in run.items() if k!='targets'} for run in runs],
        metric_aggregation='Per-side F1 on all original semantic samples; equal mean across target sides. Missing-anchor failures prohibit admission.',
        boundary_method='Symmetric nearest interface distances on identical local geometric neighbour graph; median per side then median of sides; null denotes no interface prediction',
        validation_leakage_policy='Only union geometry, sampling mass and ICA/MCA root enter target Geometry. Internal native labels used only after predictions.',
        source_code_hashes={str(f.relative_to(ROOT)):sha256(f) for f in codefiles},production_permission=chosen['passed'],
        decision_basis='User section 100: FGW must improve on nearest-neighbor baseline; no baseline drop allowed')
    write_json(output/'pilot_summary.json',summary)
    print(summary['status'],flush=True)
    return 0 if chosen['passed'] else 2

if __name__=='__main__':raise SystemExit(main())
