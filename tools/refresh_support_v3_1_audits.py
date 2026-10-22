"""Re-audit saved actual toolpaths after validator changes; no re-slicing."""
from pathlib import Path
import sys,json,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from vascular_processing.support_v3_1_pipeline import preparation,production
from vascular_processing.support_v3_1_validation import semantic_audit,removal_audit,gates
from vascular_processing.support_v3_1_closure import compact
from vascular_processing.support_v3_1_review import generate
from vascular_processing.support_inputs import save,verify,digest,ROOT


def main():
    began=time.monotonic();c=preparation();base,out,inputs,profile=c[:4]
    a=json.loads((out/'phase_a_result.json').read_text());records=[]
    audit_hashes={n:digest(ROOT/'vascular_processing'/n) for n in ('support_v3_1_validation.py','support_v3_1_release.py')}
    for old in a['records']:
        folder=Path(old['folder']);r=json.loads((folder/'candidate.json').read_text())
        assert folder.relative_to(base).parts[0]=='inverted_study'
        if a.get('audit_code_sha256')!=audit_hashes:
            semantic_audit(r,inputs,profile);removal_audit(r,inputs)
        r['failed_gates']=gates(r,True)
        r['status']='NEEDS_ADJUSTMENT' if r['failed_gates'] else 'READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED'
        save(folder/'candidate.json',r);records.append(r);print('REFRESH',r['candidate_id'],flush=True)
    passed=[r for r in records if r['candidate_id']!='UPRIGHT_O000' and not r['failed_gates']]
    if passed:raise RuntimeError('NEW_PHASE_A_PASS_REQUIRES_RESELECTING_BEFORE_PHASE_B_REPORT')
    a['records']=records;a['status']='INVERTED_STUDY_FAILED';a['audit_code_sha256']=audit_hashes;save(out/'phase_a_result.json',a);save(out/'full_candidates.json',records)
    if '--phase-a-only' in sys.argv:return
    bpath=base/'v3_1/closure_summary.json';b=json.loads(bpath.read_text());iterations=[];final=None
    for item in b['iterations']:
        folder=Path(item['folder']);r=json.loads((folder/'candidate.json').read_text());assert folder.relative_to(base).parts[0]=='v3_1'
        semantic_audit(r,inputs,profile);r['failed_gates']=gates(r,False);save(folder/'candidate.json',r)
        iterations.append(compact(r));final=r
        path=folder/'iteration_complete.json';done=json.loads(path.read_text());done['metrics']=compact(r);save(path,done)
        print('REFRESH',folder.name,flush=True)
    b['iterations']=iterations;b['winner']=final if not final['failed_gates'] else None
    if b['winner']:production(final,base/'v3_1/production_selected')
    b['status']='READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED' if b['winner'] else 'NEEDS_ADJUSTMENT';save(bpath,b)
    protection=verify(json.loads((out/'validation/protected_before.json').read_text()));save(out/'validation/source_protection.json',protection)
    result=generate(c,records,None,b,protection,time.monotonic()-began);print(result['status'])


if __name__=='__main__':main()
