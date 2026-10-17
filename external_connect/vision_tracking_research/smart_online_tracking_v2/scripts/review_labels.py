"""Audit actual human revisions; freeze only on explicit post-annotation invocation."""
from pathlib import Path
import sys,json,argparse,datetime
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from annotation_store import load_annotation,sha
from runtime import identity
from annotation_regions import reviewed_boxes,reviewed_counts
from annotation_snapshot import freeze_reviewed,export_frozen,dataset_readiness

def audit():
    rows=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows'];review=[]
    for r in rows:
        if sha(ROOT/r['image'])!=r['image_sha256']:raise ValueError('Candidate image changed: '+r['id'])
        meta,masks=load_annotation(ROOT/'data/annotations',r)
        review.append({'id':r['id'],'split':r['split'],'state':meta['status'] if meta else 'UNREVIEWED',
                       'expert_confirmed':bool(meta and meta['expert_confirmed']),
                       'review_scope':meta.get('review_scope','FULL_FRAME') if meta else None,
                       'reviewed_boxes_yxyx':reviewed_boxes(meta,masks) if meta else [],
                       'reviewed_counts':reviewed_counts(meta,masks) if meta else {},
                       'reviewer':meta['reviewer'] if meta else None,'counts':meta['instance_counts'] if meta else {},
                       'annotation_revision':meta['revision'] if meta else None,
                       'mask_file':meta['mask_file'] if meta else None,
                       'mask_sha256':meta['mask_sha256'] if meta else None,
                       'metadata_sha256':sha(ROOT/'data/annotations'/r['id']/(meta['revision']+'.json')) if meta else None})
    return rows,review

def main():
    p=argparse.ArgumentParser();p.add_argument('--freeze-name');p.add_argument('--export-train',action='store_true');a=p.parse_args()
    rows,review=audit();confirmed=[r for r in review if r['expert_confirmed']]
    summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'candidates':len(rows),
             'expert_confirmed_frames':len(confirmed),
             'region_reviewed_frames':sum(bool(r['reviewed_boxes_yxyx']) and not r['expert_confirmed'] for r in review),'draft_frames':sum(r['state']=='DRAFT_UNREVIEWED' for r in review),
             'confirmed_by_split':{k:sum(r['split']==k for r in confirmed) for k in ['train','validation','test','QUARANTINE_UNASSIGNED']},
             'status':'EXPERT_GROUND_TRUTH_PENDING' if not any(r['reviewed_boxes_yxyx'] for r in review) else 'PARTIAL_EXPERT_LABELS_REVIEW_REQUIRED',
             'training_run':False,'accuracy_evaluated':False,'frames':review}
    out=ROOT/'results/annotation_audit_current.json';out.write_text(json.dumps(summary,indent=2))
    if a.freeze_name:
        if not any(r['reviewed_boxes_yxyx'] for r in review):raise RuntimeError('No reviewed images or regions; no freeze/training allowed')
        if Path(a.freeze_name).name!=a.freeze_name:raise ValueError('Freeze name must be a simple directory name')
        freeze=ROOT/'results/expert_freezes'/a.freeze_name
        frozen=freeze_reviewed(ROOT,freeze,identity(ROOT))
        (freeze/'readiness.json').write_text(json.dumps(dataset_readiness(ROOT),indent=2))
        if a.export_train:
            export_frozen(freeze)
        print('Frozen',len(frozen['rows']),'expert frames in',freeze)
    else:print(json.dumps({k:v for k,v in summary.items() if k!='frames'},indent=2))
if __name__=='__main__':main()
