"""Audit actual human revisions; freeze only on explicit post-annotation invocation."""
from pathlib import Path
import sys,json,argparse,shutil,datetime
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from annotation_store import load_annotation,sha,export_cellpose
from runtime import identity

def audit():
    rows=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows'];review=[]
    for r in rows:
        if sha(ROOT/r['image'])!=r['image_sha256']:raise ValueError('Candidate image changed: '+r['id'])
        meta,masks=load_annotation(ROOT/'data/annotations',r)
        review.append({'id':r['id'],'split':r['split'],'state':meta['status'] if meta else 'UNREVIEWED',
                       'expert_confirmed':bool(meta and meta['expert_confirmed']),
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
             'expert_confirmed_frames':len(confirmed),'draft_frames':sum(r['state']=='DRAFT_UNREVIEWED' for r in review),
             'confirmed_by_split':{k:sum(r['split']==k for r in confirmed) for k in ['train','validation','test','QUARANTINE_UNASSIGNED']},
             'status':'EXPERT_GROUND_TRUTH_PENDING' if not confirmed else 'PARTIAL_EXPERT_LABELS_REVIEW_REQUIRED',
             'training_run':False,'accuracy_evaluated':False,'frames':review}
    out=ROOT/'results/annotation_audit_current.json';out.write_text(json.dumps(summary,indent=2))
    if a.freeze_name:
        if not confirmed:raise RuntimeError('No confirmed expert labels; no freeze/training allowed')
        if Path(a.freeze_name).name!=a.freeze_name:raise ValueError('Freeze name must be a simple directory name')
        freeze=ROOT/'results/expert_freezes'/a.freeze_name;freeze.mkdir(parents=True,exist_ok=False)
        selected=[]
        for r in rows:
            meta,masks=load_annotation(ROOT/'data/annotations',r)
            if not meta or not meta['expert_confirmed'] or r['split']=='QUARANTINE_UNASSIGNED':continue
            source=ROOT/'data/annotations'/r['id']/meta['mask_file']
            target=freeze/source.name;shutil.copy2(source,target)
            selected.append({'record':r,'annotation':meta,'frozen_mask_file':target.name})
        (freeze/'frozen.json').write_text(json.dumps({'execution':identity(ROOT),'protocol':json.loads((ROOT/'protocol.json').read_text()),'rows':selected},indent=2))
        if a.export_train:
            export_cellpose(ROOT,rows,freeze/'cellpose_train',split='train')
            export_cellpose(ROOT,rows,freeze/'cellpose_validation',split='validation')
        print('Frozen',len(selected),'expert frames in',freeze)
    else:print(json.dumps({k:v for k,v in summary.items() if k!='frames'},indent=2))
if __name__=='__main__':main()
