"""Versioned expert annotation I/O, with explicit draft and approval semantics."""
from pathlib import Path
import json,hashlib,uuid
from datetime import datetime,timezone
import numpy as np
from annotation_regions import normalize_box,region_digest,valid_regions,reviewed_boxes,reviewed_counts

CLASSES=('SINGLE_BUBBLE','AGGREGATE','BACKGROUND_ARTIFACT','UNCERTAIN','OUT_OF_FOCUS')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def validate(masks,shape,confirm=False):
    overlap=np.zeros(shape,np.uint8)
    for name in CLASSES:
        a=masks[name]
        if a.shape!=tuple(shape) or a.dtype.kind not in 'ui':raise ValueError('Invalid instance mask '+name)
        overlap+=(a>0)
    conflicts=int((overlap>1).sum())
    if confirm and conflicts:raise ValueError(f'{conflicts} pixels belong to multiple classes; resolve before confirmation')
    return conflicts

def save_annotation(directory,record,masks,reviewer='',confirm=False,track_ids=None,notes='',*,
                    review_scope='FULL_FRAME',review_regions=None,confirm_region=None,
                    workspace_roi=None,proposal_sources=None):
    shape=(record['resolution_wh'][1],record['resolution_wh'][0])
    if review_scope not in ('FULL_FRAME','REGIONS'):raise ValueError('Unknown review scope')
    if confirm_region is not None and review_scope!='REGIONS':raise ValueError('Region confirmation requires REGIONS scope')
    conflicts=validate(masks,shape,confirm and review_scope=='FULL_FRAME')
    if (confirm or confirm_region is not None) and (not reviewer.strip() or reviewer.startswith('SOFTWARE_TEST')):
        raise ValueError('Expert confirmation requires the actual human reviewer name')
    if confirm and review_scope=='REGIONS' and confirm_region is None:raise ValueError('请先选择需要审核的区域')
    if workspace_roi is not None:workspace_roi=normalize_box(workspace_roi,shape)
    regions=valid_regions({'review_regions':review_regions or []},masks)
    invalidated=[r.get('id') for r in (review_regions or []) if r.get('id') not in {v.get('id') for v in regions}]
    if confirm_region is not None:
        box=normalize_box(confirm_region,shape);y0,x0,y1,x1=box
        validate({k:v[y0:y1,x0:x1] for k,v in masks.items()},(y1-y0,x1-x0),True)
        previous=next((r for r in regions if r['box_yxyx']==box),{})
        regions=[r for r in regions if r['box_yxyx']!=box]
        regions.append({'id':previous.get('id',uuid.uuid4().hex),'box_yxyx':box,
                        'confirmed':True,'reviewer':reviewer.strip(),
                        'reviewed_utc':datetime.now(timezone.utc).isoformat(),'mask_sha256':region_digest(masks,box)})
    full=bool(confirm and review_scope=='FULL_FRAME')
    folder=Path(directory)/record['id'];folder.mkdir(parents=True,exist_ok=True)
    revision=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8]
    maskfile=folder/f'{revision}.npz'
    with maskfile.with_suffix('.tmp').open('wb') as f:
        np.savez_compressed(f,**{k:np.asarray(v,dtype=np.uint32) for k,v in masks.items()})
    maskfile.with_suffix('.tmp').replace(maskfile)
    meta={'schema':2,'candidate_id':record['id'],'image_sha256':record['image_sha256'],
          'source_sha256':record['source_sha256'],'split':record['split'],'sequence_id':record.get('sequence_id'),
          'revision':revision,'utc':datetime.now(timezone.utc).isoformat(),'mask_file':maskfile.name,
          'mask_sha256':sha(maskfile),'reviewer':reviewer.strip(),'expert_confirmed':full,
          'status':'EXPERT_CONFIRMED' if full else 'REGIONS_CONFIRMED' if regions else 'DRAFT_UNREVIEWED',
          'review_scope':review_scope,'review_regions':regions,'invalidated_region_ids':invalidated,
          'workspace_roi_yxyx':workspace_roi,'proposal_sources':proposal_sources or [],
          'classes':list(CLASSES),'class_overlap_pixels':conflicts,
          'instance_counts':{k:int(len(np.unique(v))-int(0 in v)) for k,v in masks.items()},
          'track_ids':track_ids or {},'notes':notes,'units':'pixels','physical_fusion_claim':False}
    meta['reviewed_instance_counts']=reviewed_counts(meta,masks)
    (folder/f'{revision}.json').write_text(json.dumps(meta,indent=2))
    pointer=folder/'latest.tmp';pointer.write_text(json.dumps(meta,indent=2));pointer.replace(folder/'latest.json')
    return meta

def load_annotation(directory,record):
    folder=Path(directory)/record['id'];p=folder/'latest.json'
    if not p.exists():return None,None
    meta=json.loads(p.read_text())
    if meta['image_sha256']!=record['image_sha256']:raise ValueError('Image identity mismatch')
    path=folder/meta['mask_file']
    if sha(path)!=meta['mask_sha256']:raise ValueError('Annotation checksum mismatch')
    with np.load(path,allow_pickle=False) as z:masks={k:z[k].copy() for k in CLASSES}
    validate(masks,(record['resolution_wh'][1],record['resolution_wh'][0]),meta['expert_confirmed'])
    if len(valid_regions(meta,masks))!=len(meta.get('review_regions',[])):
        raise ValueError('Regional review checksum mismatch')
    return meta,masks

def export_cellpose(root,records,output,split='train'):
    """Export confirmed whole frames or reviewed crops, never their unreviewed surroundings.

    A crop containing ignore pixels is excluded: standard Cellpose loss does not
    mask these pixels. Boundary instances retain their visible clipped pixels;
    use the recorded review box, not these crops, for whole-object morphology.
    """
    from PIL import Image
    import tifffile,shutil
    root=Path(root);output=Path(output);output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for record in records:
        if record['split']!=split:continue
        meta,masks=load_annotation(root/'data/annotations',record)
        if meta is None:continue
        for region_index,box in enumerate(reviewed_boxes(meta,masks)):
            y0,x0,y1,x1=box;local={k:a[y0:y1,x0:x1] for k,a in masks.items()}
            name=record['id'] if meta['expert_confirmed'] else f'{record["id"]}_roi{region_index:03d}'
            if np.any(local['UNCERTAIN']) or np.any(local['OUT_OF_FOCUS']):
                rows.append({'id':name,'status':'EXCLUDED_IGNORE_REGIONS_FOR_STANDARD_CELLPOSE_LOSS','box_yxyx':box});continue
            combined=np.zeros_like(local['SINGLE_BUBBLE']);n=0;class_map={}
            for cls in ['SINGLE_BUBBLE','AGGREGATE']:
                for label in np.unique(local[cls]):
                    if not label:continue
                    n+=1;combined[local[cls]==label]=n;class_map[str(n)]=cls
            with Image.open(root/record['image']) as im:im.crop((x0,y0,x1,y1)).save(output/f'{name}.png')
            tifffile.imwrite(output/f'{name}_masks.tif',combined)
            rows.append({'id':name,'source_id':record['id'],'status':'EXPORTED','class_map':class_map,
                         'review_scope':meta.get('review_scope','FULL_FRAME'),'box_yxyx':box,
                         'annotation_revision':meta['revision'],'split':record['split']})
    (output/'export_manifest.json').write_text(json.dumps(rows,indent=2))
    return rows
