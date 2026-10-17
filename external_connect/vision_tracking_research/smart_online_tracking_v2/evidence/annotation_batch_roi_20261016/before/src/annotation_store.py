"""Versioned expert annotation I/O, with explicit draft and approval semantics."""
from pathlib import Path
import json,hashlib,uuid
from datetime import datetime,timezone
import numpy as np

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

def save_annotation(directory,record,masks,reviewer='',confirm=False,track_ids=None,notes=''):
    shape=(record['resolution_wh'][1],record['resolution_wh'][0])
    conflicts=validate(masks,shape,confirm)
    if confirm and (not reviewer.strip() or reviewer.startswith('SOFTWARE_TEST')):
        raise ValueError('Expert confirmation requires the actual human reviewer name')
    folder=Path(directory)/record['id'];folder.mkdir(parents=True,exist_ok=True)
    revision=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8]
    maskfile=folder/f'{revision}.npz'
    with maskfile.with_suffix('.tmp').open('wb') as f:
        np.savez_compressed(f,**{k:np.asarray(v,dtype=np.uint32) for k,v in masks.items()})
    maskfile.with_suffix('.tmp').replace(maskfile)
    meta={'schema':1,'candidate_id':record['id'],'image_sha256':record['image_sha256'],
          'source_sha256':record['source_sha256'],'split':record['split'],'sequence_id':record.get('sequence_id'),
          'revision':revision,'utc':datetime.now(timezone.utc).isoformat(),'mask_file':maskfile.name,
          'mask_sha256':sha(maskfile),'reviewer':reviewer.strip(),'expert_confirmed':bool(confirm),
          'status':'EXPERT_CONFIRMED' if confirm else 'DRAFT_UNREVIEWED',
          'review_scope':'FULL_FRAME','classes':list(CLASSES),'class_overlap_pixels':conflicts,
          'instance_counts':{k:int(len(np.unique(v))-int(0 in v)) for k,v in masks.items()},
          'track_ids':track_ids or {},'notes':notes,'units':'pixels','physical_fusion_claim':False}
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
    return meta,masks

def export_cellpose(root,records,output,split='train'):
    """Export only explicitly confirmed images; uncertain pixels are recorded separately.

    Training excludes frames with ignore regions by default: Cellpose's standard
    loss does not mask those regions, so exporting them as background is unsafe.
    """
    from PIL import Image
    import tifffile,shutil
    root=Path(root);output=Path(output);output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for record in records:
        if record['split']!=split:continue
        meta,masks=load_annotation(root/'data/annotations',record)
        if meta is None or not meta['expert_confirmed']:continue
        if np.any(masks['UNCERTAIN']) or np.any(masks['OUT_OF_FOCUS']):
            rows.append({'id':record['id'],'status':'EXCLUDED_IGNORE_REGIONS_FOR_STANDARD_CELLPOSE_LOSS'});continue
        combined=np.zeros_like(masks['SINGLE_BUBBLE']);n=0;class_map={}
        for cls in ['SINGLE_BUBBLE','AGGREGATE']:
            for label in np.unique(masks[cls]):
                if not label:continue
                n+=1;combined[masks[cls]==label]=n;class_map[str(n)]=cls
        name=record['id'];shutil.copyfile(root/record['image'],output/f'{name}.png')
        tifffile.imwrite(output/f'{name}_masks.tif',combined)
        rows.append({'id':name,'status':'EXPORTED','class_map':class_map,'annotation_revision':meta['revision']})
    (output/'export_manifest.json').write_text(json.dumps(rows,indent=2))
    return rows
