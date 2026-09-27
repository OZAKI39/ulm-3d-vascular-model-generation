"""Cached native MCA mask display; TopBrain is never converted into a tree."""
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import nibabel as nib
import numpy as np
from .topbrain_geometry import mask_surface
from .topbrain_dataset import save_mask
from .topbrain_mevo import TopBrainROI, extract_rois
from .topbrain_qc import TopBrainError, sha256, write_json


@dataclass
class TopBrainOptions:
    context_surface: bool = True

    def validate(self):
        if not isinstance(self.context_surface, bool):
            raise ValueError('context_surface must be boolean')


def cache_key(case, options):
    payload = dict(algorithm_version='topbrain-native-mask-2', image=case.provenance['image_sha256'],
                   label=case.provenance['label_sha256'], labelmap=case.label_map.report(),
                   options=asdict(options), annotation_source=case.provenance['annotation_source'],
                   code={p.name:sha256(p) for p in sorted(Path(__file__).parent.glob('topbrain_*.py'))},
                   packages={p:importlib.metadata.version(p) for p in ['nibabel','scikit-image','pyvista']})
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest(), payload


def process_case(case, output_root, *, options=None, recompute=False, cache_root=None):
    options = options or TopBrainOptions(); options.validate()
    key, payload = cache_key(case, options)
    parent = Path(output_root).resolve()/f'topcow_{case.paths.modality}_{case.paths.case_id}'
    run = parent/key[:16]
    if recompute or (run.exists() and not (run/'case_manifest.json').exists()):
        run = parent/(key[:16]+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))
    manifest_path = run/'case_manifest.json'
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        bad = [p for p,h in manifest['artifact_hashes'].items() if not (run/p).is_file() or sha256(run/p)!=h]
        if bad: raise TopBrainError('CACHE_CORRUPTED',f'Use --recompute; changed artifacts: {bad}')
        rois = []
        for name,item in manifest['rois'].items():
            roi = TopBrainROI(name,None,item['side'],item['territory'],item['requested_segments'],
                             item['native_labels_used'],item['anatomical_status'],item['status'],item['warnings'],item['qc'])
            roi.output_dir = run/name
            if 'mask' in item['outputs']:
                roi.mask = np.asanyarray(nib.load(item['outputs']['mask']).dataobj).astype(bool)
            rois.append(roi)
        return run,rois,{**manifest,'cache_hit':True}
    run.mkdir(parents=True,exist_ok=False); write_json(run/'cache_inputs.json',payload)
    rois = extract_rois(case)
    manifest = {**case.provenance,'cache_key':key,'cache_hit':False,'rois':{},'options':asdict(options)}
    if options.context_surface:
        values = [v for v,n in case.label_map.values.items() if n in {'BA','Acom'} or n.startswith(('R-','L-','3rd-'))]
        context = np.isin(case.labels,values)
        if context.any():
            path = run/'native_arterial_context.vtp'; mask_surface(context,case).save(path)
            manifest['context_surface'] = str(path)
    for roi in rois:
        folder = run/roi.name; folder.mkdir(); roi.output_dir = folder; outputs = {}
        if roi.mask is not None and roi.mask.any():
            path = folder/f'{roi.name}_roi.nii.gz'; save_mask(path,roi.mask,case); outputs['mask']=str(path)
            path = folder/f'{roi.name}_topbrain_mask_surface.vtk'; mask_surface(roi.mask,case).save(path)
            outputs['mask_surface'] = str(path)
        item = dict(roi_name=roi.name,side=roi.side,territory=roi.territory,requested_segments=roi.segments,
                    native_labels_used=roi.native_labels,anatomical_status=roi.anatomical_status,
                    status=roi.status,warnings=roi.warnings,qc=roi.qc,outputs=outputs,
                    centerline_export_status='RETIRED_TOPBRAIN_TREE_ROUTE')
        manifest['rois'][roi.name]=item
        write_json(folder/f'{roi.name}_manifest.json',item); write_json(folder/f'{roi.name}_qc.json',roi.qc)
    if any(sha256(getattr(case.paths,k)) != case.provenance[k+'_sha256'] for k in ['image','label']):
        raise TopBrainError('SOURCE_CHANGED','Original NIfTI changed')
    manifest.update(original_nifti_preserved=True,status='PASS' if any(r.status=='PASS' for r in rois) else 'NEEDS_MANUAL_REVIEW',
                    artifact_hashes={str(p.relative_to(run)):sha256(p) for p in run.rglob('*') if p.is_file()})
    write_json(manifest_path,manifest)
    if cache_root: write_json(Path(cache_root)/case.paths.case_id/case.paths.modality/key/'index.json',{'manifest':str(manifest_path)})
    return run,rois,manifest
