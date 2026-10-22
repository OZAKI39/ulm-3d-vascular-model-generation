"""Scientific invariants and actual adapters for TopBrain, with labelled synthetic NIfTI."""
from pathlib import Path
import json
import shutil
import subprocess
import sys

import networkx as nx
import nibabel as nib
import numpy as np
import pytest

from vascular_processing.topbrain_dataset import discover,load_case,save_mask
from vascular_processing.topbrain_labels import OFFICIAL_2025_ARTERIES,parse_itksnap,load_labelmap
from vascular_processing.topbrain_qc import TopBrainError,components,physical_points,sha256
from vascular_processing.topbrain_mevo import EXACT,extract_rois
from vascular_processing.topbrain_geometry import mask_surface
from vascular_processing.topbrain_pipeline import TopBrainOptions,process_case,cache_key

ROOT=Path(__file__).resolve().parents[1]


def make_dataset(root, *, modality='mr', case_id='001', affine=None, labels=None, values=None):
    root=Path(root)
    values=values or OFFICIAL_2025_ARTERIES
    names={v:k for k,v in values.items()}
    if labels is None:
        labels=np.zeros((64,33,96),dtype=np.int16)
        x,y,z=np.indices(labels.shape)
        for side,cx in [('R',18),('L',46)]:
            cylinder=((x-cx)**2+(y-16)**2<=2.1**2)&(z>=5)&(z<=88)
            for name,low,high in [('M1',5,25),('M2',26,55),('M3',56,88)]:
                labels[cylinder&(z>=low)&(z<=high)]=names[f'{side}-{name}']
    if affine is None:
        theta=.31
        affine=np.array([[np.cos(theta)*.45,-np.sin(theta)*.7,0,17],
                         [np.sin(theta)*.45,np.cos(theta)*.7,0,-12],
                         [0,0,1.1,33],[0,0,0,1.]])
    for kind,array in [('images',np.zeros_like(labels,dtype=np.int16)),('labels',labels)]:
        folder=root/f'{kind}Tr_topbrain_{modality}'
        folder.mkdir(parents=True,exist_ok=True)
        name=f'topcow_{modality}_{case_id}'+('_0000' if kind=='images' else '')+'.nii.gz'
        image=nib.Nifti1Image(array,affine)
        image.header.set_xyzt_units('mm')
        nib.save(image,folder/name)
    folder=root/'itksnap_labelmap_txt';folder.mkdir(exist_ok=True)
    path=folder/f'topbrain_{modality}.txt'
    path.write_text('# Synthetic fixture; names use official semantics\n'+''.join(f'{k} 255 0 0 1 1 1 "{v}"\n' for k,v in values.items()))
    return root


@pytest.fixture
def dataset(tmp_path):
    return make_dataset(tmp_path/'training')


@pytest.mark.parametrize('name',['R-M1','L-M1','Acom','R-A1A2','L-A1A2','R-A3','L-A3','R-M2','R-M3','L-M2','L-M3','R-P1P2','L-P1P2','R-Pcom','L-Pcom','R-P3P4','L-P3P4'])
def test_name_resolver_uses_actual_map(dataset,name):
    path=dataset/'itksnap_labelmap_txt/topbrain_mr.txt'
    values={n+100:label for n,label in OFFICIAL_2025_ARTERIES.items()}
    path.write_text(''.join(f'{n} 0 0 0 1 1 1 "{label}"\n' for n,label in values.items()))
    assert parse_itksnap(path).resolve(name)==next(n for n,v in values.items() if v==name)


def test_pairing_separates_ct_and_case_ids(dataset):
    make_dataset(dataset,modality='ct',case_id='002')
    assert [(p.modality,p.case_id) for p in discover(dataset)]==[('ct','002'),('mr','001')]
    with pytest.raises(TopBrainError,match='CASE_NOT_FOUND'):
        load_case(dataset,'001','ct')


def test_duplicate_pair_rejected(dataset):
    source=next((dataset/'labelsTr_topbrain_mr').glob('*.nii.gz'))
    shutil.copy2(source,source.with_name('topbrain_mr_001_seg.nii.gz'))
    with pytest.raises(TopBrainError,match='IMAGE_LABEL_MISMATCH'):
        discover(dataset)


def test_missing_dataset_no_bra_va_fallback(tmp_path):
    proc=subprocess.run([sys.executable,str(ROOT/'s1-3_swc_roi_generate_MeVO.py'),'--topbrain-root',str(tmp_path/'missing'),'--no-gui'],capture_output=True,text=True)
    assert proc.returncode==2 and 'TOPBRAIN_DATA_NOT_FOUND' in proc.stderr


def test_ta36_requires_real_map(tmp_path):
    with pytest.raises(TopBrainError,match='LABELMAP_NOT_FOUND'):
        load_labelmap(tmp_path,'mr','ta36')


def test_affine_world_and_mask_roundtrip(dataset,tmp_path):
    case=load_case(dataset)
    point=np.array([[18,16,30]])
    expected=(case.label_image.affine@np.r_[point[0],1])[:3]
    np.testing.assert_allclose(physical_points(point,case.affine_mm)[0],expected)
    roi=extract_rois(case)[0]
    output=tmp_path/'roi.nii.gz';save_mask(output,roi.mask,case)
    image=nib.load(output)
    np.testing.assert_array_equal(np.asanyarray(image.dataobj),roi.mask)
    np.testing.assert_allclose(image.affine,case.label_image.affine)


@pytest.mark.parametrize('change',['affine','shape','units'])
def test_pair_geometry_mismatch(dataset,change):
    path=next((dataset/'labelsTr_topbrain_mr').glob('*.nii.gz'))
    image=nib.load(path);array=np.asanyarray(image.dataobj).copy();affine=image.affine.copy()
    if change=='affine':affine[0,3]+=1
    if change=='shape':array=array[:-1]
    new=nib.Nifti1Image(array,affine,image.header)
    if change=='units':new.header.set_xyzt_units('micron')
    nib.save(new,path)
    with pytest.raises(TopBrainError,match='IMAGE_LABEL_MISMATCH'):
        load_case(dataset)


def test_exact_mca_no_m1_and_proximal_root(dataset):
    case=load_case(dataset);roi=extract_rois(case)[0]
    np.testing.assert_array_equal(roi.mask,case.mask('R-M2')|case.mask('R-M3'))
    assert not np.any(roi.mask&case.mask('R-M1')) and roi.anatomical_status==EXACT


def test_multiple_components_preserved(dataset):
    case=load_case(dataset)
    labels=case.labels.copy();labels[18,16,40]=0;labels[2,2,2]=case.label_map.resolve('R-M3')
    case.labels=labels
    roi=extract_rois(case)[0]
    assert roi.mask[2,2,2] and roi.status=='MULTIPLE_COMPONENTS_NEEDS_REVIEW'
















def test_surface_transforms_affine(dataset):
    case=load_case(dataset);roi=extract_rois(case)[0]
    mesh=mask_surface(roi.mask,case)
    vox=nib.affines.apply_affine(np.linalg.inv(case.affine_mm),mesh.points)
    pts=np.argwhere(roi.mask)
    np.testing.assert_allclose(vox.min(axis=0),pts.min(axis=0)-.5,atol=1e-5)
    np.testing.assert_allclose(vox.max(axis=0),pts.max(axis=0)+.5,atol=1e-5)


def test_cache_and_source_preservation(dataset,tmp_path):
    case=load_case(dataset);before=sha256(case.paths.label)
    opts=TopBrainOptions(context_surface=False)
    run,rois,report=process_case(case,tmp_path/'out',options=opts)
    assert report['status']=='PASS'
    run2,rois2,report2=process_case(case,tmp_path/'out',options=opts)
    assert run==run2 and report2['cache_hit']
    assert sha256(case.paths.label)==before
    assert cache_key(case,opts)[0]!=cache_key(case,TopBrainOptions())[0]
    next(run.rglob('*_topbrain_mask_surface.vtk')).write_text('corrupted')
    with pytest.raises(TopBrainError,match='CACHE_CORRUPTED'):
        process_case(case,tmp_path/'out',options=opts)




def test_human_and_renderer_protected():
    baseline=json.loads((ROOT/'outputs/topbrain_validation/protection_before.json').read_text())['files']
    allowed={'s1-3_swc_roi_generate_MeVO.py'}
    from vascular_processing.project_paths import preserved_workspace_path
    changed=[name for name,value in baseline.items() if name not in allowed and sha256(preserved_workspace_path(name))!=value]
    assert changed==[]


def test_offscreen_ui_switches_roi(dataset,tmp_path):
    from vascular_processing.topbrain_display_adapter import TopBrainViewer
    case=load_case(dataset);run,rois,manifest=process_case(case,tmp_path/'out')
    viewer=TopBrainViewer(case,rois,run,manifest,show=False)
    assert viewer.plotter.shape==(1,2)
    viewer.select(1);viewer.next_group()
    report=viewer.run_window(show=False)
    assert report['status']=='OFFSCREEN_RENDERED'
    assert 'ROI:LMCA_M2M3' in report['events']
    assert Path(report['preview']).stat().st_size>1000


def test_display_envelope_converts_mm_once(dataset,tmp_path):
    from vascular_processing.topbrain_display_adapter import TopBrainViewer
    case=load_case(dataset);run,rois,manifest=process_case(case,tmp_path/'out')
    viewer=TopBrainViewer(case,rois,run,manifest,show=False)
    try:
        record=viewer.records[0]
        original=mask_surface(rois[0].mask,case)
        np.testing.assert_allclose(record.bbox_min_um,np.asarray(original.bounds).reshape(3,2)[:,0]*1000.)
        np.testing.assert_allclose(record.bbox_max_um,np.asarray(original.bounds).reshape(3,2)[:,1]*1000.)
        np.testing.assert_allclose(viewer.surfaces[0].points.min(axis=0),record.bbox_min_um)
    finally:
        for control in viewer.controllers:
            if control is not None:control.dispose()
        viewer.plotter.close()


def test_official_missing_label_units_use_matching_image(dataset):
    path=next((dataset/'labelsTr_topbrain_mr').glob('*.nii.gz'))
    original=nib.load(path)
    image=nib.Nifti1Image(np.asanyarray(original.dataobj).copy(),original.affine,original.header)
    image.header.set_xyzt_units('unknown');nib.save(image,path)
    case=load_case(dataset)
    assert case.provenance['effective_spatial_units']=='mm'
    assert 'LABEL_UNITS_FROM_PAIRED_IMAGE' in case.provenance['unit_resolution']
    np.testing.assert_allclose(case.affine_mm,original.affine)


def test_ta36_does_not_reuse_2025_map(dataset):
    with pytest.raises(TopBrainError,match='LABELMAP_NOT_FOUND'):
        load_labelmap(dataset,'mr','ta36')


def test_ta36_loader_uses_its_actual_integer_values(tmp_path):
    root=make_dataset(tmp_path/'data',values={k+100:v for k,v in OFFICIAL_2025_ARTERIES.items()})
    (root/'labelsTr_topbrain_mr').rename(root/'labelsTr_topbrain_v2_topaneu36class')
    (root/'itksnap_labelmap_txt/topbrain_mr.txt').rename(root/'itksnap_labelmap_txt/TA36_mr.txt')
    case=load_case(root,version='ta36')
    assert case.paths.version=='ta36' and case.label_map.resolve('R-M2')==117
    assert extract_rois(case)[0].mask.any()


def test_contact_crop_identical_to_full_grid():
    from scipy import ndimage as ndi
    from vascular_processing.topbrain_qc import contact,CONNECTIVITY
    rng=np.random.default_rng(82)
    a=rng.random((25,23,29))>.8;b=np.zeros_like(a);b[2:8,5:14,1:9]=rng.random((6,9,8))>.8
    np.testing.assert_array_equal(contact(a,b),a&ndi.binary_dilation(b,structure=CONNECTIVITY))




def test_case_request_callback_only_schedules(dataset,tmp_path):
    from vascular_processing.topbrain_display_adapter import TopBrainViewer
    case=load_case(dataset);run,rois,manifest=process_case(case,tmp_path/'out')
    viewer=TopBrainViewer(case,rois,run,manifest,show=False)
    viewer.request_case(1)
    assert viewer.case_delta==1 and viewer.events[-1]=='case_request:1'
    for control in viewer.controllers:
        if control is not None:control.dispose()
    viewer.plotter.close()


def test_partial_mca_union_is_preserved_with_missing_segment(dataset):
    case=load_case(dataset)
    labels=case.labels.copy();labels[labels==case.label_map.resolve('L-M3')]=0;case.labels=labels
    roi=extract_rois(case)[1]
    assert roi.status=='MISSING_SEGMENT' and roi.qc['native_voxel_counts']['L-M3']==0
    np.testing.assert_array_equal(roi.mask,case.mask('L-M2'))






@pytest.mark.parametrize('modality',['mr','ct'])
def test_real_official_case_001(modality):
    data=ROOT/'data/TopBrain'
    if not data.exists():pytest.skip('Official dataset not installed; run explicit provisioning tool')
    case=load_case(data,'001',modality,version='2025')
    assert case.label_map.path.is_relative_to(data.resolve())
    assert case.provenance['dataset_release']['archive_md5']=='b703ea31cd1f0e7115a5d3e6e61f59b3'
    rois=extract_rois(case)
    for roi in rois[:2]:
        assert roi.anatomical_status==EXACT
        if roi.mask is not None:
            np.testing.assert_array_equal(roi.mask,case.mask(f'{roi.side}-M2')|case.mask(f'{roi.side}-M3'))
    assert rois[0].qc['voxel_count']>0
