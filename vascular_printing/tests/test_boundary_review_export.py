"""High-value review integrity tests; no registration, classification or model runs."""
from collections import Counter
import importlib.util
import json
from pathlib import Path
import zipfile

import pytest
from PIL import Image

from vascular_processing.boundary_review_export import (ROOT,HUMAN_SHA,ReviewSourceMismatch,load_cache,
    detect_events,stable_id,build_evidence,separation_pairs,component_summary,protect_inputs,prohibit_recomputation)
from vascular_processing.topbrain_qc import sha256

INPUT=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA'


@pytest.fixture(scope='module')
def cache():
    with prohibit_recomputation():
        value=load_cache(INPUT)
        yield value
        value.verify_unchanged()


def test_boundary_transition_detection_keeps_unknown_and_reverse_events():
    sequence=[('a',None,'M1',[1,2]),('b','a','MeVO',[2,3]),('c','b','UNKNOWN',[3,4]),
              ('d','c','MeVO',[4,5]),('e','d','M1',[5,6]),('f','b','MeVO',[3,7])]
    branches={k:dict(branch_id=k,parent_branch=p,label=label,node_ids=nodes) for k,p,label,nodes in sequence}
    before=json.dumps(branches,sort_keys=True)
    events=detect_events(branches,{1:{},2:{}},{'b':1,'f':1,'d':2},{6,7})
    counts=Counter(e['boundary_type'] for e in events)
    assert counts=={'PROXIMAL_M1_TO_MEVO':1,'DISTAL_MEVO_TO_UNKNOWN':1,'UNKNOWN_TO_MEVO':1,
        'MEVO_TO_M1_REVERSE':1,'DISTAL_MEVO_TO_TERMINAL':1,'MEVO_COMPONENT_ROOT':2,'DETACHED_MEVO_COMPONENT':1}
    assert all(e['review_status']=='UNREVIEWED' for e in events)
    assert next(e for e in events if e['boundary_type']=='MEVO_TO_M1_REVERSE')['priority']=='CRITICAL'
    assert json.dumps(branches,sort_keys=True)==before


def test_five_component_summary_preserves_bifurcating_roots_and_model_failure(cache):
    events,details,local=build_evidence(cache);pairs=separation_pairs(cache)
    rows=component_summary(cache,events,pairs)
    assert len(rows)==5 and sum(r['node_count'] for r in rows)==687 and sum(r['branch_count'] for r in rows)==44
    part3=rows[2]
    assert part3['proximal_branch_id']=='112_113;112_128'
    assert part3['proximal_parent_label']=='UNKNOWN' and part3['separation_reason']=='UNKNOWN_GAP'
    assert part3['VascularMD_status']=='VASCULARMD_MODEL_FAILED'
    assert part3['semantic_candidate_status']=='ANATOMICAL_TRANSFER_CANDIDATE'
    assert len(events)==44 and len({e['boundary_node_id'] for e in events})==31
    assert len(pairs)==10 and all(not p['modification_performed'] for p in pairs)


def test_stable_boundary_id_independent_of_row_order(cache):
    terminal={n for n in cache.node_graph if cache.node_graph.out_degree(n)==0}
    a=detect_events(cache.branches,cache.components,cache.branch_component,terminal)
    b=detect_events(dict(reversed(list(cache.branches.items()))),cache.components,cache.branch_component,terminal)
    assert {e['boundary_id'] for e in a}=={e['boundary_id'] for e in b}
    assert len({e['boundary_id'] for e in a})==len(a)
    assert stable_id('a','b','UNKNOWN_TO_MEVO')!=stable_id('a','b','MEVO_COMPONENT_ROOT')


def test_read_only_write_and_unlink_are_rejected_before_mutation(tmp_path):
    source=tmp_path/'frozen.swc';source.write_text('frozen evidence\n');digest=sha256(source)
    with protect_inputs({str(source):digest}):
        with pytest.raises(PermissionError,match='REVIEW_READ_ONLY_INPUT'):source.write_text('changed')
        with pytest.raises(PermissionError,match='REVIEW_READ_ONLY_INPUT'):source.unlink()
    assert sha256(source)==digest


def exporter():
    spec=importlib.util.spec_from_file_location('review_export_cli',ROOT/'tools/export_bg001_rmca_review_bundle.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_source_hash_mismatch_blocks_before_package_creation(tmp_path):
    source=tmp_path/'BG001/RMCA';source.mkdir(parents=True)
    manifest=json.loads((INPUT/'roi_manifest.json').read_text());manifest['source']['raw_sha256']='0'*64
    (source/'roi_manifest.json').write_text(json.dumps(manifest))
    (source.parent/'input_provenance.json').write_bytes((INPUT.parent/'input_provenance.json').read_bytes())
    out=tmp_path/'review'
    with pytest.raises(ReviewSourceMismatch,match='REVIEW_SOURCE_MISMATCH'):exporter().generate(source,out)
    assert not out.exists() and not list(tmp_path.glob('*.zip'))


def test_review_png_and_local_geometry_smoke(cache,tmp_path):
    from vascular_processing.boundary_review_figures import boundary_figure,write_local_vtp
    import pyvista as pv
    events,details,local=build_evidence(cache);event=next(e for e in events if e['boundary_type']=='UNKNOWN_TO_MEVO')
    key=event['boundary_id'];path=tmp_path/'review.png'
    boundary_figure(cache,event,details[key],local[key],path)
    with Image.open(path) as im:
        assert im.width>=1800 and im.height>=1200
        assert max(im.convert('L').getextrema())-min(im.convert('L').getextrema())>100
    write_local_vtp(cache,event,local[key],tmp_path/'review.vtp')
    mesh=pv.read(tmp_path/'review.vtp')
    assert set(mesh.cell_data['role'])=={'native_branch','boundary_marker','direction_arrow'}
    assert mesh.n_verts==1


def test_zip_manifest_completeness_and_actual_bundle_status():
    """The real export is the integration fixture, not another production run."""
    import hashlib
    archive=INPUT/'BG001_RMCA_external_review_bundle.zip'
    assert archive.is_file(),'Run the real read-only bundle export before this integration test'
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        prefix='BG001_RMCA_external_review/'
        manifest=json.loads(z.read(prefix+'review_manifest.json'))
        assert {prefix+item['file'] for item in manifest['files']}==set(z.namelist())
        for item in manifest['files']:
            if item['file']=='review_manifest.json':continue
            assert hashlib.sha256(z.read(prefix+item['file'])).hexdigest()==item['SHA256']
        self_entry=next(f for f in manifest['files'] if f['file']=='review_manifest.json');expected=self_entry['SHA256'];self_entry['SHA256']=None
        assert hashlib.sha256(json.dumps(manifest,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()==expected
        summary=json.loads(z.read(prefix+'bundle_generation_summary.json'))
        assert summary['total_bundle_size']==archive.stat().st_size<100*1024**2
        assert summary['boundary_count']==44 and summary['PNG_count']==53 and summary['VTP_count']==45
        assert summary['boundary_VTP_count']==44 and summary['file_count']==len(z.namelist())
        assert summary['errors']==[] and summary['recomputation_count']==0 and summary['all_boundaries_unreviewed']


def test_human_hash_unchanged(cache):
    assert sha256(ROOT/'s1-2_swc_roi_generate_human.py')==HUMAN_SHA
    cache.verify_unchanged()
