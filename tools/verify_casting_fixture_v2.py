#!/usr/bin/env python3
"""Read-only source/archive audit; writes evidence only under the V2 output."""
import sys,json,zipfile,shutil
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from vascular_processing.support_inputs import ROOT,digest,save,verify
from vascular_processing.casting_fixture_v2_pipeline import default_output,gates
from vascular_processing.casting_fixture_geometry import rigid,mesh_qc
from vascular_processing.support_slicer_adapter import production_geometry
from vascular_processing.casting_fixture_v2_audit import manifold


def run():
    out=default_output();s=json.loads((out/'internal_wall_supported_fixture_summary.json').read_text());before=json.loads((out/'protection/before.json').read_text());allowed={str(ROOT/'s1-9_vascular_casting_fixture.py'),str(ROOT/'tests/test_s1_9_vascular_casting_fixture.py')};unchanged={p:h for p,h in before.items() if p not in allowed};changes=[]
    for p in allowed:
        original=out/'protection/original'/Path(p).name
        assert digest(original)==before[p]
        if p.endswith('test_s1_9_vascular_casting_fixture.py'):assert Path(p).read_bytes().startswith(original.read_bytes())
        changes.append(dict(path=p,before_sha256=before[p],after_sha256=digest(p),reason='Current entry routes to V2; permanent V2 tests appended, V1 tests preserved verbatim'))
    print('V2 verifying frozen upstream / historical V1 files',flush=True)
    legacy=verify(unchanged);upstream=verify(json.loads((ROOT/'research_retrospective/bcf_protected_sources.json').read_text()))
    original=trimesh.load_mesh(s['vascular']['source']);matrix=np.array(json.loads((out/'assembly/casting_core_transform.json').read_text())['transform_4x4']);rigid(matrix);original.apply_transform(matrix);placed=trimesh.load_mesh(out/'vascular_core/assembly_placed_ABS_core.stl');err=max(cKDTree(original.vertices).query(placed.vertices)[0].max(),cKDTree(placed.vertices).query(original.vertices)[0].max());assert err<1e-5
    assert digest(Path(s['vascular']['source']))==digest(Path(s['vascular']['independent_frozen_copy']))==s['vascular']['source_hash']
    archives={}
    for group in ['frame_slice','retainer_slice']:
        r=s['bambu'][group];folder=Path(r['archive']).parent
        with zipfile.ZipFile(r['archive']) as z:assert z.testzip() is None
        a=production_geometry(r['archive'],folder/'print_input.stl');assert a['passed'];assert r['support_free_slice_passed'];archives[group]=a
    parts=[out/'frame/closed_wall_bottomless_frame_with_support_blades.stl',*(out/'retainers').glob('port_*.stl')];qc={str(p.relative_to(out)):mesh_qc(trimesh.load_mesh(p)) for p in parts};assert all(x['passed'] for x in qc.values())
    diagnostic_qc={}
    for name in ['PDMS_before_trim','PDMS_after_trim']:
        m=trimesh.load_mesh(out/'diagnostic'/(name+'.stl'));assert m.is_watertight and m.volume>0
        boundary_shells=m.split(only_watertight=False);volumes=[float(x.volume) for x in boundary_shells]
        # A fully enclosed vascular void is a negatively oriented inner shell,
        # not a second detached PDMS body. After trim it connects to the outer
        # boundary via the four exposed ports.
        assert sum(v>0 for v in volumes)==1
        if name=='PDMS_after_trim':assert len(boundary_shells)==1
        actual=manifold(m);expected=s['trim']['before_pdms_volume_mm3' if name=='PDMS_before_trim' else 'after_pdms_volume_mm3']
        assert abs(actual.volume()-expected)<.05
        diagnostic_qc[name]=dict(watertight=True,signed_boundary_shell_volumes_mm3=volumes,material_positive_shell_count=1,mesh_volume_mm3=m.volume,boolean_volume_mm3=actual.volume(),stl_volume_tolerance_mm3=.05)
    protected=dict(legacy_and_sources=legacy,all_retained_upstream_files=upstream,authorized_implementation_changes=changes,all_unchanged=legacy['all_unchanged'] and upstream['all_unchanged'])
    save(out/'protection/final.json',protected)
    result=dict(status='PASS',protection=protected,frozen_placed_max_vertex_distance_mm=float(err),archive_geometry=archives,individual_mesh_qc=qc,original_core_sha256=s['vascular']['source_hash'],s1_2_sha256=digest(ROOT/'s1-2_swc_roi_generate_human.py'),s1_3_sha256=digest(ROOT/'s1-3_swc_roi_generate_MeVO.py'),no_printer_upload_or_start=True)
    result['PDMS_diagnostic_meshes']=diagnostic_qc
    assert gates(s)[0]=='READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_SEAL_PENDING'
    save(out/'reports/artifact_verification.json',result)
    files=[ROOT/'s1-9_vascular_casting_fixture.py',ROOT/'config/internal_wall_supported_fixture_v2.yaml',ROOT/'tests/test_s1_9_vascular_casting_fixture.py',Path(__file__).resolve(),*sorted((ROOT/'vascular_processing').glob('casting_fixture_v2_*.py'))]
    save(out/'reports/implementation_manifest.json',{str(p.relative_to(ROOT)):digest(p) for p in files})
    print('V2 artifact verification PASS',len(unchanged)+upstream['checked_files'],'protected checks',flush=True)
    return result

if __name__=='__main__':run()
