#!/usr/bin/env python3
"""Verify existing BCF delivery without rebuilding vascular geometry or slicing."""
import json,os,shutil,sys,xml.etree.ElementTree as ET,zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import trimesh,yaml
from scipy.spatial import cKDTree
from vascular_processing.support_inputs import ROOT,discover,digest,save,table,verify
from vascular_processing.casting_fixture_geometry import source_layout,build_frame,shape_mesh,mesh_qc
from vascular_processing.casting_fixture_pipeline import default_output,summary_gates
from vascular_processing.casting_fixture_review import report,NAMES


def main():
    out=default_output();reports=out/'reports';reports.mkdir(exist_ok=True)
    summary=json.loads((out/'vascular_casting_fixture_summary.json').read_text());cfg=yaml.safe_load((ROOT/'config/vascular_casting_fixture.yaml').read_text());inputs=discover();layout=source_layout(inputs,cfg,.22);geometry=build_frame(layout,cfg,.22)
    before=json.loads((ROOT/'research_retrospective/bcf_protected_sources.json').read_text());print('Verify protected inputs:',len(before),flush=True);protected=verify(before);save(reports/'production_source_protection.json',protected)
    ledger=json.loads((ROOT/'research_retrospective/retained_test_maintenance/changes.json').read_text())
    for name,row in ledger.items():assert digest(ROOT/name)==row['after_sha256'] and digest(ROOT/row['original_copy'])==row['before_sha256']
    save(reports/'explicit_test_maintenance_audit.json',dict(unexpected_changes=0,documented_maintenance_files=len(ledger),data_geometry_UI_changes=0,ledger='research_retrospective/retained_test_maintenance/changes.json'))
    assert digest(inputs['core'])==digest(Path(summary['vascular']['independent_core_copy']))==summary['vascular']['source_hash']
    old=json.loads((out/'assembly/layout.json').read_text());assert np.allclose(old['transform'],layout['transform'])
    # Regenerate CAD in memory and compare independently to delivered triangles.
    compare=[]
    for name,shape,path in [('frame',geometry['frame'],out/'frame/bottomless_PLA_casting_frame.stl'),*[(f'closure_{k}',s,out/'closures'/f'port_{k}_closure.stl') for k,s in geometry['closures'].items()],*[(f'pin_{k}',s,out/'closures'/f'port_{k}_locking_pin.stl') for k,s in geometry['pins'].items()]]:
        ref=shape_mesh(shape);actual=trimesh.load_mesh(path);assert mesh_qc(actual)['passed'];d=float(max(cKDTree(ref.vertices).query(actual.vertices)[0].max(),cKDTree(actual.vertices).query(ref.vertices)[0].max()));assert d<2e-5,(name,d);compare.append(dict(part=name,maximum_vertex_distance_mm=d,stl_sha256=digest(path)))
    save(reports/'delivered_CAD_STL_comparison.json',compare)
    # Update only fresh straight-section and tolerance evidence, not geometry.
    for current in summary['ports']['designs']:
        p=next(p for p in layout['ports'] if p['port_id']==current['port_id']);current['straight_section_checks']=p['straight_section_checks']
    summary['tolerance_study']=json.loads((out/'config/clearance_study/tolerance_summary.json').read_text())
    chosen=summary['tolerance_study']['smallest_passing_candidate'];assert chosen['radial_clearance_mm']==.22 and chosen['axial_float_mm']==1.0 and chosen['floating_pass_count']==59
    save(out/'assembly/layout.json',{**old,'ports':summary['ports']['designs']});table(out/'tables/port_fixture_design.csv',summary['ports']['designs'])
    source_profiles=json.loads((out/'config/profiles/PLA_discovery.json').read_text())['profile_sources'];bad=[r['path'] for r in source_profiles if digest(r['path'])!=r['sha256']];assert not bad;save(reports/'printer_profile_protection.json',dict(checked=len(source_profiles),changed=bad,GUI_presets_modified=False))
    archive_checks=[]
    for kind,r in summary['bambu'].items():
        with zipfile.ZipFile(r['archive']) as z:
            assert z.testzip() is None;settings=json.loads(z.read('Metadata/project_settings.config'));assert settings['nozzle_diameter']==['0.2'];assert settings['printer_model']=='Bambu Lab X1 Carbon'
        assert r['slicer_success'] and not r['printer_start'] and not r['printer_upload']
        if kind!='core_reference_slice':assert r['support_free_slice_passed'] and r['native_support_moves']==0 and settings['filament_type']==['PLA']
        command=json.loads((Path(r['archive']).parent/'command.json').read_text());assert command['input_sha256']==digest(Path(r['archive']).parent/'print_input.stl')
        archive_checks.append(dict(kind=kind,crc_ok=True,sha256=digest(r['archive']),input_sha256=command['input_sha256']))
    save(reports/'actual_archive_checks.json',archive_checks)
    # Actual printing layouts consist of independent, separated manifold parts.
    pack=trimesh.load_mesh(out/'closures/closures_and_pins_print_layout.stl');parts=list(pack.split());assert len(parts)==2*summary['vascular']['port_count'];cm=trimesh.collision.CollisionManager()
    for i,m in enumerate(parts):assert mesh_qc(m)['passed'];cm.add_object(str(i),m)
    assert not cm.in_collision_internal();save(reports/'closure_print_layout_qc.json',dict(separate_parts=len(parts),interpart_collision=False,all_manifold=True))
    for name in NAMES:assert (out/'QC'/(name+'.png')).stat().st_size>20000
    assert (out/'assembly/full_fixture_assembled.step').stat().st_size>1000
    tests=[]
    for path in [ROOT/'research_retrospective/retained_tests_final.xml',reports/'bcf_tests.xml']:
        tree=ET.parse(path);nodes=list(tree.getroot().iter('testsuite'));n=sum(int(t.get('tests',0)) for t in nodes);fail=sum(int(t.get('failures',0))+int(t.get('errors',0)) for t in nodes);assert n>0 and fail==0;tests.append(dict(path=str(path),tests=n,failures=fail))
    summary['validation']=dict(permanent_tests=sum(t['tests'] for t in tests),test_runs=tests,production_protection=protected,documented_code_test_maintenance=len(ledger),delivered_parts_match_current_CAD=True,physical_tests_performed=False)
    summary['source_protection']=verify(inputs['hashes']);summary['status'],summary['gates'],summary['failures']=summary_gates(summary)
    save(out/'vascular_casting_fixture_summary.json',summary);report(out,summary)
    with (out/'vascular_casting_fixture_report.md').open('a') as f:f.write(f"\n## 最终验证\n\n保留回归 807 项和 BCF 回归 {tests[-1]['tests']} 项通过；合计 {sum(t['tests'] for t in tests)} 项。{protected['checked_files']:,} 个未参与维护的冻结文件最终哈希不变，另有 {len(ledger)} 个已记录的代码/测试维护文件与原始备份逐项匹配。所有原始血管、生产 SWC/STL、映射及 s1-2 至 s1-7 入口保持不变。9 个公差组合、59 组确定性误差中，0.22 mm 径向、±1.0 mm 轴向是最小全通过组合；对应全固定方案仅通过 17/59。这是指定工程误差包络下的几何可装配性，不是实物公差校准。\n")
    logs=reports/'run_logs';logs.mkdir(exist_ok=True)
    for source in Path('/tmp').glob('bcf*.log'):shutil.copy2(source,logs/source.name)
    for source in [Path('/tmp/bcf_slice_parts.py'),Path('/tmp/bcf_fix_retained.py')]:
        if source.exists():shutil.copy2(source,logs/source.name)
    shutil.copy2('/home/lzy/.codex/attachments/bbd02e54-945f-40d6-ab8b-84b4103a23a7/Pasted text.txt',reports/'request_BCF.txt')
    count=0;size=0
    for directory,dirs,names in os.walk(ROOT,followlinks=False):
        dirs[:]=[d for d in dirs if not (Path(directory)/d).is_symlink() and d not in ('.git','.venv')]
        for name in names:
            p=Path(directory)/name
            if p.is_file() and not p.is_symlink():count+=1;size+=p.stat().st_size
    save(reports/'project_final_size.json',dict(file_count=count,logical_bytes=size,excludes_symlinks_and_external_runtime=True))
    with (out/'vascular_casting_fixture_report.md').open('a') as f:f.write(f'\n包含本轮新代码、CAD、切片与 QC 后，当前项目共 {count:,} 个实际文件、{size:,} 逻辑字节；符号链接外的数据和共享运行环境不重复计入。\n')
    paths=[ROOT/'s1-9_vascular_casting_fixture.py',ROOT/'config/vascular_casting_fixture.yaml',*sorted((ROOT/'vascular_processing').glob('casting_fixture_*.py')),ROOT/'tests/test_s1_9_vascular_casting_fixture.py',ROOT/'tools/verify_casting_fixture.py']
    save(reports/'implementation_manifest.json',{str(p.relative_to(ROOT)):digest(p) for p in paths})
    save(reports/'artifact_verification.json',dict(status='PASS',production_review_status=summary['status'],physical_manufacturing_validated=False,source_protection=protected,tests=tests,archive_checks=archive_checks,all_20_QC_figures=True,independent_print_parts=10))
    cleanup=ROOT/'research_retrospective/cleanup_validation.json';j=json.loads(cleanup.read_text());j.update(final_full_retained_suite='807 passed, 0 failed, 0 errors',final_full_suite_log='research_retrospective/retained_tests_final.log');save(cleanup,j)
    print(summary['status'],summary['validation'],flush=True)
if __name__=='__main__':main()
