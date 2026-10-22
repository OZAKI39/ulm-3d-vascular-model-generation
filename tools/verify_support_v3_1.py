"""Verify completed V3.1 artifacts without writing to any historical output."""
from pathlib import Path
import sys,json,csv,zipfile,subprocess
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from vascular_processing.support_inputs import ROOT,discover,verify,save,digest


def main():
    base=discover()['root'].parent/'vascular_support_planner';a=base/'inverted_study';b=base/'v3_1'
    summary_path=b/'support_orientation_and_closure_summary.json'
    if not summary_path.exists():summary_path=a/summary_path.name
    s=json.loads(summary_path.read_text());checks=[]
    def check(name,value,detail=None):checks.append(dict(check=name,passed=bool(value),detail=detail))
    protection=verify(json.loads((a/'validation/protected_before.json').read_text()))
    check('upstream_v1_v2_v3_sha256',protection['all_unchanged'],protection)
    rows=json.loads((a/'orientation_candidates.json').read_text());check('14_orientation_matrices',len(rows)==14 and all(r['label_matrix_matches'] for r in rows))
    ids=json.loads((a/'full_candidate_ids.json').read_text());check('required_controls_retained','INV_0' in ids and 'UPRIGHT_O000' in ids and len(ids)>=7)
    for name in ['orientation_matrix','surrogate_metrics','full_candidate_metrics','direct_buildplate_support','box_ceiling_burden','plate_release_results','bambu_results','inverted_vs_O000']:
        p=a/(name+'.csv');check('csv_'+name,p.is_file() and p.stat().st_size>0)
    from PIL import Image
    for path,count in [(a/'QC_INVERTED',14)]+([(b/'QC_V3_1',12)] if s['phase_b']['executed'] else []):
        files=list(path.glob('*.png'));check(str(path.relative_to(base))+'_count',len(files)==count)
        for p in files:
            with Image.open(p) as image:
                image.verify()
            check('image_'+p.name,p.stat().st_size>1000)
    folders=[a/'candidates'/i for i in ids]
    if s['phase_b']['executed']:
        import pickle
        previous_ids=None
        for item in s['phase_b']['iterations']:
            folder=Path(item['folder']);folders.append(folder)
            with (folder/'rib_geometry_cache.pkl').open('rb') as f:ids_now={r['rib_id'] for r in pickle.load(f)}
            check(folder.name+'_baseline_ribs_preserved',previous_ids is None or ids_now==previous_ids)
            previous_ids=ids_now
            candidate=json.loads((folder/'candidate.json').read_text())
            if candidate.get('repair_iteration'):
                with (folder/'repair_actions.csv').open() as f:actions=list(csv.DictReader(f))
                check(folder.name+'_accepted_actions_match_export',candidate['repair_iteration']['accepted']==sum(r['accepted']=='True' for r in actions))
            for name in ['custom_support.stl','actual_coverage.csv','vascular_residuals.csv','support_self_residuals.csv','repair_actions.csv','bambu_validation.json','QC/actual_coverage.png']:
                check(folder.name+'/'+name,(folder/name).is_file())
        check('bounded_repair_count',s['phase_b']['iteration_count']<=5)
        check('residual_map',(b/'residual_failure_map.json').is_file())
        with Image.open(b/'actual_coverage_before_after.png') as image:image.verify()
    for folder in folders:
        r=json.loads((folder/'candidate.json').read_text());result=r['slice']
        for label,file in [('full','X1C_project.3mf'),('support','bambu_support_only/sliced.3mf')]:
            with zipfile.ZipFile(folder/file) as z:
                settings=json.loads(z.read('Metadata/project_settings.config'));bad=z.testzip()
                check(folder.name+'_'+label+'_zip_crc',bad is None)
                check(folder.name+'_'+label+'_auto_off',str(settings.get('enable_support')) in ('0','false','False'))
                check(folder.name+'_'+label+'_actual_gcode',any(p.endswith('.gcode') for p in z.namelist()))
        check(folder.name+'_native_support_zero',result['native_support_toolpaths']==0)
        check(folder.name+'_scale_one',r['scale_factor']==1.)
        check(folder.name+'_four_semantic_classes',set(r['semantic_coverage'])=={'VASCULAR','PORT','BOX_BOTTOM','BOX_WALLS'})
        check(folder.name+'_actual_slice',result['slicer_success'] and result['support_only_slicer_success'])
    check('no_false_production',s['production_selected'] is None if s['status']=='NEEDS_ADJUSTMENT' else s['production_selected'] is not None)
    check('no_printer_commands',s['printer_commands_sent'] is False)
    validation=summary_path.parent/'validation';validation.mkdir(exist_ok=True)
    test=subprocess.run([str(ROOT/'run_python.sh'),'-m','pytest','-q','tests/test_s1_8_vascular_support_planner.py','--junitxml='+str(validation/'pytest.xml')],cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (validation/'pytest.log').write_text(test.stdout);check('persistent_regression',test.returncode==0,test.stdout.splitlines()[-1:])
    protection=verify(json.loads((a/'validation/protected_before.json').read_text()))
    result=dict(artifact_workflow_passed=all(c['passed'] for c in checks),production_status=s['status'],checks=checks,source_protection=protection,
        code_sha256={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'vascular_processing').glob('support_v3_1*.py'))})
    save(validation/'artifact_verification.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ('checks','code_sha256')},ensure_ascii=False,indent=2))
    if not result['artifact_workflow_passed']:raise SystemExit(1)


if __name__=='__main__':main()
