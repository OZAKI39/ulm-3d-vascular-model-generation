#!/usr/bin/env python3
"""Read-only completion audit; writes only its V3 verification record."""
import csv,json,sys,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from PIL import Image
from vascular_processing.support_inputs import ROOT,save,verify,digest
from vascular_processing.support_v3_review import STEMS


def main():
    csv.field_size_limit(16*1024*1024)  # Full interval polygons and swept-path evidence.
    out=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_support_planner/v3'
    summary=json.loads((out/'vascular_support_v3_summary.json').read_text());r=summary['diagnostic_selected']
    tables=['v3_orientation_candidates','v3_critical_intervals','v3_rib_segments','v3_roof_segments','v3_cut_points','v3_fragments','v3_removal_paths','v3_actual_coverage','v3_bambu_warnings','v2_vs_v3']
    counts={name:len(list(csv.DictReader((out/'tables'/f'{name}.csv').open()))) for name in tables}
    images=[]
    for i,stem in enumerate(STEMS,1):
        p=out/'QC'/f'{i:02}_{stem}.png'
        with Image.open(p) as image:image.verify()
        images.append(str(p.relative_to(out)))
    archives=[]
    for p in sorted((out/'candidates').glob('*/X1C_project.3mf')):
        with zipfile.ZipFile(p) as z:
            crc=z.testzip();settings=json.loads(z.read('Metadata/project_settings.config'))
        archives.append(dict(candidate=p.parent.name,crc_passed=crc is None,auto_support_off=str(settings['enable_support'])=='0',sha256=digest(p)))
    protection=verify(json.loads((out/'validation/protected_history_before.json').read_text()))
    suites=list(ET.parse(out/'test_results.xml').getroot().iter('testsuite'))
    tests=dict(count=sum(int(s.get('tests',0)) for s in suites),failures=sum(int(s.get('failures',0))+int(s.get('errors',0)) for s in suites))
    old=json.loads((out.parent/'v2/vascular_support_v2_summary.json').read_text())
    report=(out/'vascular_support_v3_report.md').read_text()
    checks=dict(eight_orientations=counts['v3_orientation_candidates']==8,all_tables_nonempty=all(counts.values()),
        twenty_qc=len(images)==20,top4_real_slices=summary['top4_actual_slices']==4,
        archive_crc_and_auto_off=all(a['crc_passed'] and a['auto_support_off'] for a in archives),
        historical_tests_retained=tests['count']>=158 and tests['failures']==0,upstream_and_history_unchanged=protection['all_unchanged'],
        v2_production_still_empty=old.get('production_selected') is None,
        failed_v3_not_published=not r.get('failed_gates') or (summary['production_selected'] is None and not (out/'production_selected').exists()),
        physical_limitations_reported='NO PHYSICAL COUPON AVAILABLE' in report and 'NOT PHYSICALLY VALIDATED' in report,
        main_stl_closed=r['exported_geometry_audit']['passed'])
    result=dict(artifact_workflow_passed=all(checks.values()),checks=checks,table_rows=counts,images=images,archives=archives,tests=tests,
        protection=protection,production_software_status=summary['status'],production_selected=summary['production_selected'],
        note='Artifact/test workflow completion is separate from manufacturing software gates and physical validation.')
    save(out/'validation/artifact_verification.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['images','archives','protection']},ensure_ascii=False,indent=2))
    if not result['artifact_workflow_passed']:raise SystemExit(1)


if __name__=='__main__':main()
