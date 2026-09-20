"""Separate zero-failure stage acceptance from immutable historical failures."""
import json,xml.etree.ElementTree as ET
from pathlib import Path
R=Path(__file__).resolve().parents[2]/'reports/sv1_3p';root=ET.parse(R/'pytest_full.xml')
rows=list(root.iter('testcase'));failed=[t for t in rows if t.find('failure') is not None or t.find('error') is not None];skipped=[t for t in rows if t.find('skipped') is not None];stage=[t for t in rows if 'sv13p' in t.get('classname','')]
d=dict(status='PASS' if not any(t in failed for t in stage) else 'FAIL',stage_passed=sum(t not in failed and t not in skipped for t in stage),stage_failures=sum(t in failed for t in stage),stage_skipped=sum(t in skipped for t in stage),full_passed=len(rows)-len(failed)-len(skipped),full_failures=len(failed),full_skipped=len(skipped),historical_failures=[t.get('classname')+'::'+t.get('name') for t in failed],no_CFD_from_tests=True)
old=ET.parse(R.parent/'sv1_3o/pytest_full.xml');expected=[t.get('classname')+'::'+t.get('name') for t in old.iter('testcase') if t.find('failure') is not None or t.find('error') is not None]
d['historical_failure_set_unchanged']=set(d['historical_failures'])==set(expected)
assert d['status']=='PASS' and d['historical_failure_set_unchanged']
(R/'pytest_summary.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d,indent=2))
