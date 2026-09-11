"""CPU export regression using explicitly labelled archived smoke data, never a new simulation."""
import sys,shutil,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common import *
from rbc_report import analyze,render_review
from cloud_geometry import read_off,geometry
base=Path('/home/lzy/projects/cloud_compute');out=base/'tests/outputs'/('cpu-export-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'));out.mkdir(parents=True,exist_ok=False)
raw=out/'raw';raw.mkdir();source=Path('/home/lzy/projects/cloud_results/cloud-smoke-20260910T221827Z');shutil.copytree(source/'C_repair',raw/'simulation')
s=read(source/'repair_spec.json');original=read('/home/lzy/projects/mirheo_starter/runs/single_rbc_repair/rbc_repair_20260910T131105Z/gpu/A6_continuous_preparation_half_dt/spec.json')
s['mesh']=original['mesh'];s['steps']=400000;s['prep_sample_steps']=25;s['sample_steps']=10000
s['preparation_criteria']=read('/home/lzy/projects/mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/material_matching.json')['common_preparation_criteria']
v,f=read_off(s['mesh']);s['reference_radius']=geometry(v,f)['a']
execution=dict(status='STOPPED',job_id='CPU_EXPORT_FIXTURE_ARCHIVED_SMOKE',exit_code=0,stop_reason='CPU export test of old 250-step smoke, not this full run',solver_process_wall_s=None)
write(raw/'execution.json',execution);write(raw/'run_status.json',execution);write(raw/'full_spec.json',s);write(raw/'actual_parameters.json',{'fixture_only':True});write(raw/'provenance.json',{'fixture_only':True,'source':str(source)})
result=analyze(raw,s,execution);assert result['CLOUD_RBC_RUN_COMPLETE']=='NOT_COMPLETE';assert result['membership']['strict_impermeability']=='NOT_VERIFIED';assert result['membership']['tested_point_frames']==192
write(raw/'disk_after.json',dict(bytes=sum(p.stat().st_size for p in raw.rglob('*') if p.is_file()),fixture_only=True))
(raw/'README_FIXTURE.txt').write_text('CPU export regression from old cloud-smoke-20260910T221827Z. No new simulation.\n')
seal(raw,execution['job_id'],'STOPPED');before=sha(raw/'results_manifest.json');html=render_review(raw,out/'review');assert validate_result(raw)['manifest_sha256']==before
text=Path(html).read_text();assert 'cdn.plot.ly' not in text.split('<script src=')[-1] if '<script src=' in text else True
assert 'NOT_COMPLETE' in text and 'NOT_VERIFIED' in text
summary=dict(status='PASS',fixture_only=True,no_solver_started=True,output=str(out),html=html,raw_manifest_sha256=before,geometry_frames=len(result['geometry'].get('unreadable_frames',[])),sampled_member_point_frames=result['membership']['tested_point_frames'])
write(base/'records/rbc_full_export_test.json',summary);print(json.dumps(summary))
