from pathlib import Path
import json,shutil,hashlib,subprocess,os
S=Path(__file__).resolve().parents[1];B=S.parent/'bcflux_20260916_224835';P=S/'provenance'
paths=json.loads((P/'STAGE_PATHS.json').read_text());bp=json.loads((P/'BCFLUX_STAGE_PATHS.json').read_text())
for p in (S/'configs').iterdir():
 if p.is_file():p.write_text(p.read_text().replace(bp['remote_stage'],paths['remote_stage']))
p=S/'scripts/run_case.py';p.write_text(p.read_text().replace('_v0_bcflux_','_v0_wallslide_').replace('timeout=600','timeout=1800'))
cfg=(S/'configs/LOW.cfg').read_text().replace('max_time 0.25','max_time 0.8').replace('number_flux 18.046472132892898','number_flux 6.0')
(S/'configs/LONG_TRANSPORT.cfg').write_text(cfg)
# cfg is the authoritative executed configuration; preserve YAML donor as lineage.
for p in (S/'configs').glob('*.yaml'):p.rename(p.with_suffix('.donor.yaml'))
probe=json.loads((P/'BCFLUX_REMOTE_INPUTS_TO_PRESERVE.json').read_text())
for rel,h in json.loads((P/'BCFLUX_IMMUTABLE_FILES.json').read_text()).items():probe['inputs'][bp['remote_stage']+'/'+rel]={'sha256':h}
(P/'REMOTE_INPUTS_TO_PRESERVE.json').write_text(json.dumps(probe,indent=2)+'\n')
shutil.copy2(P/'BCFLUX_BASELINE_CPU_IDENTITY_TO_RECHECK.json',P/'BASELINE_CPU_IDENTITY_TO_RECHECK.json')
git={}
for path in json.loads((P/'BCFLUX_GIT_BEFORE.json').read_text()):
 def g(args):return subprocess.check_output(['git','-C',path]+args,text=True,env={**os.environ,'GIT_OPTIONAL_LOCKS':'0'})
 git[path]={'branch':g(['branch','--show-current']),'commit':g(['rev-parse','HEAD']),'status':g(['status','--porcelain=v1','--untracked-files=all'])}
(P/'GIT_BEFORE.json').write_text(json.dumps(git,indent=2)+'\n')
contract={'model':'KINEMATIC_WALL_CONSTRAINT','hydrodynamics':'R_bulk + R_BB_excess, twist pending; unchanged raw solve','normal':'(X - closest WALL_ONLY point) / distance','activation':'raw stage segment unsafe AND dot(V_raw,n)<0','stage_bases':{'START':'X_n -> X_n + dt/2 V_raw(X_n)','MIDPOINT':'X_n -> X_n + dt V_raw(X_mid)'},'wall_margin_m':1e-10,'wall_validation_tolerance_m':2e-14,'pair_swept_tolerance_m':1e-12,'velocity_projection':'V_used=V_raw-min(dot(V_raw,n),0)n only if active','Omega':'bitwise unchanged','position_correction_policy':'pending synthetic velocity-only evidence','no_new_physical_contact_threshold':True,'production_concentration':'UNSPECIFIED','flux':'TEST_ONLY; NOT EXPERIMENTAL CONCENTRATION','flow':'ENGINEERING_TRANSIENT_FIELD_ONLY','position_support':'WORKFLOW_V0_APPROXIMATION','scientific_production_ready':False}
(S/'contracts/WALL_SLIDING_CONTRACT.json').write_text(json.dumps(contract,indent=2)+'\n')
print('Protected BCFLUX files:',545,'remote inputs:',len(probe['inputs']))
