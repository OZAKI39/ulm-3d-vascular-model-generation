"""Run 115 portable legacy checks plus 35 new and 15 legacy inlet checks."""
from pathlib import Path
import argparse,json,os,subprocess,sys,time,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[2]

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',VTK_SMP_MAX_THREADS='1',PYTHONPATH=str(ROOT/'particle_3d/src'),SONOVUE_ROOT=str(ROOT/'sonovue_size_distribution_v0'))
 start=time.perf_counter()
 with (out/'portable_driver.log').open('w') as f:
  code=subprocess.run([sys.executable,str(ROOT/'sync_metadata/network_h0_particle_20260925/run_checks.py'),'--output',str(out/'portable')],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 if code:raise SystemExit(code)
 with (out/'p9a4_and_legacy.log').open('w') as f:
  code=subprocess.run([sys.executable,'-m','pytest','-q','-p','no:cacheprovider','--junitxml='+str(out/'p9a4_and_legacy.xml'),'particle_3d/tests/particle9a4_population_inlet','particle_3d/tests/particle9a2_inlet','particle_3d/tests/particle8_2a/test_open_inlet_not_solid_wall.py'],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 old=json.loads((out/'portable/summary.json').read_text());x=ET.parse(out/'p9a4_and_legacy.xml').getroot().find('testsuite')
 result=dict(portable=old,new_and_legacy={k:int(x.get(k)) for k in ['tests','failures','errors','skipped']},tests=old['tests']+int(x.get('tests')),failures=old['failures']+int(x.get('failures'))+int(x.get('errors')),skipped=old['skipped']+int(x.get('skipped')),runtime_s=time.perf_counter()-start)
 (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
 if code or result['failures']:raise SystemExit(1)
if __name__=='__main__':main()
