#!/usr/bin/env python3
"""Run permanent checks and bind their result to source and rendered artifacts."""
from pathlib import Path
import sys, subprocess, xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle8_replay import REPO,write,digest
from particle_3d.particle8_full3d_data import OUTPUT,verify_preservation


def source_files():
    paths=[]
    for pattern in ['PARTICLE8_FULL3D_README.md','scripts/*particle8_full3d.py',
                    'src/particle_3d/particle8_full3d_*.py','tests/particle8_full3d/*.py']:
        paths.extend((REPO/'particle_3d').glob(pattern))
    return {str(p.relative_to(REPO)):digest(p) for p in sorted(paths)}


def artifact_files():
    return {str(p.relative_to(OUTPUT)):digest(p) for folder in ['data','frames','mp4','gif','keyframes','storyboard']
            for p in sorted((OUTPUT/folder).glob('*')) if p.is_file() and p.name not in
            ['test_run.json','visual_inspection.json','frame_audit.csv','active_particles.csv']}


def main():
    before=source_files();artifacts=artifact_files();results=[]
    for suite in ['particle8_full3d','particle8']:
        name='full3d_tests' if suite.endswith('full3d') else 'upstream_particle8'
        xml=OUTPUT/'logs'/f'{name}.xml';log=OUTPUT/'logs'/f'{name}.log'
        cmd=[sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider',f'--junitxml={xml}',f'particle_3d/tests/{suite}']
        with log.open('w') as stream:result=subprocess.run(cmd,cwd=REPO,stdout=stream,stderr=subprocess.STDOUT)
        suites=ET.parse(xml).getroot().findall('testsuite')
        counts={key:sum(int(s.get(key,0)) for s in suites) for key in ['tests','failures','errors','skipped']}
        results.append(dict(suite=suite,command=cmd,exit_code=result.returncode,**counts,
            junit=str(xml.relative_to(OUTPUT)),log=str(log.relative_to(OUTPUT)),junit_sha256=digest(xml),log_sha256=digest(log)))
        print(suite,counts,flush=True)
    unchanged=before==source_files() and artifacts==artifact_files()
    passed=unchanged and all(r['exit_code']==0 and r['tests']>0 and not (r['failures']+r['errors']+r['skipped']) for r in results)
    write(OUTPUT/'data/test_run.json',dict(suites=results,all_pass=passed,inputs_unchanged_during_tests=unchanged,
        source_sha256=before,artifact_sha256=artifacts,upstream_preserved_files=verify_preservation(OUTPUT)))
    return 0 if passed else 1


if __name__=='__main__':sys.exit(main())
