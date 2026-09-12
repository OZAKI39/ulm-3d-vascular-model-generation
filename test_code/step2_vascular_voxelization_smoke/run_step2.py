#!/usr/bin/python3
"""One new geometry-only run, no automatic retries or refinement."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from prepare_port_contract import EXPECTED_SHA, write_json
from verify_artifacts import snapshot

def run_once():
    source=Path(__file__).resolve().parent;hc=source.parent.parent
    step1=Path('/home/lzy/projects/compre_output/step1/20260912_215759')
    contract=step1/'geometry_contract'
    stl=contract/'geometry/cfd_surface_axis_aligned_inlet_m.stl'
    old=Path('/home/lzy/projects/ulm_3D_vascular')
    run=Path('/home/lzy/projects/compre_output/step2')/datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    run.mkdir(parents=True,exist_ok=False)
    for d in ['logs','diagnostics','build']:(run/d).mkdir()
    print('RUN_DIR='+str(run),flush=True)
    codes={}
    try:
        before=''
        for args in [['status','--short'],['rev-parse','HEAD'],['diff','--stat']]:
            before+='git '+' '.join(args)+'\n'+subprocess.check_output(['git','-C',str(hc),*args],text=True)
        (run/'git_status_before.txt').write_text(before)
        digest=hashlib.sha256(stl.read_bytes()).hexdigest()
        write_json(run/'diagnostics/input_integrity.json',dict(path=str(stl),sha256=digest,passed=digest==EXPECTED_SHA))
        assert digest==EXPECTED_SHA,'STEP1_INPUT_INTEGRITY=FAIL; stopped before voxelization'
        assert subprocess.check_output(['git','-C',str(hc),'rev-parse','HEAD'],text=True).strip()=='5a410848bd5c57d5ae1c171112e78eab4a82e650','HemoCell HEAD changed'
        for key,root in [('hemocell',hc),('step1',step1),('old_project',old)]:
            write_json(run/'diagnostics'/f'{key}_before_manifest.json',snapshot(root))
        tracked=subprocess.check_output(['git','-C',str(hc),'ls-files','-z']).decode().split('\0')
        write_json(run/'diagnostics/hemocell_tracked_before.json',{x:hashlib.sha256((hc/x).read_bytes()).hexdigest() if (hc/x).is_file() else None for x in tracked if x})
        write_json(run/'diagnostics/step1_hashes_before.json',{str(p.relative_to(step1)):hashlib.sha256(p.read_bytes()).hexdigest() for p in step1.rglob('*') if p.is_file()})
        checks=[]
        for line in (contract/'provenance/SHA256SUMS').read_text().splitlines():
            sha,name=line.split(maxsplit=1);checks.append(dict(path=name,passed=hashlib.sha256((contract/name.lstrip('*')).read_bytes()).hexdigest()==sha))
        write_json(run/'diagnostics/step1_package_verification.json',checks)
        assert all(c['passed'] for c in checks),'Frozen package SHA mismatch'
        prov=json.loads((contract/'provenance/provenance.json').read_text());bad=[]
        for item in prov['source_absolute_paths_and_sha256']:
            path=Path(item['absolute_path']);sha=item.get('sha256')
            if sha and (not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=sha):bad.append(str(path))
        write_json(run/'diagnostics/provenance_review.json',dict(source_evidence_count=len(prov['source_absolute_paths_and_sha256']),current_source_hash_mismatches=bad,unknowns=prov['unknowns']))
        assert not bad,'Step1 provenance source content changed'
        archive=hc/'build/libhemocell.a'
        write_json(run/'diagnostics/archive_link_provenance.json',dict(archive=str(archive),sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),size=archive.stat().st_size,read_only=True))
        env=dict(os.environ,PATH='/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(run/'build'))
        def stage(name,cmd,key):
            started=time.time()
            with (run/'logs'/(name+'.log')).open('w') as log:
                log.write('COMMAND '+repr(cmd)+'\n');log.flush()
                result=subprocess.run(cmd,cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT)
            codes[key]=result.returncode
            if key=='SMOKE_RC':codes['smoke_wall_seconds']=time.time()-started
            write_json(run/'logs/return_codes.json',codes)
            print(key+'='+str(result.returncode),flush=True)
            assert result.returncode==0,f'{name} failed; see log; no automatic retry'
        py=['/usr/bin/python3','-B']
        stage('prepare_ports',py+[str(source/'prepare_port_contract.py'),str(contract),str(run),str(old/'configs/cfd_flow.yaml')],'PREPARE_RC')
        stage('configure',['/usr/bin/cmake','-S',str(source),'-B',str(run/'build'),'-DCMAKE_CXX_COMPILER=/usr/bin/mpicxx','-DCMAKE_MAKE_PROGRAM=/usr/bin/make','-DCMAKE_BUILD_TYPE=Release'],'CONFIGURE_RC')
        stage('build',['/usr/bin/cmake','--build',str(run/'build'),'--parallel','1'],'BUILD_RC')
        s=json.loads((run/'run_settings.json').read_text())
        stage('smoke_mpi1',['/usr/bin/time','-v','/usr/bin/mpirun','-np','1',str(run/'build/closed_voxelizer'),str(stl),str(s['ref_dir']),str(s['ref_dir_n']),str(run)],'SMOKE_RC')
        stage('port_mapping',py+[str(source/'map_ports_and_diagnose.py'),str(run)],'PORT_MAPPING_RC')
        stage('artifact_verification',py+[str(source/'verify_artifacts.py'),str(run)],'ARTIFACT_VERIFICATION_RC')
        stage('reports',py+[str(source/'write_reports.py'),str(run)],'REPORT_RC')
    except Exception as error:
        (run/'EXECUTION_STOPPED.md').write_text('# 本次新运行停止\n\nSTEP2_STATUS=PARTIAL\nSTEP2_AUTO_CHECK=FAIL\n\n'+str(error)+'\n\n未自动重跑、加密、修几何或进入Step3。\n')
        raise

if __name__=='__main__':run_once()
