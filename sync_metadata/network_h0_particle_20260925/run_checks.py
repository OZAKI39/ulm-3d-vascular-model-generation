"""Run current regressions in a relocated temporary copy, preserving originals.

Only historical absolute paths in five JSON manifests are relocated. Generated
Python bytecode entries are omitted from the temporary source hash contract;
every original .py/.cpp source digest remains mandatory. No equations, test
assertions, tolerances, result arrays or source files are changed.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[2]
META=Path(__file__).resolve().parent
SV='formal_3D_flow_solver/FEM_SimVascular'
REPORT='particle_3d/reports/network_derived_flow_mb_validation_v1'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True,help='Directory for test logs and relocation receipts.')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    path_map=json.loads((META/'original_path_map.json').read_text())
    with tempfile.TemporaryDirectory(prefix='ulm-portable-checks-') as tmp:
        work=Path(tmp);copied=set()
        def copy(relative):
            source=ROOT/relative;target=work/relative
            if not source.exists():raise FileNotFoundError('Missing selected dependency: '+str(relative))
            if source.is_dir():
                for item in source.rglob('*'):
                    if item.is_file() and '__pycache__' not in item.parts:copy(str(item.relative_to(ROOT)))
            elif str(relative) not in copied:
                target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target);copied.add(str(relative))
        for relative in [
            'particle_3d/src','particle_3d/tests/network_flow_mb_validation_v1','particle_3d/tests/particle9a1_2mmps',
            REPORT,'particle_3d/reports/particle9a1_2mmps/reference',
            'vascular_network/network_1d0d','vascular_network/tests/network_1d0d','vascular_network/tests/network_h0','vascular_network/tests/conftest.py',
            'vascular_network/reports/a_network_1d0d_boundary_v1','vascular_network/reports/a_network_1d0d_boundary_v2_idealized',
            SV+'/src',SV+'/scripts/flow_2mmps',SV+'/scripts/sv13q/flow_parser.py',SV+'/tests/conftest.py',SV+'/tests/flow_2mmps',
            SV+'/frozen_reference',SV+'/upstream_stage_q_reference',SV+'/flow_cases/mean-2p0-mmps',
            SV+'/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/frozen_flow',SV+'/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/SV_MESH',
            SV+'/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/reports',SV+'/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/run/solver.xml',
            SV+'/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/policy.json',
            SV+'/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/run/PETSC_OPTIONS.txt','sonovue_size_distribution_v0']:
            copy(relative)
        for relative in path_map.values():copy(relative)

        prefixes=[
            ('/home/lzy/projects/ulm_flow_mean_2p0_mmps/'+SV+'/frozen_reference',SV+'/upstream_stage_q_reference'),
            ('/home/lzy/projects/ulm_flow_mean_2p0_mmps/'+SV,SV),
            ('/home/lzy/projects/ulm_particle_3d_particle0',''),
            ('/home/lzy/projects/ulm_3D_vascular','vascular_network'),
            ('/home/lzy/projects/formal_3D_flow_solver','formal_3D_flow_solver')]
        def remap(s):
            if s in path_map:return str(work/path_map[s])
            for old,new in prefixes:
                if s==old or s.startswith(old+'/'):return str(work/new)+s[len(old):]
            return s
        def relocate(value):
            if isinstance(value,str):return remap(value)
            if isinstance(value,list):return [relocate(x) for x in value]
            if isinstance(value,dict):return {remap(k):relocate(v) for k,v in value.items()}
            return value
        receipt=[]
        for relative in [
            'vascular_network/reports/a_network_1d0d_boundary_v1/data/protected_input_hashes.json',
            'vascular_network/reports/a_network_1d0d_boundary_v1/data/a_network_provenance.json',
            'vascular_network/reports/a_network_1d0d_boundary_v2_idealized/data/3d_case_config_diff.json',
            REPORT+'/data/old_new_flow_contract.json',REPORT+'/data/particle_science_protection_manifest.json']:
            source=ROOT/relative;target=work/relative;value=json.loads(source.read_text());cache_entries=[]
            if relative.endswith('particle_science_protection_manifest.json'):
                cache_entries=[k for k in value['files'] if Path(k).suffix in {'.pyc','.pyo'}]
                value['files']={k:v for k,v in value['files'].items() if k not in cache_entries}
                assert len(value['files'])>=114
            result=relocate(value);target.write_text(json.dumps(result,indent=2)+'\n')
            receipt.append(dict(file=relative,original_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                                temporary_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),omitted_generated_bytecode=len(cache_entries)))
        (a.output/'relocation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        env=dict(os.environ)
        env.update(PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VTK_SMP_MAX_THREADS='1',
                   SONOVUE_ROOT=str(work/'sonovue_size_distribution_v0'),
                   PYTHONPATH=os.pathsep.join(str(work/x) for x in ['particle_3d/src',SV+'/src','vascular_network']))
        groups=[('particle',work,['particle_3d/tests/network_flow_mb_validation_v1','particle_3d/tests/particle9a1_2mmps']),
                ('network',work/'vascular_network',['tests/network_1d0d','tests/network_h0']),
                ('flow',work/SV,['tests/flow_2mmps'])]
        results=[]
        for name,cwd,selectors in groups:
            if name=='flow':
                shutil.rmtree(work/SV/'frozen_reference')
                shutil.copytree(work/SV/'upstream_stage_q_reference',work/SV/'frozen_reference')
            log=a.output/(name+'.log');xml=a.output/(name+'.xml')
            cmd=[sys.executable,'-m','pytest','-q','-p','no:cacheprovider','--junitxml='+str(xml),*selectors]
            with log.open('w') as out:r=subprocess.run(cmd,cwd=cwd,env=env,stdout=out,stderr=subprocess.STDOUT)
            tests=ET.parse(xml).getroot().find('testsuite')
            results.append(dict(group=name,returncode=r.returncode,**{k:int(tests.attrib[k]) for k in ['tests','failures','errors','skipped']}))
            print(name,results[-1],flush=True)
        summary=dict(groups=results,temporary_copy=True,science_source_files_modified=False,
                     tests=sum(r['tests'] for r in results),failures=sum(r['failures']+r['errors'] for r in results),
                     skipped=sum(r['skipped'] for r in results))
        (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        if any(r['returncode'] for r in results):raise SystemExit(1)

if __name__=='__main__':main()
