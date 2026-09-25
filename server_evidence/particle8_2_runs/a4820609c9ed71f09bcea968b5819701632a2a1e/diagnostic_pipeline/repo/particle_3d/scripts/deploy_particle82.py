#!/usr/bin/env python3
"""Deploy exact committed source and manifest-required immutable inputs by SSH."""
from pathlib import Path
import argparse,hashlib,json,shlex,subprocess,tempfile

REPO=Path(__file__).resolve().parents[2]
OUT=REPO/'particle_3d/reports/particle8_2'
SSH=['ssh','-p','4159','-i',str(Path.home()/'.ssh/vast_step3b_ed25519'),'-o','IdentitiesOnly=yes',
     '-o','BatchMode=yes','-o','ServerAliveInterval=30','root@50.115.148.16']


def main(run_id):
    if not run_id.replace('_','').replace('-','').isalnum():raise ValueError('Invalid run ID')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
    changed=subprocess.check_output(['git','diff','HEAD','--name-only'],cwd=REPO,text=True).splitlines()
    if changed:raise ValueError('Deploy only from committed tracked source')
    base=f'/root/particle8_2_runs/{commit}/{run_id}';remote=base+'/repo'
    def ssh(cmd,**kwargs):return subprocess.run(SSH+[cmd],check=True,**kwargs)
    ssh('test ! -e '+shlex.quote(base)+' && mkdir -p '+shlex.quote(remote))
    fem=Path('formal_3D_flow_solver/FEM_SimVascular')
    extras=[str(fem/line.split('  ',1)[1]) for line in (REPO/fem/'frozen_reference/SHA256SUMS.txt').read_text().splitlines() if not line.split('  ',1)[1].startswith('frozen_reference/')]
    paths=['particle_3d/src','particle_3d/scripts','particle_3d/tests','particle_3d/contracts','particle_3d/pyproject.toml',str(fem/'frozen_reference'),*extras]
    archive=Path(tempfile.gettempdir())/f'particle82_{commit}.tar'
    with archive.open('wb') as f:subprocess.run(['git','archive',commit,*paths],cwd=REPO,stdout=f,check=True)
    with archive.open('rb') as f:ssh('tar -xf - -C '+shlex.quote(remote),stdin=f)
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    files=[p for p in (REPO/fem/'frozen_reference').rglob('*') if p.is_file()]+[REPO/p for p in extras]
    frozen={str(p.relative_to(REPO)):digest(p) for p in files}
    sources=subprocess.check_output(['git','ls-files','particle_3d/src','particle_3d/scripts','particle_3d/tests','particle_3d/contracts'],cwd=REPO,text=True).splitlines()
    source_sha={p:digest(REPO/p) for p in sources}
    host=json.loads((OUT/'REMOTE_HOST_PROVENANCE.json').read_text())
    host.update(source_git_commit=commit,frozen_input_sha256=frozen,source_sha256=source_sha)
    ssh('cat > '+shlex.quote(base+'/host_provenance.json'),input=json.dumps(host,indent=2),text=True)
    ssh('cat > '+shlex.quote(remote+'/SOURCE_COMMIT'),input=commit+'\n',text=True)
    verify="from pathlib import Path;import json,hashlib;r=Path("+repr(remote)+");h=json.load(open("+repr(base+'/host_provenance.json')+"));assert all(hashlib.sha256((r/p).read_bytes()).hexdigest()==v for k in ['frozen_input_sha256','source_sha256'] for p,v in h[k].items());print('EXACT_SOURCE_AND_FROZEN_HASH_MATCH')"
    ssh('python3 -c '+shlex.quote(verify))
    old_path=OUT/'REMOTE_RUN_MANIFEST.json';old=json.loads(old_path.read_text()) if old_path.exists() else {}
    old_deployments=old.get('deployments',[{k:v for k,v in old.items() if k!='frozen_input_sha256'}] if old else [])
    manifest=dict(remote_base=base,remote_repo=remote,source_commit=commit,remote_hostname=host['hostname'],
        remote_endpoint=host['ssh_endpoint'],snapshot_archive_sha256=digest(archive),frozen_input_sha256=frozen,source_sha256=source_sha,
        existing_remote_results='PRESERVED; each deployment is a new directory',deployments=old_deployments)
    old_deployments.append({k:v for k,v in manifest.items() if k not in ['deployments','frozen_input_sha256','source_sha256']})
    old_path.write_text(json.dumps(manifest,indent=2)+'\n')
    print(base)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run-id',required=True);main(p.parse_args().run_id)
