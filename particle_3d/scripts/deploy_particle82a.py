#!/usr/bin/env python3
"""Snapshot only committed scientific source into a new isolated remote run."""
from pathlib import Path
import hashlib, json, shlex, subprocess, tempfile

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/'particle_3d/reports/particle8_2a'
REMOTE='/workspace/particle8_2a_20260922'
SSH=['ssh','-p','4159','-i',str(Path.home()/'.ssh/vast_step3b_ed25519'),
     '-o','IdentitiesOnly=yes','-o','BatchMode=yes','-o','ConnectTimeout=15','root@50.115.148.16']


def main():
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    remote=REMOTE+'/'+commit+'/repo';base=REMOTE+'/'+commit
    paths=['particle_3d/src','particle_3d/tests','particle_3d/scripts','particle_3d/contracts','particle_3d/pyproject.toml']
    changed=subprocess.check_output(['git','diff','HEAD','--name-only','--','*82a*','particle_3d/tests/particle8_2a'],cwd=ROOT,text=True)
    if changed:raise ValueError('Commit P8.2A source before deployment')
    old=json.loads((ROOT/'particle_3d/reports/particle8_2/REMOTE_RUN_MANIFEST.json').read_text())
    subprocess.run(SSH+['mkdir -p '+shlex.quote(remote)],check=True)
    archive=Path(tempfile.gettempdir())/('particle82a_'+commit+'.tar')
    with archive.open('wb') as f:subprocess.run(['git','archive',commit,*paths],cwd=ROOT,stdout=f,check=True)
    with archive.open('rb') as f:subprocess.run(SSH+['tar -xf - -C '+shlex.quote(remote)],stdin=f,check=True)
    command='test -e '+shlex.quote(remote+'/formal_3D_flow_solver')+' || cp -al '+shlex.quote(old['remote_repo']+'/formal_3D_flow_solver')+' '+shlex.quote(remote+'/formal_3D_flow_solver')
    subprocess.run(SSH+[command],check=True)
    host=json.loads((REPORT/'SERVER_COMPUTE_PROVENANCE.json').read_text())
    host.update(source_git_commit=commit,frozen_input_sha256=old['frozen_input_sha256'])
    source_files=subprocess.check_output(['git','ls-files',*paths],cwd=ROOT,text=True).splitlines()
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    sono=Path('/home/lzy/projects/sonovue_size_distribution_v0')
    host['external_scientific_input_sha256']={str(sono/p):digest(sono/p) for p in
        ['SHA256SUMS']+[line.split('  ',1)[1] for line in (sono/'SHA256SUMS').read_text().splitlines()]}
    host['source_sha256']={p:digest(ROOT/p) for p in source_files
        if 'particle82_figures.py' not in p and 'particle82_visuals.py' not in p and 'finalize_particle82_report.py' not in p}
    for file,value in [(base+'/provenance.json',json.dumps(host,indent=2)),(remote+'/SOURCE_COMMIT',commit+'\n')]:
        subprocess.run(SSH+['cat > '+shlex.quote(file)],input=value,text=True,check=True)
    verify='from pathlib import Path;import json,hashlib;p=Path('+repr(remote)+');h=json.load(open('+repr(base+'/provenance.json')+'));assert all(hashlib.sha256((p/k).read_bytes()).hexdigest()==v for typ in ["source_sha256","frozen_input_sha256","external_scientific_input_sha256"] for k,v in h[typ].items());print("SOURCE_AND_FROZEN_INPUTS_VERIFIED")'
    subprocess.run(SSH+['python3 -c '+shlex.quote(verify)],check=True)
    manifest=dict(remote_base=base,remote_repo=remote,source_commit=commit,remote_hostname=host['hostname'],
        archive_sha256=digest(archive),frozen_input_sha256=host['frozen_input_sha256'],source_sha256=host['source_sha256'])
    (REPORT/'DEPLOYMENT.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(base)


if __name__=='__main__':main()
