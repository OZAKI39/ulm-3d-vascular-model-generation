"""Verify the published handoff branch via GitHub API and an independent clone.

Read-only network actions. This records the commit being verified, never predicts
the SHA of a later commit containing this verification receipt. No CFD.
"""
import base64
import hashlib
import json
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
META = json.loads((ROOT/'sync_metadata/main_base.json').read_text())
BRANCH = META['branch']
PREFIX = 'formal_3D_flow_solver/FEM_SimVascular/'
GIT = ['git','-c','credential.helper=','-c','credential.helper=!gh auth git-credential']


def git(*args):
    return subprocess.check_output([*GIT,*args],cwd=REPO,text=True).strip()


def main():
    git('fetch','origin')
    head = git('rev-parse','HEAD')
    remote = git('rev-parse','origin/'+BRANCH)
    assert head == remote
    remote_main = git('ls-remote','--heads','origin','main').split()[0]
    assert remote_main == META['base_commit']
    files = []
    for name in ['PARTICLE_HANDOFF.md','FROZEN_FEM_BASELINE.md',
                 'frozen_reference/flow/flow_field_manifest.json',
                 'frozen_reference/mesh_manifest.json','frozen_reference/boundary_manifest.json']:
        endpoint = f'repos/{META["repository"]}/contents/{PREFIX}{name}?ref={quote(BRANCH,safe="")}'
        d = json.loads(subprocess.check_output(['gh','api',endpoint],text=True))
        assert d['encoding']=='base64'
        payload = base64.b64decode(d['content'])
        assert payload == (ROOT/name).read_bytes(),name
        files.append(dict(path=name,sha256=hashlib.sha256(payload).hexdigest(),github_blob_sha=d['sha'],bytes=len(payload)))
    with tempfile.TemporaryDirectory(prefix='fem_handoff_github_fresh_') as tmp:
        checkout = Path(tmp)/'repo'
        subprocess.run([*GIT,'clone','--depth','1','--single-branch','--branch',BRANCH,META['remote_url'],str(checkout)],check=True)
        clone_head = subprocess.check_output(['git','rev-parse','HEAD'],cwd=checkout,text=True).strip()
        assert clone_head == head
        fem = checkout/PREFIX
        validation = subprocess.check_output([sys.executable,'-I','-B',str(fem/'scripts/fem_freeze_sync/validate_frozen.py')],cwd=tmp,text=True)
        inventory = subprocess.check_output([sys.executable,'-I','-B',str(fem/'scripts/fem_freeze_sync/verify_manifest.py')],cwd=tmp,text=True)
        assert 'PASS:' in validation and 'PASS:' in inventory
        assert not subprocess.check_output(['git','status','--porcelain'],cwd=checkout)
    receipt = dict(status='PASS',repository=META['repository'],branch=BRANCH,verified_commit=head,
        local_HEAD=head,remote_HEAD=remote,remote_main=remote_main,main_unchanged=True,
        verified_at=datetime.now(timezone.utc).isoformat(),github_content_reread=files,
        independent_fresh_clone=True,fresh_clone_validation=validation.strip(),fresh_clone_manifest=inventory.strip(),
        fresh_clone_worktree_clean=True,no_CFD=True,scope='Verification anchors the published snapshot. The following evidence-only commit must be pushed and its final HEAD rechecked separately.')
    (ROOT/'sync_metadata/remote_verification.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(f'PASS: remote/local {head}; five GitHub files byte-identical; independent remote clone fresh-read and full inventory verified; main unchanged.')


if __name__=='__main__':
    main()
