import subprocess
from fem_sync_support import ROOT, read


def test_branch_descends_from_recorded_latest_main_without_rewriting_it():
    repo=ROOT.parents[1];meta=read('sync_metadata/main_base.json')
    git=lambda *a:subprocess.check_output(['git',*a],cwd=repo,text=True).strip()
    assert git('rev-parse','origin/main')==meta['base_commit']
    assert git('merge-base',meta['base_commit'],'HEAD')==meta['base_commit']
    assert git('branch','--show-current')==meta['branch'] and meta['branch']!='main'
    for name in git('diff','--name-only',meta['base_commit']).splitlines():
        assert name=='README.md' or name.startswith('formal_3D_flow_solver/FEM_SimVascular/')


def test_existing_main_files_are_byte_identical_except_readme_append():
    repo=ROOT.parents[1];base=read('sync_metadata/main_base.json')['base_commit']
    names=subprocess.check_output(['git','ls-tree','-r','--name-only',base],cwd=repo,text=True).splitlines()
    for name in names:
        original=subprocess.check_output(['git','show',base+':'+name],cwd=repo)
        if name=='README.md':assert (repo/name).read_bytes().startswith(original)
        else:assert (repo/name).read_bytes()==original,name
