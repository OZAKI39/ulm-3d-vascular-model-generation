"""Freeze V2 dependencies without mutating any historical scientific file."""
from common import *
import shutil,subprocess,platform,sys

def main():
    (O/'provenance').mkdir(parents=True,exist_ok=True)
    previous=json.loads((V2_O/'SOURCE_FREEZE.json').read_text())
    checks={}
    for name,want in previous['implementation_sha256'].items():
        actual=sha(R/name);assert actual==want,name;checks[name]=actual
    for name,want in previous['input_sha256'].items():
        assert sha(V2_O/name)==want,name
    for name in ['BIPHASIC_FIELD.npz','BIPHASIC_MATERIAL_SPEC.json','GEOMETRY_FREEZE.json','FROZEN_SHEAR_DRIVE.json','SHEAR_DRIVE.npz','MULTIFRAGMENT_OBJECTIVE.json','provenance/reference.npz','provenance/reference.native','provenance/body.data']:
        dest=O/name;dest.parent.mkdir(parents=True,exist_ok=True)
        if not dest.exists():shutil.copy2(V2_O/name,dest)
        assert sha(dest)==sha(V2_O/name)
    dump(O/'V2_DEPENDENCY_FREEZE.json',dict(source_freeze_sha256=sha(V2_O/'SOURCE_FREEZE.json'),implementation_sha256=checks,input_sha256=previous['input_sha256'],formal_cases=1,new_uniform=False,V2_remains_read_only=True))
    dump(O/'ENVIRONMENT.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,native_library_sha256=sha(LIB),current_git_sha=subprocess.check_output(['git','-c',f'safe.directory={R}','rev-parse','HEAD'],cwd=R,text=True).strip()))
    (O/'provenance/git_before.diff').write_bytes(subprocess.check_output(['git','-c',f'safe.directory={R}','diff'],cwd=R))
    (O/'provenance/git_before_status.txt').write_bytes(subprocess.check_output(['git','-c',f'safe.directory={R}','status','--short'],cwd=R))
    print('INPUTS_PRESERVED',sha(O/'BIPHASIC_FIELD.npz'))

if __name__=='__main__':main()
