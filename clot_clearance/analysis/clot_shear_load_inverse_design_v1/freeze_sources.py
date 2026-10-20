from common import *
import subprocess,datetime,difflib
assert json.loads((O/'KEY_CHECKS.json').read_text())['status']=='PASS'
assert not (O/'SOURCE_FREEZE.json').exists()
for p,q in json.loads((O/'BASELINE_MANIFEST.json').read_text())['files'].items():assert sha(p)==q['sha256'],p
sources={str(p.relative_to(R)):sha(p) for p in A.glob('*.py')}
git=['git','-c','safe.directory='+str(R)]
head=subprocess.check_output(git+['rev-parse','HEAD'],cwd=R,text=True).strip()
status=subprocess.check_output(git+['status','--short','--untracked-files=no'],cwd=R,text=True)
checks=json.loads((O/'KEY_CHECKS.json').read_text())
dump(O/'SOURCE_FREEZE.json',dict(status='FROZEN_BEFORE_MAIN_RUNS',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),git_sha=head,existing_tracked_changes=status,
     native_lammps_version=checks['native_version'],native_library_sha256=sha(LIB),material_card='MATERIAL_AND_SPECIMEN_FREEZE.json',specimen_sha256=sha(O/'provenance/reference.npz'),
     input_sha256={n:sha(O/n) for n in ['BASELINE_MANIFEST.json','MATERIAL_AND_SPECIMEN_FREEZE.json','SHEAR_PARAMETERIZATION.json','SHEAR_FIELD_LIBRARY.json','CURRENT_SURFACE_SPEC.json','MULTIFRAGMENT_OBJECTIVE.json']},
     implementation_sha256=sources,
     old_entrypoints=[str(PA/'native.py'),str(PA/'run_shear_cases.py'),str(PA/'shear_loading.py')],
     new_entrypoints=['shear_field_library.py','current_surface.py','native.py','run_cases.py','multifragment_shear_objective.py'],
     change_scope='New isolated field/surface/external callback and diagnostics; original native VV, material kernel, contact, MPI and restart implementation unchanged. Only two main cases. No optimization launches by this freeze.'))
(O/'SOURCE_DIFF_AT_FREEZE.patch').write_text(''.join(''.join(difflib.unified_diff([],p.read_text().splitlines(True),fromfile='/dev/null',tofile='b/'+str(p.relative_to(R)))) for p in sorted(A.glob('*.py'))))
print('SOURCE_FROZEN',head,len(sources))
