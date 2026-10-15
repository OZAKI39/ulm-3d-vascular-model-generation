"""Archive a setup-only failure; preserve physical design and refreeze syntax.
Executed once before any accepted dynamic step. Not a numerical retry.
"""
from common import *
import numpy as np,shutil,difflib,time
f=O/'setup_failure_1';src=Path(__file__).parent
freeze=json.loads((f/'TIMESTEP_DESIGN_FREEZE.json').read_text())
assert json.loads((f/'dt_0p5us/DRIVER_ERROR.json').read_text())['last_accepted']==0
assert len((f/'dt_0p5us/native_scalars.csv').read_text().splitlines())==2
old=(src/'start_managed.py').read_text().replace("+C['run']+'_formal'","+C['run']").replace('managed_formal_run.sh','managed_run.sh')
assert hashlib.sha256(old.encode()).hexdigest()==freeze['diagnostic_source_sha256']['start_managed.py']
(f/'start_managed.py').write_text(old)
for n in ['driver.py','start_managed.py']:
 assert sha(f/n)==freeze['diagnostic_source_sha256'][n]
 (f/(n+'.diff')).write_text(''.join(difflib.unified_diff((f/n).read_text().splitlines(True),(src/n).read_text().splitlines(True),fromfile='failed/'+n,tofile='corrected/'+n)))
tr=freeze['loading']['ramp_s'];A=freeze['loading']['amplitude_m'];t=np.arange(201)*.25e-6
u=A*.5*(1-np.cos(np.pi*np.minimum(t/tr,1)));u2=(t<tr)*A*.5*(1-np.cos(np.pi*t/tr))+(t>=tr)*A
assert np.max(abs(u-u2)) <= 16*np.finfo(float).eps*A
record={'reason':'Locked LAMMPS treats min() as special reduction, not a binary scalar function','accepted_dynamic_steps':0,'initial_state_recorded':True,'newton_or_LPS_algorithm_changed':False,'original_freeze_sha256':sha(f/'TIMESTEP_DESIGN_FREEZE.json'),'function_equivalence_max_abs_m':float(abs(u-u2).max()),'native_instances_in_failed_attempt':1,'continued_dynamic_case_budget':2,'archive_tooling_note':'Two archive precondition assertions were corrected before this script: last_accepted was 0 (initial state), not -1; quote replacement in old managed-script reconstruction was corrected. Neither invoked LAMMPS.'}
dump(f/'SETUP_FAILURE_AND_CORRECTION.json',record)
q=O/'dt_0p5us';q.mkdir();shutil.copy2(O/'historical_static/body.data',q/'body.data');shutil.copy2(O/'historical_static/T0.npz',q/'initial_T0.npz');dump(q/'case.json',freeze['formal_cases'][0])
freeze['implementation_revision']={'epoch':time.time(),'physical_design_unchanged':True,'reason':record['reason'],'failed_attempt_archive':str(f)}
freeze['diagnostic_source_sha256']={q.name:sha(q) for q in src.iterdir() if q.is_file()}
dump(O/'TIMESTEP_DESIGN_FREEZE.json',freeze);print(record)
