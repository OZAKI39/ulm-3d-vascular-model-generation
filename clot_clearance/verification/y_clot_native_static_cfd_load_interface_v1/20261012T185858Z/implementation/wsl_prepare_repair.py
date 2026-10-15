from pathlib import Path
import json,hashlib
C=json.loads(Path(__file__).with_name('config.json').read_text());P=Path(C['local'])/'analysis/static_cfd_v1';O=Path(C['local'])/'verification'/C['task']/C['run'];before=(P/'interface.py').read_text();after=before.replace("assert sum(d0['nonidentity_permutation_counts'])>0","assert sum(d1['nonidentity_permutation_counts'])>0")
(O/'interface_before_repair.py').write_text(before);(P/'interface.py').write_text(after)
# Complete audits from already saved full native states. No force evaluation repeated.
block=after.split('if rank==0:\n',1)[1].split('\nn.close()',1)[0]
head="""from common import *
import numpy as np
from types import SimpleNamespace
p=O/'native_interface';r=np.load(O/'REFERENCE_AND_LOAD.npz');X=r['X'];f=r['force_N'];n=SimpleNamespace(ref=r)
base=np.load(p/'reference_unloaded.npz');loaded0=np.load(p/'reference_loaded.npz');loaded1=np.load(p/'reference_loaded_after_migration.npz')
db=json.loads((p/'reference_unloaded.json').read_text());d0=json.loads((p/'reference_loaded.json').read_text());d1=json.loads((p/'reference_loaded_after_migration.json').read_text())
delta=loaded0['force']-base['force'];delta1=loaded1['force']-base['force'];eps=float(np.linalg.norm(delta-f,axis=1).max());eps1=float(np.linalg.norm(delta1-f,axis=1).max());moved=np.flatnonzero(loaded0['owner']!=loaded1['owner'])+1
assert len(moved)>0 and eps<1e-18 and eps1<1e-18 and sum(d1['nonidentity_permutation_counts'])>0
assert np.array_equal(base['x'],X) and np.array_equal(loaded0['x'],X) and np.array_equal(loaded1['x'],X)
"""
(P/'finish_interface_audit.py').write_text(head+'\n'.join(l[1:] if l.startswith(' ') else l for l in block.splitlines())+'\n')
launch=(P/'launch.py').read_text().replace("with (O/'CAMPAIGN_LAUNCHED.json').open('x')", "with (O/'CAMPAIGN_REPAIR_LAUNCHED.json').open('x')")
launch=launch.replace("('interface',['mpiexec','--allow-run-as-root','-n','2',sys.executable,str(Path(__file__).with_name('interface.py'))],fr['budget']['interface_wall_seconds']),", "('finish_interface_audit',[sys.executable,str(Path(__file__).with_name('finish_interface_audit.py'))],30),")
launch=launch.replace("O/'CAMPAIGN_PERFORMANCE.json'","O/'REPAIR_AND_SOLVE_PERFORMANCE.json'")
(P/'launch_repair1.py').write_text(launch)
s=(P/'start_managed.py').read_text().replace("not (O/'CAMPAIGN_LAUNCHED.json').exists()","not (O/'CAMPAIGN_REPAIR_LAUNCHED.json').exists()").replace("+'_formal'","+'_repair1'").replace('managed_formal_run.sh','managed_repair1_run.sh').replace('analysis/static_cfd_v1/launch.py','analysis/static_cfd_v1/launch_repair1.py').replace('MANAGED_START.txt','MANAGED_REPAIR_START.txt').replace('/supervisor.log','/supervisor_repair1.log')
(P/'start_repair1.py').write_text(s)
import difflib
(O/'IMPLEMENTATION_REPAIR_1.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='interface.py.before',tofile='interface.py.after')))
(O/'IMPLEMENTATION_REPAIR_1.json').write_text(json.dumps({'reason':'First interface assertion demanded a nonidentity sort at reference_loaded although reference positions already sorted and unchanged. Saved initial permutation has [1800,1684] and loaded after true ownership migration [658,2830]. Applied-force errors are zero in both loaded states. Check after migration instead.','native_calls_repeated':0,'recovery':'Finish audit from full precision saved native states, then first and only formal static solve','physical_parameters_changed':False,'gates_changed':False,'freeze_sha256':hashlib.sha256((O/'STATIC_LOAD_DESIGN_FREEZE.json').read_bytes()).hexdigest(),'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in P.glob('*.py')}},indent=2))
