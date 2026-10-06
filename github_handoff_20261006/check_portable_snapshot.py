"""Check the unchanged small suites in an isolated copy of selected files only.

No scientific simulations, package installs, Git mutations, or output rewrites.
The temporary copy remains local; it is not a backup of unsynchronized artifacts.
"""
from pathlib import Path
import ast, csv, hashlib, json, os, shutil, subprocess, sys, tempfile, time

R=Path(__file__).resolve().parents[1];OUT=R/'github_handoff_20261006'
inventory=list(csv.DictReader((OUT/'SCIENTIFIC_FILE_INVENTORY.csv').open()))
selected=[r['path'] for r in inventory if r['decision'] in ('INCLUDE','INHERITED')]
dest=Path(tempfile.mkdtemp(prefix='clot-pd2-handoff-check-'))
for rel in selected:
    src=R/rel;target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,target)
known=set(selected)
missing=[];parsed=0
# Resolve imports recursively in all project source, scripts and focused suites.
for rel in selected:
    if not rel.endswith('.py'):continue
    p=R/rel
    try:tree=ast.parse(p.read_text())
    except (UnicodeError,SyntaxError) as e:
        missing.append({'file':rel,'parse_error':str(e)});continue
    parsed+=1
    parts=Path(rel).relative_to('clot_clearance').with_suffix('').parts
    package=list(parts[:-1])
    for n in ast.walk(tree):
        if not isinstance(n,(ast.Import,ast.ImportFrom)):continue
        if isinstance(n,ast.Import):mods=[a.name for a in n.names]
        elif n.level:
            base=package[:len(package)-n.level+1];mods=['.'.join(base+([n.module] if n.module else []))]
        else:mods=[n.module or '']
        for mod in mods:
            if not mod.startswith(('pd_clot.','scripts.','verification.','vendor.')):continue
            q=Path('clot_clearance',*mod.split('.'));options=[str(q.with_suffix('.py')),str(q/'__init__.py')]
            exists=[x for x in options if (R/x).is_file()]
            if exists and not any(x in known for x in exists):missing.append({'file':rel,'import':mod,'existing_but_not_selected':exists})
env=os.environ.copy();env.update(PYTHONPATH='.',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
jobs=[
 ['verification/higher_order_baqp_energy_force_candidate_v1/tests/test_variational_candidate.py','-k','not test_01_cached_PD2_kinematics_preserved'],
 ['verification/higher_order_baqp_boundary_reintegration_v1/tests/test_boundary_reintegration.py','-k','not test_01_pd2_mechanics_cache_and_force_regression_unchanged and not test_04_nonaffine_nodewise_decomposition and not test_12_classical_shear_target_residual'],
]
records=[];logs=[]
for args in jobs:
    cmd=[sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider',*args];t=time.monotonic()
    p=subprocess.run(cmd,cwd=dest/'clot_clearance',env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    records.append({'command':cmd,'returncode':p.returncode,'wall_seconds':time.monotonic()-t,'output':p.stdout});logs.append('COMMAND: '+repr(cmd)+'\n'+p.stdout)
result={'snapshot':str(dest),'copied_files':len(selected),'new_equilibrium_solves':0,'source_parse_count':parsed,'local_import_missing_files':missing,'tests':records,'PASS':not missing and all(x['returncode']==0 for x in records),'scope':'20 unchanged existing tests; other tests require explicitly excluded large/server-only evidence'}
(OUT/'PORTABLE_TARGETED_TESTS.json').write_text(json.dumps(result,indent=2)+'\n')
(OUT/'PORTABLE_TARGETED_TESTS.log').write_text('\n'.join(logs))
print(json.dumps(result,indent=2));sys.exit(0 if result['PASS'] else 1)
