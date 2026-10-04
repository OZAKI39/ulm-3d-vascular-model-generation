"""Verify retained payload, references, syntax and imports without simulations."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from audit_inventory import AUDIT, PROJECT, ROOT, dump_json, extract_references, hash_file


def verify_payload():
    rows = json.loads((AUDIT/'CLASSIFIED_INVENTORY.json').read_text())
    kept = [r for r in rows if r['candidate_category'] != 'DELETE_SAFE']
    removed = [r for r in rows if r['candidate_category'] == 'DELETE_SAFE']
    missing = [r['path'] for r in kept if not os.path.lexists(ROOT/r['path'])]
    assert not missing, missing[:10]
    with ThreadPoolExecutor(max_workers=4) as pool:
        checks = list(pool.map(hash_file, kept))
    changed = [r['path'] for r,h in zip(kept,checks) if r['hash'] != h]
    unexpected = [r['path'] for r in removed if os.path.lexists(ROOT/r['path'])]
    papers = []
    for path in Path('/home/lzy/projects').glob('*.pdf'):
        if path.name in ('2604.09215v1.pdf', '1-s2.0-S0045782510001283-main.pdf', '1-s2.0-S2666496826000567-main.pdf', 'An adaptive dynamic relaxation method for quasi-static simulations using the peridynamic theory.pdf'):
            papers.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'size_bytes':path.stat().st_size})
    paper = next(r for r in papers if r['path'].endswith('/2604.09215v1.pdf'))
    assert paper['sha256']=='a77d29c88bd45bbd1cba78f1e823c09a9157b3be436b5eb6b3855c97a260a157'
    dump_json('PRESERVATION_CHECK.json', dict(passed=not changed and not missing and not unexpected,
                                           kept_entries_verified=len(kept), unchanged_hashes=len(kept)-len(changed),
                                           changed=changed,missing=missing,deleted_entries_reappeared=unexpected,
                                           protected_external_papers=papers,hash_workers=4))
    assert not changed and not unexpected, (changed[:10],unexpected[:10])
    refs, unresolved, expressions=extract_references(rows)
    before=json.loads((AUDIT/'REFERENCE_SCAN_BEFORE.json').read_text())
    broken=[r for r in unresolved if r['kind']=='markdown_link']
    lost=[r for r in before['resolved'] if r['within_repository'] and not os.path.lexists(ROOT/r['target'])]
    dump_json('REFERENCE_INTEGRITY.json',dict(broken_referenced_path_count=len(broken)+len(lost),
              preexisting_literal_markdown_broken_count=sum(r['kind']=='markdown_link' for r in before['unresolved']),
              new_broken_references=len(lost),broken=broken,lost=lost,resolved_reference_count=len(refs),
              unresolved_contextual_mentions=unresolved,expression_count=len(expressions),
              scope='Literal local Markdown links and resolved artifact paths in first-party docs and STATUS/SUMMARY JSON. Historical manifest listings, pseudocode, DOI/URLs and external-checkout labels are not live artifact dependencies.',
              restoration='Exact canonical originals restored using saved snapshot SHA256; no report links rewritten.'))
    assert not broken and not lost
    print('Preserved hashes:',len(kept),'broken references:',len(broken)+len(lost),flush=True)


def syntax_imports():
    temp=Path(json.loads((AUDIT/'VALIDATION_TEMP.json').read_text())['path'])
    python='/home/lzy/projects/computation_examples/_shared/python_environment/bin/python'
    env=os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE='1',NUMBA_CACHE_DIR=str(temp/'numba'),MPLCONFIGDIR=str(temp/'mpl'),
               OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1',PYTHONPATH=str(PROJECT))
    modules=['pd_clot.mechanics','pd_clot.streaming','pd_clot.regularization.equilibrium_adr',
             'pd_clot.pfpd_candidate.baqp','pd_clot.pfpd_candidate.kernel','pd_clot.pfpd_candidate.normalization',
             'pd_clot.pfpd_candidate.svk_spectral_split','pd_clot.pfpd_candidate.at1_evolution',
             'pd_clot.pfpd_candidate.petsc_vi_candidate']
    commands=[('compileall',[python,'-X','pycache_prefix='+str(temp/'compileall'),'-m','compileall','-q','pd_clot','scripts','tests'],env),
              ('python_imports',[python,'-B','-c','import importlib,json; names='+repr(modules)+'; [importlib.import_module(n) for n in names]; print(json.dumps(names))'],env)]
    petsc_env=env.copy()
    petsc_env['PYTHONPATH']=str(PROJECT/'verification/pfpd_at1_petsc_vi_solver_v1/external/python_deps')+os.pathsep+str(PROJECT)
    commands.append(('petsc_runtime',['/home/lzy/.cache/codex-dolfinx-env/bin/python','-B','-c',
        'from petsc4py import PETSc; import pd_clot.pfpd_candidate.petsc_vi_candidate; import numba, vtk; print(PETSc.Sys.getVersion()); print(numba.__version__); print(vtk.vtkVersion.GetVTKVersion())'],petsc_env))
    julia_env=env.copy();julia_env['JULIA_DEPOT_PATH']=str(PROJECT/'verification/pfpd_candidate_v1/external/julia_depot')
    commands.append(('julia_runtime',[str(PROJECT/'verification/pfpd_candidate_v1/external/julia-1.10.12/bin/julia'),
                     '--startup-file=no','--compiled-modules=no','-e','print(VERSION)'],julia_env))
    reports=[]
    for name,cmd,e in commands:
        r=subprocess.run(cmd,cwd=PROJECT,env=e,capture_output=True,text=True)
        (AUDIT/(name.upper()+'.log')).write_text(r.stdout+r.stderr)
        reports.append(dict(check=name,command=cmd,returncode=r.returncode,log=name.upper()+'.log'))
        print(name,r.returncode,flush=True)
        if r.returncode:
            dump_json('SYNTAX_IMPORT_SMOKE.json',dict(passed=False,checks=reports));raise RuntimeError(name)
    dump_json('SYNTAX_IMPORT_SMOKE.json',dict(passed=True,checks=reports))


if __name__=='__main__':
    verify_payload()
    if '--imports' in sys.argv:
        syntax_imports()
