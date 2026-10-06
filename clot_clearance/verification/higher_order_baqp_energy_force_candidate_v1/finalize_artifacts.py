"""Seal this round's outputs after independent file/format/integrity checks."""
import ast,csv,hashlib,json,re,subprocess
from pathlib import Path
from PIL import Image
from pypdf import PdfReader

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
RESULTS=ROOT/'results/higher_order_baqp_energy_force_candidate_v1'
CODE=ROOT/'pd_clot/pfpd_candidate/higher_order_rk'
NEWCODE=[CODE/n for n in ['variational_energy.py','variational_force.py','adjoint_accumulation.py','mechanics_candidate.py','validation_states.py']]
MANIFEST=OUT/'OUTPUT_MANIFEST.json'
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def read(n):return json.loads((OUT/n).read_text())
def rows(n):
    with (OUT/n).open() as f:return list(csv.DictReader(f))
def dump(n,v):(OUT/n).write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n')

def main():
    protected=read('PROTECTED_INPUT_HASHES.json')
    changed=[p for p,h in protected.items() if not (ROOT/p).exists() or sha(ROOT/p)!=h]
    assert not changed,changed
    before=read('STARTING_GIT_STATE.json')
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()
    status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True)
    assert head==before['HEAD'] and branch==before['branch']
    oldset={x.strip() for x in before['status_short'].splitlines()};newset={x.strip() for x in status.splitlines()}
    added=sorted(newset-oldset);lost=sorted(oldset-newset)
    assert not lost,lost
    allowed={'?? verification/higher_order_baqp_energy_force_candidate_v1/','?? results/higher_order_baqp_energy_force_candidate_v1/'}
    assert set(added)<=allowed,added
    required=['CANDIDATE_CONFIGURATION.json','C1_ENERGY_REGRESSION.json','C1_FORCE_REGRESSION.npz',
        'C1_FORCE_REGRESSION_SUMMARY.json','OBJECTIVITY_RESULTS.json','GLOBAL_BALANCE_RESULTS.csv',
        'ENERGY_FORCE_DIRECTIONAL_DERIVATIVES.csv','SELECTED_DOF_FORCE_ORACLE.csv','AFFINE_ENERGY_COMPARISON.json',
        'QUADRATIC_ENERGY_DECOMPOSITION.csv','PERFORMANCE_SUMMARY.json','REVIEW_HIGHER_ORDER_BAQP_ENERGY_FORCE_V1.md']
    assert all((OUT/n).is_file() for n in required)
    assert not (OUT/'HARD_STOP.json').exists()
    assert not (OUT/'OPTIONAL_MANUFACTURED_RECOVERY.json').exists()
    report=(OUT/'REVIEW_HIGHER_ORDER_BAQP_ENERGY_FORCE_V1.md').read_text()
    assert re.findall(r'^## (\d+)\. ',report,re.M)==[str(n) for n in range(1,21)]
    figures={}
    for name in ['c1_variational_force_regression','pd2_energy_force_directional_derivative',
                 'pd2_objectivity_and_global_balance','c1_pd2_energy_decomposition']:
        with Image.open(RESULTS/'figures'/(name+'.png')) as im:
            dpi=im.info['dpi'];size=im.size;assert all(abs(d-600)<.02 for d in dpi)
        doc=PdfReader(RESULTS/'figures'/(name+'.pdf'));assert len(doc.pages)==1
        page=doc.pages[0];assert not list(page.images)
        drawings=sum(op in {b'S',b's',b'f',b'F',b'f*',b'B',b'B*',b'b',b'b*'} for args,op in page.get_contents().operations)
        assert drawings>0
        text=page.extract_text();assert len(text)>100 and not re.search('[\u4e00-\u9fff]',text)
        assert not read('FIGURE_LAYOUT_QA.json')[name]['text_outside_canvas']
        figures[name]=dict(pixels=size,dpi=dpi,vector_PDF=True,embedded_rasters=0,
            vector_drawing_operations=drawings,English_text=True,visual_inspection='PASS: readable labels and no observed text/plot overlap')
    assert len(list((RESULTS/'figures').glob('*.png')))==4
    assert len(list((RESULTS/'figures').glob('*.pdf')))==4
    cfg=read('CANDIDATE_CONFIGURATION.json');cache=read('PD2_CACHE_IDENTITY.json')
    assert sha(ROOT/cfg['primary_PD2_cache'])==cfg['primary_PD2_cache_sha256']==cache['primary_cache_sha256']
    assert sha(OUT/'N7_PD2_REFERENCE_CACHE.npz')==cache['small_verification_cache']['cache_sha256']
    fd=rows('ENERGY_FORCE_DIRECTIONAL_DERIVATIVES.csv');selected=rows('SELECTED_DOF_FORCE_ORACLE.csv')
    assert len(fd)==18 and len({(r['state'],r['direction']) for r in fd})==6 and len(selected)==12
    for row in fd:
        num=(float(row['U_plus_J'])-float(row['U_minus_J']))/(2*float(row['eps_m']))
        assert abs(num-float(row['D_fd_N']))<1e-25
        err=abs(num-float(row['D_force_N']))/float(row['denominator_N'])
        assert err<=cfg['FD_all_steps_gross_gate']
        if row['middle_step']=='True':assert err<=cfg['FD_middle_gate']
    for row in selected:assert float(row['normalized_error'])<=cfg['selected_DOF_gate']
    tests=read('TARGETED_TESTS.json');assert tests['passed']==12 and tests['failed']==0 and not tests['full_project_suite_run']
    flags=read('FINAL_STATUS.json');assert flags['HIGHER_ORDER_VARIATIONAL_MECHANICS_CANDIDATE']=='VERIFIED'
    for value in flags.values():assert value not in ['FAIL','REGRESSION_FAILED','NOT_VERIFIED']
    calls=[];imports=[]
    for p in NEWCODE:
        tree=ast.parse(p.read_text())
        for node in ast.walk(tree):
            if isinstance(node,ast.ImportFrom):imports.append(node.module or '')
            if isinstance(node,ast.Call):
                name=node.func.attr if isinstance(node.func,ast.Attribute) else node.func.id if isinstance(node.func,ast.Name) else ''
                if name in {'inv','pinv','hessian','degraded_svk_split_response','build_weights','minimize'}:calls.append([p.name,node.lineno,name])
    assert not calls,calls
    assert all('baqp'!=module and 'at1' not in module and 'petsc' not in module for module in imports)
    pending={OUT/'FINAL_ARTIFACT_QA.json',MANIFEST};bad=[]
    for p in list(OUT.rglob('*.md'))+[RESULTS/'README.md']:
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if target.startswith(('http://','https://','#')):continue
            dest=(p.parent/target.split('#')[0]).resolve()
            if not dest.exists() and dest not in pending:bad.append([p.name,target])
    assert not bad,bad
    perf=read('PERFORMANCE_SUMMARY.json');assert perf['minimum_available_RAM_fraction']>=.2
    qa=dict(PASS=True,protected_files=len(protected),protected_changed=changed,production_BAQP_unchanged=True,
        old_PD2_oracle_code_and_outputs_unchanged=True,old_manifests_unchanged=True,active_context_unchanged=True,
        branch=branch,HEAD=head,existing_git_status_preserved=True,new_git_entries=added,
        required_outputs=required,report_numbered_sections=20,figures=figures,local_document_links='PASS',
        PD2_cache_identity='IDENTICAL',six_FD_checks_all_18_values_verified=True,selected_DOF_count=12,
        code_review=dict(no_dense_Jacobian=True,no_Hessian=True,no_old_Z_PD2_path=True,
            no_moment_solve_in_energy_or_force_call=True,reference_weights_frozen=True,
            same_constitutive_energy_and_stress=True,adjoint_work='O(B + directed_neighbors + N)',forbidden_calls=calls),
        targeted_tests=tests,optional_recovery='NOT_RUN',next_task_started=False,numerical_recomputation_during_finalization=False)
    dump('FINAL_ARTIFACT_QA.json',qa)
    ce=read('C1_ENERGY_REGRESSION.json');cf=read('C1_FORCE_REGRESSION_SUMMARY.json');di=read('DELTA_FTILDE_IDENTITY.json')
    obj=read('OBJECTIVITY_RESULTS.json');bal=rows('GLOBAL_BALANCE_RESULTS.csv');aff=read('AFFINE_ENERGY_COMPARISON.json')
    fds=read('ENERGY_FORCE_FD_SUMMARY.json');quad=rows('QUADRATIC_ENERGY_DECOMPOSITION.csv')
    peak=max(perf['main_peak_RSS_bytes'],perf['peak_process_tree_RSS_bytes'])
    lines=[f'branch: {branch}',f'HEAD: {head}',
        'geometry: primary N19 = 6859 nodes / 690020 bonds; small N7 = 343 nodes / 21356 bonds',
        'h: 62.5 um','delta: 251.25 um','delta/h: 4.02','PD2 weight-cache identity: IDENTICAL',
        f"C1 energy max relative difference: {max(v['relative_difference'] for v in ce.values()):.9e}",
        f"C1 force max / weighted RMS normalized difference: {max(v['max_normalized'] for v in cf.values()):.9e} / {max(v['RMS_normalized'] for v in cf.values()):.9e}",
        f"delta-Ftilde identity max normalized residual: {di['max_normalized_residual']:.9e}",
        f"rigid translation energy / max force: {obj['translation']['energy_J']:.9e} J / {obj['translation']['max_force_N']:.9e} N",
        f"rigid rotation energy / max force: {obj['rotation']['energy_J']:.9e} J / {obj['rotation']['max_force_N']:.9e} N",
        f"global force maximum normalized residual: {max(float(v['force_residual']) for v in bal):.9e}",
        f"global current-moment maximum normalized residual: {max(float(v['moment_residual']) for v in bal):.9e}",
        f"energy-force FD max normalized error, all 18 / middle: {max(float(v['normalized_error']) for v in fd):.9e} / {max(v['middle_normalized_error'] for v in fds):.9e}"]
    for r in fds:lines.append(f"FD {r['state']}, direction {r['direction']+1}: {r['status']}; middle error {r['middle_normalized_error']:.9e}")
    lines.append(f"selected-DOF force oracle maximum normalized error: {max(float(r['normalized_error']) for r in selected):.9e}")
    for name,r in aff.items():lines.append(f"{name} C1-PD2 relative energy difference: {r['relative_energy_difference']:.9e}")
    for r in quad:
        lines.append(r['field']+' U_C1 / U_PD2 / U_exactkin / U_continuum [J]: '+ ' / '.join(f'{float(r[k]):.12e}' for k in ['U_C1_J','U_PD2_J','U_exactkin_currentweights_J','U_continuum_J']))
    lines+=['optional recovery: NOT_RUN','HIGHER_ORDER_VARIATIONAL_MECHANICS_CANDIDATE: VERIFIED',
        'production BAQP: UNCHANGED','boundary re-integration: NOT_REINTEGRATED_THIS_ROUND','phase force: NOT_IMPLEMENTED',
        f"targeted tests: 12 passed in {tests['wall_seconds']:.2f} s",'CPU workers: 3; BLAS threads: 1',
        f"wall time: {perf['wall_seconds']:.3f} s (final implementation numerical audit)",
        f'peak RAM: {peak/2**30:.3f} GiB (process high-water / sampled tree maximum)',
        'GPU: GPU_NOT_USED_BY_DESIGN',f"REVIEW: {OUT/'REVIEW_HIGHER_ORDER_BAQP_ENERGY_FORCE_V1.md'}",f"figures: {RESULTS/'figures'}"]
    summary='\n'.join(lines)+'\n';(OUT/'FINAL_TERMINAL_SUMMARY.txt').write_text(summary)
    paths=NEWCODE+[p for folder in [OUT,RESULTS] for p in folder.rglob('*') if p.is_file() and p!=MANIFEST]
    manifest={str(p.relative_to(ROOT)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(paths)}
    dump('OUTPUT_MANIFEST.json',dict(task=cfg['task'],scope='Only five new source modules plus this round verification/results; old files protected separately',
        excludes='This output manifest itself',files=manifest))
    loaded=read('OUTPUT_MANIFEST.json')['files']
    assert len(loaded)==len(paths) and all(sha(ROOT/p)==v['sha256'] for p,v in loaded.items())
    print(summary,end='')

if __name__=='__main__':main()
