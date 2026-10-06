"""Read-only evidence checks, then seal only this task's new artifacts."""
import ast
import hashlib
import json
import re
import subprocess
from pathlib import Path
import numpy as np
from PIL import Image
from pypdf import PdfReader

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
CODE=ROOT/'pd_clot/pfpd_candidate/higher_order_rk'
RESULTS=ROOT/'results/higher_order_baqp_rk_oracle_v1'
MANIFEST=OUT/'OUTPUT_MANIFEST.json'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def read(name):return json.loads((OUT/name).read_text())
def dump(name,value):(OUT/name).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')

def main():
    protected=read('PROTECTED_INPUT_HASHES.json')
    changed=[p for p,h in protected.items() if not (ROOT/p).is_file() or sha(ROOT/p)!=h]
    assert not changed,changed
    package=read('JULIA_PACKAGE_SOURCE_HASHES.json')
    package_changed=[p for p,h in package.items() if sha(Path(p))!=h]
    assert not package_changed,package_changed
    git=read('STARTING_GIT_STATE.json')
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()
    status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True)
    assert head==git['HEAD'] and branch==git['branch']
    before={s.strip() for s in git['status_short'].splitlines()}
    after={s.strip() for s in status.splitlines()}
    new=sorted(after-before);lost=sorted(before-after)
    allowed={'?? pd_clot/pfpd_candidate/higher_order_rk/','?? verification/higher_order_baqp_rk_oracle_v1/',
             '?? results/higher_order_baqp_rk_oracle_v1/'}
    assert not lost,lost
    assert set(new)<=allowed,new
    required=['PD2_CONFIGURATION.json','PD2_NODE_CONDITIONING.csv','PD2_POLYNOMIAL_REPRODUCTION.csv',
        'PD2_AXIAL_QUADRATIC.csv','PD2_SHEAR_QUADRATIC.csv','PD2_CROSS_TERM.csv','PD2_FTILDE_DIAGNOSTIC.csv',
        'PD2_PSI_PLUS_DIAGNOSTIC.csv','C1_PD2_SUMMARY.json','PERIDYNAMICSJL_CROSSCHECK.json','GMLS_TRIGGER_DECISION.md',
        'REVIEW_HIGHER_ORDER_BAQP_RK_ORACLE_V1.md','FINAL_STATUS.json','TARGETED_TESTS.json']
    assert all((OUT/p).is_file() for p in required)
    review=(OUT/'REVIEW_HIGHER_ORDER_BAQP_RK_ORACLE_V1.md').read_text()
    headings=re.findall(r'^## (\d+)\. ',review,re.M)
    assert headings==[str(n) for n in range(1,20)]
    figqa={}
    for name in ['c1_vs_pd2_boundary_bias','pd2_polynomial_reproduction','pd2_conditioning',
                 'pd2_ftilde_and_energy_bias','python_vs_peridynamicsjl_pd2']:
        with Image.open(RESULTS/'figures'/(name+'.png')) as im:
            dpi=im.info['dpi'];size=im.size
            assert all(abs(v-600)<.02 for v in dpi)
        doc=PdfReader(RESULTS/'figures'/(name+'.pdf'))
        assert len(doc.pages)==1
        page=doc.pages[0];images=list(page.images)
        drawings=sum(op in {b'S',b's',b'f',b'F',b'f*',b'B',b'B*',b'b',b'b*'}
                     for args,op in page.get_contents().operations)
        text=page.extract_text()
        assert not images and drawings>0,(name,len(images),drawings)
        assert not re.search('[\u4e00-\u9fff]',text)
        assert len(text)>100
        assert not read('FIGURE_LAYOUT_QA.json')[name]['text_outside_canvas']
        figqa[name]=dict(pixels=size,dpi=dpi,pdf_pages=1,pdf_embedded_rasters=0,
            pdf_vector_drawings=drawings,English_text=True,text_within_canvas=True,
            visual_inspection='PASS: readable labels, no observed text/plot overlap')
    provenance=read('CACHE_PROVENANCE.json');cachecheck={}
    for name,entry in provenance.items():
        file=OUT/(name+'_REFERENCE_CACHE.npz')
        assert sha(file)==entry['cache_sha256']
        assert sha(OUT/'PD2_CONFIGURATION.json')==entry['configuration_sha256']
        with np.load(file) as a:
            assert hashlib.sha256(a['g'].tobytes()).hexdigest()==entry['g_array_sha256']
        cachecheck[name]='UNCHANGED'
    papers=read('references/SOURCE_PROVENANCE.json')['used']
    for item in papers:assert sha(OUT/'references'/(item['key']+'.pdf'))==item['sha256']
    forbidden=[]
    for p in list(CODE.glob('*.py'))+[OUT/'run_audit.py']:
        tree=ast.parse(p.read_text())
        for node in ast.walk(tree):
            if isinstance(node,ast.Call):
                name=node.func.attr if isinstance(node.func,ast.Attribute) else node.func.id if isinstance(node.func,ast.Name) else ''
                if name in {'inv','pinv','evaluate','svk_split_second_pk','force','stress_divergence'}:
                    forbidden.append((str(p.relative_to(ROOT)),node.lineno,name))
    assert not forbidden,forbidden
    # Check ordinary Markdown local links. The two sealed products are created below.
    pending={MANIFEST,OUT/'FINAL_ARTIFACT_QA.json'}
    badlinks=[]
    for p in list(OUT.glob('*.md'))+[RESULTS/'README.md']:
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if target.startswith(('https://','http://','#')):continue
            dest=(p.parent/target.split('#')[0]).resolve()
            if not dest.exists() and dest not in pending:badlinks.append((p.name,target))
    assert not badlinks,badlinks
    tests=read('TARGETED_TESTS.json')
    assert tests['passed']==12 and tests['failed']==0 and not tests['full_project_suite_run']
    assert read('PERIDYNAMICSJL_CROSSCHECK.json')['status']=='PASS'
    report=dict(PASS=True,protected_files=len(protected),protected_changed=changed,
        installed_Julia_package_files=len(package),Julia_source_changed=package_changed,
        branch=branch,HEAD=head,existing_git_status_preserved=True,new_git_status_entries=new,
        old_manifests_modified=False,production_modified=False,active_context_modified=False,
        work_namespaces=[str(p.relative_to(ROOT)) for p in [CODE,OUT,RESULTS]],
        required_outputs_present=required,review_numbered_sections=19,figures=figqa,
        local_document_links='PASS',reference_cache_hashes=cachecheck,primary_source_PDF_hashes='PASS',
        forbidden_calls_in_primary_candidate=forbidden,
        targeted_tests=tests,scientific_recomputation_during_finalization=False,
        next_task_started=False)
    dump('FINAL_ARTIFACT_QA.json',report)
    s=read('C1_PD2_SUMMARY.json');c=read('CONDITIONING_SUMMARY.json')['PD2']
    poly=read('POLYNOMIAL_SUMMARY.json');a=read('PD2_AFFINE_REGRESSION.json')
    cpu=read('RESOURCE_USAGE.json');jcpu=read('JULIA_RESOURCE_USAGE.json')
    psi=[s[f][model][k] for f in s for model in ['PD2'] for k in ['Psi_Fbar','Psi_Ftilde']]
    text=f'''branch: {branch}
HEAD: {head}
N19: 19 x 19 x 19 = 6859 nodes; 690020 physical bonds
h: 62.5 um
delta: 251.25 um
delta/h: 4.02
kernel: (1-r/delta)^2; unchanged
PD2 basis dimension: 9
moment min rank: {c['min_rank']}
moment sigma_min minimum: {c['sigma_min_min']:.9e} m^3
condition median / P95 / maximum: {c['condition_median']:.9f} / {c['condition_P95']:.9f} / {c['condition_max']:.9f}
rank-deficient count: {c['rank_deficient_count']}
polynomial max / worst per-polynomial volume-weighted RMS: {max(v['max_scaled'] for v in poly.values()):.9e} / {max(v['weighted_RMS_scaled'] for v in poly.values()):.9e}
C1 axial max F bias / 0.004: {s['axial']['C1']['Fbar']['max_scaled']:.9e}
C1 shear max F bias / 0.004: {s['shear']['C1']['Fbar']['max_scaled']:.9e}
PD2 affine max F error / 0.004: {max(v['max_scaled'] for v in a.values()):.9e}
'''
    for field in s:
        v=s[field]['PD2']['Fbar']
        text+=f"PD2 {field} max / volume-weighted RMS F error / 0.004: {v['max_scaled']:.9e} / {v['weighted_RMS_scaled']:.9e}\n"
    text+=f"RK2_DIAG XZ: FAIL (expected control); max / 0.004 = {s['cross']['RK2_DIAG']['Fbar']['max_scaled']:.9e}\n"
    for field in s:
        v=s[field]['PD2']['Ftilde']
        text+=f"PD2 Ftilde {field} max / energy-weighted RMS / 0.004: {v['max_scaled']:.9e} / {v['weighted_RMS_scaled']:.9e}\n"
    text+=f'''Psi_plus C1 max relative bias (all fields, nodes/bonds): {max(s[f]['C1'][k]['max_scaled'] for f in s for k in ['Psi_Fbar','Psi_Ftilde']):.9e}
Psi_plus PD2 max relative bias (all fields, nodes/bonds): {max(v['max_scaled'] for v in psi):.9e}
Psi_plus PD2 worst field/sample weighted RMS relative bias: {max(v['weighted_RMS_scaled'] for v in psi):.9e}
Peridynamics.jl: PASS
GMLS_REQUIRED: NO
HIGHER_ORDER_KINEMATIC_CANDIDATE: VERIFIED
production force: NOT_IMPLEMENTED
production BAQP: UNCHANGED
targeted tests: 12 passed in {tests['wall_seconds']:.2f} s; full suite NOT_RUN
CPU workers: Python 8; BLAS 1; Julia 1
wall time: Python {cpu['wall_seconds']:.3f} s + Julia {jcpu['wall_seconds']:.3f} s = {cpu['wall_seconds']+jcpu['wall_seconds']:.3f} s (computation stages only)
peak RAM: {max(cpu['peak_process_tree_RSS_bytes'],jcpu['peak_process_tree_RSS_bytes'])/2**30:.3f} GiB (monitored process-tree RSS)
GPU: GPU_NOT_USED_BY_DESIGN
REVIEW: {OUT/'REVIEW_HIGHER_ORDER_BAQP_RK_ORACLE_V1.md'}
figures: {RESULTS/'figures'}
'''
    (OUT/'FINAL_TERMINAL_SUMMARY.txt').write_text(text)
    files=sorted(p for directory in [CODE,OUT,RESULTS] for p in directory.rglob('*') if p.is_file() and p!=MANIFEST)
    content={str(p.relative_to(ROOT)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in files}
    MANIFEST.write_text(json.dumps(dict(task='HIGHER_ORDER_BAQP_RK_ORACLE_V1',
        excludes='This output manifest only',files=content),indent=2)+'\n')
    saved=json.loads(MANIFEST.read_text())['files']
    assert len(saved)==len(files)
    assert all(sha(ROOT/p)==v['sha256'] and (ROOT/p).stat().st_size==v['bytes'] for p,v in saved.items())
    print(text,end='')

if __name__=='__main__':main()
