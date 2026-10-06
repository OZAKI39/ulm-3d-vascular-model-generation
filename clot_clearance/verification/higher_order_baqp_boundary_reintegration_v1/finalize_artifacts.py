"""Seal this round only; do not modify historical outputs or manifests."""
import csv,json,re,subprocess,difflib
from pathlib import Path
import numpy as np
from PIL import Image
from pypdf import PdfReader
from run_audit import OUT,ROOT,CFG,KINDS,load,save,sha,integrity

R=ROOT/'results/higher_order_baqp_boundary_reintegration_v1'
FIGURES=[
('pd2_boundary_reintegration_overview','Exact integration architecture','Preserved PD2 internal mechanics and the external, frozen boundary conversion. The diagram describes the implementation; it is not a fit to a traction layer.'),
('pd2_affine_boundary_consistency','Affine force, moment, work and solved fields','Both affine solves start from a nonzero perturbation. Errors and predeclared gates use the historical normalization; Psi+ is an unchanged diagnostic.'),
('pd2_nonaffine_boundary_localization','Localization and proven closed-center checks','N19 layer maxima/means and four N27 lateral-depth samples are different sampling sets. Only the 13h N27 center has a fully closed dependency graph. The dashed threshold applies to that closed center.'),
('c1_vs_pd2_boundary_load_support','C1 versus PD2 conversion support','Magnitude-weighted d50/d90/d95/d99 and cumulative load fractions. A zero PD2 d50 means at least half the magnitude lies at the physical surface. These are conversion-support measures, not physical BC thickness.'),
('c1_vs_pd2_boundary_field_quality','Near-boundary field quality','Accepted historical C1 target-state diagnostics versus newly solved PD2 nonlocal-manufactured states. Node and bond-midpoint references use matched sampling locations. Normalized error in gradients and unchanged Psi+ diagnostic; no C1 campaign was rerun.'),
('pd2_manufactured_recovery','Four required perturbed-start recoveries','Target and solved gradients are compared with analytic fields, with actual mechanical residuals. Classical-loaded shear is intentionally separate because its local target has a nonzero bulk-model residual.'),
('pd2_classical_vs_nonlocal_shear','Classical shear control','Maximum deviation from the analytic target at each Z plane for displacement, nodal F and nodal Psi+. The classical load retains e_model=q_bulk+b_local. The NL recovery curve exposes numerical resolution; tiny shifts are not all assigned physical meaning.')]

def main():
    differences=[]
    for name in ['mechanics_adapter.py','localized_conversion.py']:
        before=OUT/'source_before'/name;after=ROOT/'pd_clot/pfpd_candidate/nonlocal_boundary'/name
        differences.extend(difflib.unified_diff(before.read_text().splitlines(True),after.read_text().splitlines(True),fromfile=str(before.relative_to(ROOT)),tofile=str(after.relative_to(ROOT))))
    (OUT/'BOUNDARY_SOURCE_EXTENSIONS.patch').write_text(''.join(differences))
    checks=[]
    for name,title,caption in FIGURES:
        png=R/'figures'/(name+'.png');pdf=R/'figures'/(name+'.pdf')
        with Image.open(png) as im:dim=im.size;dpi=im.info.get('dpi')
        reader=PdfReader(pdf);text='\n'.join(p.extract_text() for p in reader.pages)
        images=sum(len(p.images) for p in reader.pages)
        checks.append(dict(figure=name,png_dimensions=dim,dpi=dpi,pdf_pages=len(reader.pages),pdf_embedded_images=images,pdf_text_characters=len(text),PASS=bool(dpi and all(abs(x-600)<.1 for x in dpi) and images==0 and len(text)>50)))
    assert all(a['PASS'] for a in checks)
    doc=['# PD2 intact boundary reintegration figures\n',
         'Seven figures from recorded numerical data: six required figures and the classical-shear control. English labels, 600 dpi PNG and vector PDF.\n',
         '[Scientific review](../../verification/higher_order_baqp_boundary_reintegration_v1/REVIEW_HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1.md) · [Data dictionary](../../verification/higher_order_baqp_boundary_reintegration_v1/DATA_DICTIONARY.md)\n']
    for j,(name,title,caption) in enumerate(FIGURES,1):doc.append(f'## {j}. {title}\n\n{caption}\n\n[PNG](figures/{name}.png) · [PDF](figures/{name}.pdf)\n')
    (R/'README.md').write_text('\n'.join(doc))
    required=['PD2_BOUNDARY_CONFIGURATION.json','PD2_FORCE_STENCIL_CLOSURE.json','PD2_BULK_ORACLE.csv','PD2_NORMAL_AFFINE.json','PD2_SHEAR_AFFINE.json','PD2_NONAFFINE_AXIAL.json','PD2_NONAFFINE_SHEAR.json','PD2_CLASSICAL_SHEAR.json','PD2_NODEWISE_BOUNDARY_DECOMPOSITION.npz','PD2_BOUNDARY_SUPPORT.csv','C1_PD2_BOUNDARY_COMPARISON.csv','PD2_FORCE_MOMENT_WORK.csv','PD2_EQUILIBRIUM_RECOVERY.csv']
    assert all((OUT/n).is_file() for n in required)
    assert '12 passed' in (OUT/'TARGETED_TESTS.log').read_text()
    assert len(re.findall(r'^## \d+\.',(OUT/'REVIEW_HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1.md').read_text(),re.M))==20
    integ=integrity();assert integ['PASS'];save('HISTORICAL_INTEGRITY.json',integ)
    schema={}
    for p in sorted(OUT.glob('*.npz')):
        with np.load(p,allow_pickle=False) as a:schema[p.name]={k:dict(shape=a[k].shape,dtype=str(a[k].dtype)) for k in a.files}
    save('ARRAY_SCHEMA.json',schema)
    badlinks=[]
    for p in [OUT/'REVIEW_HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1.md',R/'README.md']:
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if not link.startswith(('https://','http://','#')) and not (p.parent/link.split('#')[0]).exists():badlinks.append(dict(file=str(p),target=link))
    assert not badlinks,badlinks
    save('FIGURE_QA.json',checks)
    (OUT/'GIT_STATUS_AFTER.txt').write_text(subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True))
    save('FINAL_ARTIFACT_QA.json',dict(PASS=True,required_outputs=len(required),figures=len(checks),targeted_tests=12,protected_files=integ['protected_files'],old_protected_changed=0,broken_local_links=badlinks,production_modified=False,accepted_PD2_modified=False,old_manifests_modified=False,continued_science_after_task=False))
    code=['pd_clot/pfpd_candidate/nonlocal_boundary/'+n for n in ['mechanics_backend.py','pd2_bulk_reference.py','mechanics_adapter.py','localized_conversion.py']]
    paths=[p for d in [OUT,R] for p in d.rglob('*') if p.is_file() and p.name!='OUTPUT_MANIFEST.json' and '__pycache__' not in p.parts and '.pytest_cache' not in p.parts]+[ROOT/p for p in code]
    save('OUTPUT_MANIFEST.json',dict(task=CFG['task'],scope='New round outputs plus explicitly extended/new boundary candidate modules; historical manifests unchanged',files={str(p.relative_to(ROOT)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(paths)}))
    print('SEALED',len(paths),'files; seven figures; twelve tests; protected inputs unchanged')

if __name__=='__main__':main()
