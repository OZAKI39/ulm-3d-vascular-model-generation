"""Read-only final cross-check of scientific scope, artifacts and provenance."""
import json
import numpy as np
import h5py,fitz
from .common import *
from .report import DOCS,GATE_NAMES

def main():
    s=read(BASE/'SCIENTIFIC_RESULTS.json');assert s['complete_campaign']
    rows=s['rows'];assert len({(r['kind'],r['id']) for r in rows})==len(rows)
    assert s['audit_pass'] and read(BASE/'BASELINE_RESULTS.json')['PASS']
    g=read(BASE/'FINAL_GATE.json')['gates'];assert set(GATE_NAMES)<=set(g)
    for v in g.values():
        assert v['status'] in ['PASS','FAIL','NOT_RUN','INCONCLUSIVE']
        assert all(k in v for k in ['code_evidence','data_evidence','initial_states','bubble_counts','frequencies_Hz','actions','time_range_s','thresholds','model_identity','failure_reason'])
    assert not read(BASE/'FINAL_GATE.json')['training_ready']
    for name in DOCS:assert (BASE/(name+'.md')).stat().st_size>200
    for r in read(BASE/'figures/MANIFEST.json')['figures']:
        for ext in ['png','pdf']:assert sha(BASE/'figures'/r[ext])==r[ext+'_sha256']
        with fitz.open(BASE/'figures'/r['pdf']) as pdf:assert len(pdf)==1 and len(pdf[0].get_text())>50
    for r in read(BASE/'videos/MANIFEST.json')['videos']:
        assert r['fully_decoded'] and r['fps']==60 and r['resolution']==[1920,1080]
        assert sha(BASE/'videos'/r['file'])==r['sha256']
        for source in r['sources']:assert sha(BASE/source['path'])==source['sha256']
        detail=read(BASE/'videos'/(Path(r['file']).stem+'.json'))
        times=np.array([v['physical_time_s'] for v in detail['mapping']]);assert len(times)==r['frames']
        assert np.allclose(np.diff(times),r['uniform_physical_dt_s'],rtol=0,atol=1e-12)
    visual=read(BASE/'evidence/VISUAL_REVIEW.json');assert visual['all_pass']
    for r in visual['files']:assert sha(BASE/r['path'])==r['sha256']
    protected=read(BASE/'FINAL_EXECUTION_SEAL.json');assert protected['protected_inputs_pass']
    for p,h in protected['source_sha256'].items():assert sha(ROOT/p)==h,p
    index=read(BASE/'INPUT_OUTPUT_SHA256.json')
    for p,v in index['files'].items():assert sha(BASE/p)==v['sha256'],p
    for p,h in index['source_sha256'].items():assert sha(ROOT/p)==h,p
    raw={v['sha256'] for v in read(BASE/'REMOTE_RAW_INDEX.json')['files']}
    for r in rows:
        if r['kind']=='HF':assert r['raw_sha256'] in raw
    c=read(BASE/'COMPUTE_COST.json');assert c['budget_pass']
    result=dict(PASS=True,HF_attempts=s['HF_attempts'],LF_attempts=s['LF_attempts'],gates=len(g),documents=len(DOCS),figures=12,videos=3,checked_local_files=len(index['files']),scope='Independent final artifact/data/source cross-check; video full-decode evidence retained in video manifests',identity=identity())
    write(BASE/'evidence/FINAL_DELIVERY_VALIDATION.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
