"""Seal completed evidence; no numerical reruns or production mutations."""
from pathlib import Path
import json,hashlib,re,datetime

R=Path('/workspace/y_clot/project')
T='lps_strength_energy_decoupled_softening_feasibility_v1'
RUN='20261013T181824Z'
O=R/'verification'/T/RUN
F=R/'results'/T/RUN
A=R/'analysis/softening_feasibility_v1'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False))

status=read(O/'FINAL_STATUS.json')
assert 'PENDING' not in json.dumps(status)
required=['REVIEW_LPS_STRENGTH_ENERGY_DECOUPLED_SOFTENING_V1.md','FINAL_STATUS.json','DESIGN_FREEZE.json','SOURCE_AND_LIBRARY_IDENTITY.json','YANG_2019_EQUATION_AND_MODEL_AUDIT.md','BASELINE_LPS_ENERGY_FORCE_CONSISTENCY.json','SOFTENING_FREE_ENERGY_DEFINITION.md','DAMAGE_EVOLUTION_DEFINITION.md','CANDIDATE_SOFTENING_IMPLEMENTATION.py','INTACT_ELASTIC_LIMIT_AUDIT.json','FIXED_DAMAGE_FORCE_GRADIENT_AUDIT.json','LOADING_UNLOADING_RELOADING_AUDIT.json','IRREVERSIBLE_DISSIPATION_AUDIT.json','STRENGTH_ENERGY_DECOUPLING_AUDIT.json','ENERGY_TARGET_FEASIBILITY.json','IMPLEMENTATION_CHECKS.json','RUN_RESOURCES.json','REPRODUCE.txt','source_additions.diff']
assert all((O/p).is_file() for p in required)
report=(O/required[0]).read_text()
assert len(re.findall(r'^\*\*\d+\.',report,re.M))==35
assert sha(O/'provenance/Yang_2019_user_fulltext.pdf')==read(O/'PAPER_SOURCE.json')['sha256']
assert read(O/'PAPER_SOURCE.json')['user_specified_path_verified_identical']
for name,value in read(O/'DESIGN_FREEZE.json')['sources'].items():
    assert sha(A/name)==sha(O/'source_snapshot'/name)==value
provenance=read(O/'FIGURE_DATA_PROVENANCE.json')
assert sha(A/'figures.py')==provenance['script_sha256']
for name,value in provenance['inputs'].items():assert sha(O/name)==value
visual=read(O/'VISUAL_REVIEW.json')
for fig in provenance['figures']:
    v=next(z for z in visual['figures'] if z['name']==fig['name'])
    for ext in ['png','svg']:
        assert sha(F/(fig['name']+'.'+ext))==fig[ext+'_sha256']==v[ext+'_sha256']
assert len(provenance['figures'])==4 and (F/'OPEN_RESULTS.html').is_file()

excluded={'PAYLOAD_MANIFEST.json','PAYLOAD_SHA256SUMS.txt','MANIFEST_SHA256.txt','WSL_DELIVERY_VERIFICATION.json','SERVER_RECEIPT.json'}
def files():
    return sorted(p for root in [O,F,A] for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in excluded)
resources=read(O/'RUN_RESOURCES.json')
inventory=[p for p in files() if p!=O/'RUN_RESOURCES.json']
resources['file_volume']={'payload_bytes_excluding_resource_record_and_seals':sum(p.stat().st_size for p in inventory),'payload_file_count_excluding_resource_record_and_seals':len(inventory),'path_state_NPZ_bytes':{p.name:p.stat().st_size for p in O.glob('*.npz')},'final_figure_bytes':sum(p.stat().st_size for p in F.iterdir() if p.is_file()),'scope':'Own isolated code, evidence, source snapshots, failure history, paper provenance, WSL backup/diff, four PNG/SVG and HTML; no old large scientific runs.'}
write(O/'RUN_RESOURCES.json',resources)
records={str(p.relative_to(R)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in files()}
manifest={'status':'SEALED','run_id':RUN,'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'server_root':str(R),'path_base':'relative to flat project root; same layout in WSL clot_clearance root','excludes':sorted(excluded)+['__pycache__'],'exclusion_reason':'Manifest/checksum self-reference and later transfer receipts; receipts carry manifest hash. Python bytecode is reproducible cache.','payload_file_count':len(records),'payload_bytes':sum(v['bytes'] for v in records.values()),'files':records}
write(O/'PAYLOAD_MANIFEST.json',manifest)
(O/'PAYLOAD_SHA256SUMS.txt').write_text(''.join(v['sha256']+'  '+k+'\n' for k,v in records.items()))
(O/'MANIFEST_SHA256.txt').write_text(sha(O/'PAYLOAD_MANIFEST.json')+'  PAYLOAD_MANIFEST.json\n')
print(json.dumps({k:manifest[k] for k in ['status','payload_file_count','payload_bytes']})+' manifest_sha256='+sha(O/'PAYLOAD_MANIFEST.json'))
