from sv13o_support import *
def test_stage_N_and_production_frozen():
 d=accepted('reference_freeze');n=json.loads((ROOT/'reports/sv1_3n/gpu_steady_candidate.json').read_text())
 assert d['GPU_baseline']==n and d['CPU_production']=='CPU_EARLY_STOP_PRODUCTION'
 assert d['PETSc']['tag']=='v3.25.5' and d['PETSc_build']['CUDA_version']=='13.2'
 assert hashlib.sha256((ROOT/'configs/sv1_3/policy.json').read_bytes()).hexdigest()==d['production_policy_sha256']
 assert d['scientific_equivalence']=='DEFERRED'
def test_source_change_is_only_final_output_condition():
 d=accepted('source_patch');assert d['changed_files']==['Code/Source/solver/main.cpp']
 changed=[n for n in d['before'] if d['before'][n]!=d['after'][n]];assert changed==d['changed_files']
 patch=(ROOT/d['patch']).read_text();assert hashlib.sha256(patch.encode()).hexdigest()==d['patch_sha256']
 assert 'save_vtu || (com_mod.saveVTK && reached_stop_time_step)' in patch
 assert all(not x.startswith(('+','-')) or x.startswith(('+++','---')) or any(s in x for s in ['if (save_vtu','//']) for x in patch.splitlines())
def test_scientific_XML_is_identical_for_all_vascular_runs():
 for p in (ROOT/'outputs/sv1_3o').glob('*/solver.xml'):
  if p.parent.name=='OFFICIAL_OUTPUT_STOP':continue
  scientific_xml_gate(ROOT/'outputs/sv1_3n/REAL_VASCULAR_GPU/solver.xml',p)

