from common import *
import time,ast,shutil
def main():
 started=time.perf_counter();assert json.loads((O/'FINAL_STATUS.json').read_text())['status']=='COMPLETED_P0_P4_NO_NEW_V3'
 assert json.loads((O/'SMALL_SPECIMEN_VALIDATION.json').read_text())['status']=='COMPLETED'
 assert (O/'WSL_DOCUMENT_UPDATE.json').exists() and (F/'FIGURE_QA.json').exists()
 expected=json.loads((B/'PAYLOAD_MANIFEST.json').read_text());bad=[]
 for name,item in expected['files'].items():
  p=R/name
  if not p.is_file() or p.stat().st_size!=item['bytes'] or sha(p)!=item['sha256']:bad.append(name)
 assert not bad,bad
 for name,h in json.loads((O/'PROTECTED_INPUTS.json').read_text()).items():assert sha(name)==h,name
 topo=np.genfromtxt(O/'derived/all_accepted_topology.csv',delimiter=',',names=True);history=np.genfromtxt(B/'V3_NATIVE_PILOT_HISTORY.csv',delimiter=',',names=True)
 for key in ['step','time_s','detached_nodes','detached_components','largest_detached_nodes']:assert np.array_equal(topo[key],history[key]),key
 geometry=R/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/V3_MATERIAL_POINTS_AND_BONDS.npz';g=load_npz(geometry);r=reference()
 for k in ['X','bonds','bond_types']:assert np.array_equal(g[k],r[k])
 sources={}
 for p in sorted(A.glob('*.py')):
  tree=ast.parse(p.read_text());sources[p.name]=dict(sha256=sha(p),entrypoints=[n.name for n in tree.body if isinstance(n,(ast.ClassDef,ast.FunctionDef))])
 dump(O/'CURRENT_CALL_CHAIN.json',dict(files=sources,chain=['prepare: verify old manifest and preserve call-chain evidence','topology_audit + spatial_audit + graph_ops: accepted connectivity and local histories','recover_states: lossless accepted before-state recovery','mixed_mode: objective kinematics and full endpoint-star virtual energy','group_audit: actual 13-node relative energy/force diagnostics','benchmark -> bench_native -> frozen V3 Native/Contact classes -> existing native library','benchmark_analysis -> figures -> finalize -> WSL insertion -> seal'],unchanged_native_library=sha(LIB)))
 for name in ['benchmark.py','bench_native.py','common.py','graph_ops.py']:
  f=json.loads((O/'benchmark_corrected_v3/DESIGN_FREEZE.json').read_text());assert sha(A/name)==f['driver_sha256'][name],name
 dump(O/'FINAL_INTEGRITY_CHECK.json',dict(status='PASS',baseline_verified_files=len(expected['files']),baseline_verified_bytes=expected['total_bytes'],baseline_unchanged=True,protected_inputs_unchanged=True,archived_surface_matches_exact_V3_geometry=True,geometry_sha256=sha(geometry),original_native_library_sha256=sha(LIB),wall_s=time.perf_counter()-started))
 payload={}
 for folder in [A,O,F]:
  for p in sorted(folder.rglob('*')):
   if not p.is_file() or '__pycache__' in p.parts or p.name in ['PAYLOAD_MANIFEST.json','LOCAL_READBACK.json'] or p.name.endswith('.tmp'):continue
   payload[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p))
 dump(O/'PAYLOAD_MANIFEST.json',dict(status='SEALED',task=TASK,run=RUN,files=payload,total_files=len(payload),total_bytes=sum(x['bytes'] for x in payload.values()),scope='New analysis, verification (including failures and corrected continuation), figures, actual WSL document evidence. Existing native binary/raw baseline retained at original paths and referenced by SHA; not duplicated.'))
 print('SEALED',len(payload),'files',sum(x['bytes'] for x in payload.values()),'bytes',flush=True)
if __name__=='__main__':main()
