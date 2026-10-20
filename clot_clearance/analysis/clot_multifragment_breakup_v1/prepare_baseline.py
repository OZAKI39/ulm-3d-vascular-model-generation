from common import *
import shutil,datetime

assert not (O/'BASELINE_MANIFEST.json').exists(), 'Preserve existing freeze'
for folder in ['provenance','build','checks','cases','offline']:
    (O/folder).mkdir(parents=True,exist_ok=True)
prior=json.loads((HISTORY/'DELIVERY_MANIFEST.json').read_text())
index={q['path']:q for q in prior['files']}
used=[]
for filename in ['reference.npz','reference.native','body.data']:
    source=HISTORY/'provenance'/filename
    expected=index[str(source.relative_to(R))]['sha256']
    assert sha(source)==expected
    shutil.copy2(source,O/'provenance'/filename)
    used.append(dict(path=str(source),sha256=expected))
for name in ['uniform','heterogeneous']:
    for filename in ['SUMMARY.json','ANALYSIS.json','CRACK_NETWORK_AUDIT.json','history.csv','damage_events.jsonl.gz','keys/initial.npz','keys/peak_indent.npz','keys/withdrawal_end.npz','keys/final.npz']:
        p=HISTORY/'cases'/name/filename;digest=sha(p)
        assert digest==index[str(p.relative_to(R))]['sha256']
        used.append(dict(path=str(p),sha256=digest))
    summary=json.loads((HISTORY/'cases'/name/'SUMMARY.json').read_text())
    assert summary['status']=='COMPLETED'
    assert summary['last_row']['detached_nodes']==(18 if name=='uniform' else 14)
ref=load_npz(O/'provenance/reference.npz')
assert int(ref['clot_count'])==1792 and len(ref['X'])==2128 and len(ref['bonds'])==95666 and ref['initial_cut'].sum()==2736
build=json.loads((HISTORY/'BUILD_AUDIT.json').read_text());assert sha(OLDLIB)==build['library_sha256']
for p in [HA/'bond_bpm_peri_variational.cpp',HA/'bond_bpm_peri_variational.h',HA/'fix_lps_variational_dynamic.cpp',HA/'fix_lps_variational_dynamic.h',HA/'indenter.py',HA/'native.py',R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',BRACKETING/'FINAL_STATUS.json']:
    used.append(dict(path=str(p),sha256=sha(p)))
dump(O/'BASELINE_MANIFEST.json',dict(status='FROZEN',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),previous_delivery_manifest_sha256=sha(HISTORY/'DELIVERY_MANIFEST.json'),baseline_library_sha256=sha(OLDLIB),used_files=used,material_nodes=1792,fixed_anchors=336,original_bonds=95666,original_precuts=2736,uniform_material=True,material=dict(G_Pa=G,K_Pa=K,rho_kg_m3=RHO,dx_m=DX,horizon_m=HORIZON,s0_by_type=S0,sc_by_type=SC),scope='Read-only baseline verification; no historical rerun, no new pre-cuts or material changes'))
print('BASELINE_FROZEN',sha(O/'BASELINE_MANIFEST.json'))
