from common import *
import shutil,re,subprocess,time,difflib,resource
assert not (O/'PREPARED.json').exists()
for p in [O,F,O/'provenance',O/'source_identity',O/'build',O/'preflight',O/'full',O/'full/checkpoints',O/'full/rejected',O/'full/keyframes']:p.mkdir(parents=True,exist_ok=True)
manifest=json.loads((BASELINE/'PAYLOAD_MANIFEST.json').read_text())
for p,z in manifest['files'].items():assert sha(R/p)==z['sha256'],p
protect=[R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so',SMALL/'build/liblammps.so',OLDLIB,Path('/workspace/y_clot/env.sh'),BASELINE/'PAYLOAD_MANIFEST.json']
dump(O/'PROTECTED_BASELINES.json',{str(p):sha(p) for p in protect})
shutil.copy2(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md')
for name in ['body.data','request.json','LOAD_PROVIDER.json','in.lammps']:shutil.copy2(HISTORY/name,O/'provenance'/name)
ref=dict(np.load(OLDV3/'C1_final/reference.npz'));X=ref['X'];b=ref['bonds'];bt=ref['bond_types'];n=int(ref['clot_count']);N=len(X)
assert n==1927 and N==3488 and len(b)==157637
assert np.array_equal(np.bincount(bt)[1:],[82270,22040,53327])
assert np.all(b[bt==1]<n) and np.all(b[bt==3]>=n) and np.all(np.sum(b[bt==2]<n,axis=1)==1)
deck=(HISTORY/'in.lammps').read_text();body=(HISTORY/'body.data').read_text()
records=np.array([[int(s) for s in line.split()[:4]] for line in body.split('\nBonds\n',1)[1].strip().splitlines() if line.strip()],dtype=np.int64)
assert np.array_equal(records[:,0],np.arange(1,len(b)+1)) and np.array_equal(records[:,1],bt) and np.array_equal(np.sort(records[:,2:4]-1,axis=1),b)
rr=np.linalg.norm(X[b[:,1]]-X[b[:,0]],axis=1);h=np.where(abs(rr-HORIZON)<=DX/2,(HORIZON+DX/2-rr)/DX,1.)
m=np.bincount(b.ravel(),weights=np.repeat(rr*h*V,2),minlength=N)
assert np.max(abs(m-ref['reference_m'])/m)<2e-12
area=np.zeros(N);tangent=np.zeros((N,3))
for line in deck.splitlines():
 if line.startswith('set atom ') and ' d_area ' in line:
  s=line.split();i=int(s[2])-1;area[i]=float(s[4]);tangent[i]=[float(s[6]),float(s[8]),float(s[10])]
assert not area[n:].any() and np.count_nonzero(area)==453
free=np.zeros((N,3),bool);free[:n]=True
np.savez_compressed(O/'provenance/reference.npz',X=X,bonds=b,bond_types=bt,original_bond_ids=records[:,0],r0=rr,taper=h,m=m,V=np.full(N,V),clot_count=n,free_dof=free,initial_cut=np.zeros(len(b),bool),areas=area,tangent=tangent)
with (O/'provenance/reference.native').open('w') as f:
 f.write(f'LPSVAR_REFERENCE_V1\n{N} {len(b)} {V:.17g}\n')
 for i,v in enumerate(m):f.write(f'{i+1} {v:.17g}\n')
 for i,(a,z) in enumerate(b):f.write(f'{i+1} {a+1} {z+1} {bt[i]} {rr[i]:.17g} {h[i]:.17g} 0 0 0\n')
req=json.loads((HISTORY/'request.json').read_text());load=json.loads((HISTORY/'LOAD_PROVIDER.json').read_text())
gates=json.loads((BASELINE/'SMALL_CRACK_DESIGN_FREEZE.json').read_text())['gates']
dump(O/'V3_TYPED_BOND_AND_ANCHOR_AUDIT.json',dict(status='PASS',material_points=n,fixed_anchors=N-n,total_points=N,original_bonds=len(b),type_counts={str(k):int(np.sum(bt==k)) for k in [1,2,3]},original_body_sha256=sha(HISTORY/'body.data'),reference_source=str(OLDV3/'C1_final/reference.npz'),reference_source_sha256=sha(OLDV3/'C1_final/reference.npz'),weighted_volume_max_relative=float(np.max(abs(m-ref['reference_m'])/m)),reference_m_uses_all_three_types=True,original_atom_and_bond_ids_exact=True,anchor_mass_excluded=True,material_mass_kg=n*MASS,material_volume_m3=n*V,load_surface_points=453,CFD_580_point_map_not_used=True,initial_bond_identity_sha256=hashlib.sha256(records.astype('<i8').tobytes()).hexdigest()))
dump(O/'V3_MATERIAL_CARD.json',dict(material_card_id='V3_ENDPOINT_STAR_TYPED_EXPLORATORY_V1',G_Pa=G,K_Pa=K,rho_kg_m3=RHO,dx_m=DX,horizon_m=HORIZON,alpha=0,bulk=dict(s0=S0[1],sc=SC[1]),interface=dict(s0=S0[2],sc=SC[2],status='EXPLORATORY_UNCALIBRATED'),anchor_scaffold=dict(s0=1.,sc=1.2,meaning='Fixed-fixed scaffold only; zero motion makes softening tail inactive'),INTERFACE_SOFTENING_ASSUMPTION='EXPLORATORY_UNCALIBRATED',historical_interface_s0=req['s00_interface'],historical_interface_onset_relative_difference=float(S0[2]/req['s00_interface']-1),GC_EFFECTIVE='UNAVAILABLE',macro_strength='UNVALIDATED',critical_load='UNAVAILABLE'))
dump(O/'PREFLIGHT_GATES.json',gates)
# Recompile every translation unit affected by the changed derived-class layout.
# The generated style_bond.cpp instantiates creator<T>; its sizeof(T) must agree.
buildold=json.loads((SMALL/'BUILD_AUDIT.json').read_text());cmd0=buildold['compile_command'];cut=cmd0.index('-c');flags=cmd0[:cut];flags[1:1]=['-I'+str(A)]
sources={'bond_bpm_peri_variational.cpp':A/'bond_bpm_peri_variational.cpp','fix_lps_variational_dynamic.cpp':A/'fix_lps_variational_dynamic.cpp','force.cpp':NATIVE/'source/src/force.cpp','style_bond.cpp':NATIVE/'build/styles/style_bond.cpp'}
cpp_before={}
for name in ['bond_bpm_peri_variational.cpp','bond_bpm_peri_variational.h','fix_lps_variational_dynamic.cpp','fix_lps_variational_dynamic.h']:
 old=NATIVE/'source_identity'/name;new=A/name
 assert old.is_file();cpp_before[name]=sha(old);shutil.copy2(old,O/'source_identity'/('before_'+name));shutil.copy2(new,O/'source_identity'/name)
 (O/'source_identity'/(name+'.diff')).write_text(''.join(difflib.unified_diff(old.read_text().splitlines(True),new.read_text().splitlines(True),fromfile='baseline/'+name,tofile='typed/'+name)))
objects=[Path(p) for p in buildold['object_sha256']]
assert len(objects)==432,len(objects)
for p in objects:assert sha(p)==buildold['object_sha256'][str(p)]
start=time.perf_counter();commands=[]
with (O/'build/build.log').open('w') as log:
 for name,src in sources.items():
  obj=O/'build'/(name+'.o');cmd=flags+['-c',str(src),'-o',str(obj)];commands.append(cmd);subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
  matches=[j for j,p in enumerate(objects) if p.name==name+'.o'];assert len(matches)==1,(name,matches);objects[matches[0]]=obj
 rsp=O/'build/objects.rsp';rsp.write_text('\n'.join(str(p) for p in objects)+'\n')
 link=buildold['link_command'].copy();link[link.index('-o')+1]=str(O/'build/liblammps.so.0')
 link=[('@'+str(rsp)) if s.startswith('@') else s for s in link];commands.append(link);subprocess.run(link,cwd=buildold['commands_cwd'],stdout=log,stderr=subprocess.STDOUT,check=True)
 LIB.symlink_to('liblammps.so.0')
dump(O/'BUILD_AUDIT.json',dict(status='PASS',compiled_translation_units=list(sources),commands=commands,link_cwd=buildold['commands_cwd'],wall_s=time.perf_counter()-start,reused_objects=len(objects)-len(sources),object_sha256={str(p):sha(p) for p in objects},old_objects_only_read=True))
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',dict(status='PASS',upstream='8de817dd79bfe4525d5d39246a212d833e6dee07',baseline_library=str(SMALL/'build/liblammps.so'),baseline_library_sha256=sha(SMALL/'build/liblammps.so'),new_library=str(LIB),new_library_sha256=sha(LIB),cpp_before_sha256=cpp_before,cpp_after_sha256={name:sha(A/name) for name in cpp_before},unchanged_elastic_kernel=True,mixed_damage_tail_by_immutable_type=True,production_replaced=False))
dump(O/'PREPARED.json',dict(status='PASS',run_id=RUN,baseline_manifest_verified_files=len(manifest['files']),original_input_sha256={p.name:sha(p) for p in (O/'provenance').iterdir()},library_sha256=sha(LIB)))
print('PREPARED_AND_BUILT',RUN,sha(LIB),flush=True)
