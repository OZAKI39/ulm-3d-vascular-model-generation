from common import *
import numpy as np,subprocess,sys,shutil,time
from scipy.spatial import cKDTree
O.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True);assert not (O/'PREFLIGHT.json').exists()
start=time.perf_counter();commands=[]
for a in [['pwd'],['git','status','--short'],['git','branch','--show-current'],['git','rev-parse','HEAD'],['ps','-eo','pid,comm,pcpu,pmem,etimes'],['df','-h','.'],['lscpu'],['free','-b'],['nvidia-smi','--query-gpu=name,memory.used,utilization.gpu','--format=csv'],['ldd',str(B/'liblammps.so')]]:
 r=subprocess.run(a,capture_output=True,text=True);commands.append({'argv':a,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
cgroup={p:Path(p).read_text().strip() for p in ['/sys/fs/cgroup/cpu/cpu.cfs_quota_us','/sys/fs/cgroup/cpu/cpu.cfs_period_us','/sys/fs/cgroup/memory/memory.limit_in_bytes','/sys/fs/cgroup/memory/memory.usage_in_bytes']}
dump(O/'PREFLIGHT.json',{'commands':commands,'python':sys.executable,'python_version':sys.version,'cgroup':cgroup,'new_native_instances':0})
meta=json.loads((MAP/'CLOT_NODAL_TRACTION_METADATA.json').read_text());sources={str(R/p):s for p,s in meta['source_SHA'].items()};sources[str(MAP/meta['file'])]=meta['force_sha256'];sources[str(MAP/'CFD_TO_LAMMPS_MAPPING_MATRIX.npz')]=meta['map_sha256']
sourcecheck={p:{'expected':s,'actual':sha(p),'bytes':Path(p).stat().st_size} for p,s in sources.items()};assert all(a['expected']==a['actual'] for a in sourcecheck.values())
a=np.load(MAP/meta['file']);g=np.load(GDIR/'V3_MATERIAL_POINTS_AND_BONDS.npz');s=np.load(CFD/'refined/CLOT_SURFACE_TRACTIONS.npz');m=np.load(MAP/'CFD_TO_LAMMPS_MAPPING_MATRIX.npz');x0=a['reference_point_m'];X=g['X'];b=g['bonds'];types=g['bond_types'];nc=int(g['clot_count']);assert nc==NC and len(X)==N and len(b)==157637
np.testing.assert_array_equal(a['material_atom_id'],np.arange(1,NC+1));np.testing.assert_array_equal(a['anchor_atom_id'],np.arange(NC+1,N+1));np.testing.assert_array_equal(a['reference_coordinates_m'],X[:NC]);np.testing.assert_array_equal(a['anchor_reference_coordinates_m'],X[NC:]);np.testing.assert_array_equal(a['support_node'],np.linalg.norm(a['total_force_N'],axis=1)>0);np.testing.assert_array_equal(x0,meta['reference_point_m']);assert a['support_node'].sum()==580
fullforce=np.zeros((N,3));fullforce[:NC]=a['total_force_N'];np.testing.assert_allclose(fullforce[:NC],a['pressure_force_N']+a['shear_force_N']+a['viscous_normal_force_N'],rtol=0,atol=1e-23)
components={};length=float(m['length_scale_m']);Q=s['quadrature_coordinates_m'];Aq=s['quadrature_area_m2'];normal=s['n_clot']
fields={'pressure':-s['quadrature_pressure_Pa'][:,:,None]*normal[:,None,:],'shear':np.repeat(s['tangential_traction_Pa'][:,None,:],3,axis=1),'viscous_normal':np.repeat(s['viscous_normal_traction_Pa'][:,None,:],3,axis=1),'total':s['quadrature_total_traction_Pa']}
for name,traction in fields.items():
 fq=traction*Aq[:,:,None];face=fq.sum(1);Fs=face.sum(0);Ms=np.cross(Q-x0,fq).sum((0,1));fn=a[name+'_force_N'];Fn=fn.sum(0);Mn=np.cross(X[:NC]-x0,fn).sum(0);scale=float(np.linalg.norm(face,axis=1).sum());fe=float(np.linalg.norm(Fn-Fs));me=float(np.linalg.norm(Mn-Ms));assert fe/scale<1e-8 and me/(scale*length)<1e-8
 components[name]={'surface_force_N':Fs.tolist(),'surface_moment_Nm':Ms.tolist(),'mapped_force_N':Fn.tolist(),'mapped_reference_moment_Nm':Mn.tolist(),'force_error_N':fe,'moment_error_Nm':me,'face_absolute_force_sum_N':scale,'length_scale_m':length,'force_scaled_error':fe/scale,'moment_scaled_error':me/(scale*length)}
r0=np.linalg.norm(X[b[:,1]]-X[b[:,0]],axis=1);h=np.where(abs(r0-DELTA)<=DX/2,(DELTA+DX/2-r0)/DX,1.);refm=np.bincount(b.ravel(),weights=np.repeat(V*r0*h,2),minlength=N)
# Read original data identities, not a newly generated geometry.
body=GDIR/'cases/C1_GAUSSIAN_A200/body.data';lines=body.read_text().splitlines();ia=next(i for i,l in enumerate(lines) if l.startswith('Atoms'));ib=next(i for i,l in enumerate(lines) if l.strip()=='Bonds');ar=np.array([[float(x) for x in l.split()] for l in lines[ia+1:ib] if l.strip()]);bb=np.array([[int(x) for x in l.split()] for l in lines[ib+1:] if l.strip()]);ar=ar[np.argsort(ar[:,0])];np.testing.assert_array_equal(ar[:,0],np.arange(1,N+1));np.testing.assert_allclose(ar[:,3:6],X,rtol=0,atol=1e-17);np.testing.assert_array_equal(bb[:,1],types);np.testing.assert_array_equal(bb[:,2:4]-1,b);assert np.all(ar[:NC,2]==1) and np.all(ar[NC:,2]==2)
# Detect nearest nonbonded pairs in a bounded radius; no contact assumptions from images.
p=cKDTree(X).query_pairs(2*DX,output_type='ndarray');pk=p[:,0]*N+p[:,1];bk=b[:,0]*N+b[:,1];p=p[~np.isin(pk,bk)];rr=np.linalg.norm(X[p[:,1]]-X[p[:,0]],axis=1)
from pd_clot.lammps_fragmentation.native_y_wall_query import NativeYWallQuery
wall=NativeYWallQuery();w=wall.query(X[:NC]);assert w['inside_lumen'].all()
np.savez_compressed(O/'REFERENCE_AND_LOAD.npz',X=X,bonds=b,bond_types=types,reference_m=refm,reference_lengths_m=r0,h=h,material_ids=a['material_atom_id'],anchor_ids=a['anchor_atom_id'],force_N=fullforce,support_node=a['support_node'],reference_point_m=x0,wall_clearance_m=w['clearance_m'],wall_normal=w['inward_normal'],near_nonbonded_pairs=p,near_nonbonded_reference_distances_m=rr)
check={'status':'PASS_SOURCE_SHA_REGISTRATION_AND_CONSERVATION','files':sourcecheck,'metadata_sha256':sha(MAP/'CLOT_NODAL_TRACTION_METADATA.json'),'source_CFD':meta['source_CFD'],'mapping_run':meta['run_id'],'physical_load_cases':1,'mapping_recomputed':False,'CFD_rerun':False,'clot_count':NC,'anchor_count':N-NC,'loaded_clot_count':580,'unloaded_clot_count':1347,'bonds':len(b),'bond_type_counts':{str(k):int(np.sum(types==k)) for k in [1,2,3]},'unit':meta['force_units'],'coordinates':meta['coordinates'],'coordinate_transform':meta['coordinate_transform'],'pressure_convention':meta['pressure'],'reference_point_m':x0.tolist(),'surface_faces':len(s['triangle_area_m2']),'surface_area_m2':float(s['triangle_area_m2'].sum()),'components':components,'nodal_force_L2_N':float(np.linalg.norm(fullforce)),'nodal_force_max_N':float(np.linalg.norm(fullforce,axis=1).max()),'unloaded_geometric_wall_min_clearance_m':float(w['clearance_m'].min()),'unloaded_nonbonded_min_distance_m_within_2dx':float(rr.min()),'wall_source_paths':str(wall.sources),'no_native_calls_yet':True}
dump(O/'CFD_LOAD_SOURCE_AUDIT.json',check)
files=[B/'liblammps.so',PROD,Path('/workspace/y_clot/env.sh'),PRIOR/'BUGFIX_SOURCE_DIFF.patch',R/'analysis/material_macro_audit_v1/energy.cpp',R/'analysis/material_macro_audit_v1/energy.so',R/'pd_clot/lammps_fragmentation/macro_lab_v1.py',V3/'observer.so',R/'analysis/final_v3_v1/observer.cpp']+[S/'src'/n for n in ['BPM/bond_bpm_peri.cpp','BPM/bond_bpm_peri.h','BPM/bond_bpm.cpp','BPM/pair_bpm_peri.cpp','fix_external.cpp','fix_setforce.cpp','verlet.cpp','balance.cpp']]+[R/'pd_clot/lammps_fragmentation'/n for n in ['lammps_writer.py','native_y_wall_contact.py','native_y_wall_query.py','early_flow_coupling.py']]
identity={str(p):sha(p) for p in files};assert identity[str(B/'liblammps.so')]==EXPECTED;assert identity[str(PROD)]=='d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a'
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',{'status':'DISK_IDENTITY_PASS_RUNTIME_PENDING','files_sha256':identity,'expected_library_sha256':EXPECTED,'production_sha256':sha(PROD),'original_C1_maps':[json.loads((V3/'C1_final'/f'LOADED_LIBRARY_rank{i}.json').read_text()) for i in [0,1]],'native_CXX_modified':False})
prov=O/'provenance';prov.mkdir();inp=O/'source_inputs';inp.mkdir()
for p in files:
 if '.so' not in p.name:shutil.copy2(p,prov/p.name)
for p in [R/'AGENTS.md',R/'PROJECT_ACTIVE_CONTEXT.md',Path('/workspace/y_clot/README_MIGRATION_ZH.txt'),Path('/etc/vast-agents-guide.md')]:shutil.copy2(p,prov/p.name)
for p in [MAP/meta['file'],MAP/'CFD_TO_LAMMPS_MAPPING_MATRIX.npz',MAP/'CLOT_NODAL_TRACTION_METADATA.json',MAP/'CONSERVATIVE_MAPPING_AUDIT.json',MAP/'INDEPENDENT_VERIFICATION.json',CFD/'refined/CLOT_SURFACE_TRACTIONS.npz',GDIR/'V3_MATERIAL_POINTS_AND_BONDS.npz',body]:shutil.copy2(p,inp/p.name)
shutil.copy2(V3/'observer.so',O/'observer.so');shutil.copy2(V3/'C1_final/in.lammps',prov/'original_C1.in.lammps')
prev=R/'verification/lps_fracture_timelayer_and_timestep_sensitivity_v1/20261011T221502Z';protected=json.loads((prev/'HISTORICAL_PROTECTION_BEFORE.json').read_text())
for base in [prev,MAP,CFD]:
 for p in base.rglob('*'):
  if p.is_file():protected[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
for p in files+list(map(Path,sources)):
 protected[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
# Verify inherited historical entries before adopting this round's protection baseline.
assert all(Path(p).is_file() and sha(p)==rec['sha256'] for p,rec in protected.items())
dump(O/'HISTORICAL_PROTECTION_BEFORE.json',protected)
print(json.dumps({'source':check['status'],'components':components,'Fmax':check['nodal_force_max_N'],'F_L2':check['nodal_force_L2_N'],'wall_clearance':check['unloaded_geometric_wall_min_clearance_m'],'nonbonded_distance':check['unloaded_nonbonded_min_distance_m_within_2dx'],'bond_types':check['bond_type_counts'],'protected_files':len(protected),'protected_bytes':sum(q['bytes'] for q in protected.values()),'wall_s':time.perf_counter()-start}))
