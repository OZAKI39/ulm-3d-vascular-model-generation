from common import *
import sys,os,time,subprocess,shutil
start=time.perf_counter()
assert not (O/'PRECHECK_FREEZE.json').exists()
O.mkdir(parents=True,exist_ok=True); F.mkdir(parents=True,exist_ok=True)
old=json.loads((SURFACE/'HISTORICAL_PROTECTION_BEFORE.json').read_text())
for p,v in old.items(): assert Path(p).is_file() and sha(p)==v['sha256'],p
protect=old.copy()
for base in [SURFACE,R/'results'/C['prior_surface'].split('/')[1]/C['prior_surface'].split('/')[2],R/'analysis/fracture_surface_v1']:
 for p in base.rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:protect[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
dump(O/'HISTORICAL_PROTECTION_BEFORE.json',protect)
commands={}
for label,cmd in [('git_status',['git','status','--short']),('git_head',['git','rev-parse','HEAD']),('cpu',['lscpu']),('memory',['free','-b']),('disk',['df','-B1','/workspace']),('GPU',['nvidia-smi','--query-gpu=name,memory.total,memory.used,utilization.gpu','--format=csv']),('processes',['ps','-eo','pid,ppid,etimes,pcpu,pmem,comm'])]:
 q=subprocess.run(cmd,capture_output=True,text=True);commands[label]={'command':cmd,'returncode':q.returncode,'stdout':q.stdout,'stderr':q.stderr}
paths=['/sys/fs/cgroup/cpu/cpu.cfs_quota_us','/sys/fs/cgroup/cpu/cpu.cfs_period_us','/sys/fs/cgroup/memory/memory.limit_in_bytes']
dump(O/'PREFLIGHT.json',{'commands':commands,'cgroup':{p:Path(p).read_text() for p in paths},'python':sys.executable,'threads':1,'GPU_used':False})
ident=json.loads((SURFACE/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text())
for p,v in ident['files'].items():assert sha(p)==v['sha256'],p
ident.update(actual_LAMMPS_instance_loaded=False,scope='File identity verified; runtime identity will be recorded only if a native check or formal case launches.',previous_task_reference=str(SURFACE))
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',ident)
rt=np.load(W/'reference_topology.npz'); assert sha(W/'reference_topology.npz')=='ef145b6d6db6c17584dce6cbea1fa1b9d95847d71a69265a73df6334d1695c03'
shutil.copy2(W/'reference_topology.npz',O/'reference_topology.npz')
shutil.copy2(W/'SAVED_CANDIDATE_SPATIAL_REGIONS.json',O/'HISTORICAL_GRIP_COMPETITION.json')
assert np.array_equal(rt['bond_reference_id'],np.arange(1,len(rt['bonds'])+1))
X=rt['X'];b=rt['bonds'];alive=rt['alive'][0].astype(bool);cross=rt['cross_z0'];m=rt['m'];d=X[b[:,1]]-X[b[:,0]];d[:,1]-=PERIOD*np.rint(d[:,1]/PERIOD);r=np.linalg.norm(d,axis=1);h=np.where(r>DELTA-DX/2,(DELTA+DX/2-r)/DX,1.);mm=np.bincount(b.ravel(),weights=np.repeat(V*r*h,2),minlength=len(X))
assert len(X)==16000 and len(b)==909392 and np.count_nonzero(~alive)==16416 and cross.sum()==55408
assert np.array_equal(~alive,cross&(rt['reference_intersection_x_m']<-.001-1e-12))
assert np.max(abs(mm-m)/m)<1e-12 and len(np.unique(b[:,0]*len(X)+b[:,1]))==len(b)
np.savez_compressed(O/'precut_and_loading_identity.npz',precut_reference_bond_ids=rt['bond_reference_id'][~alive],top_atom_ids=np.flatnonzero(rt['high'])+1,bottom_atom_ids=np.flatnonzero(rt['low'])+1,reference_m=m)
dump(O/'SPECIMEN_AND_PRECRACK_AUDIT.json',{'status':'PASS','source_sha256':sha(W/'reference_topology.npz'),'nodes':len(X),'unique_full_reference_bonds':len(b),'periodic_Y_seam_bonds':int(np.count_nonzero(rt['periodic_image_y'])),'m_definition':'Full pristine family; preserved after all initial cuts and any native breaks','m_recompute_max_relative':float(np.max(abs(mm-m)/m)),'precut_bonds':int((~alive).sum()),'initial_alive_bonds':int(alive.sum()),'all_cross_z0_bonds':int(cross.sum()),'initial_surviving_cross_z0_bonds':int((cross&alive).sum()),'precut_surviving_cross_z0_bonds':0,'dimensions_X_Y_Z_m':[.005,.0004,.001],'initial_crack_length_m':.0015,'tip_X_m':-.001,'remaining_ligament_m':.0035,'initial_precrack_area_m2':.0015*PERIOD,'boundary':'finite X,Z; periodic Y','point_indexing':'atom ID = ((ix*8+iy)*20+iz)+1','dx_m':DX,'horizon_m':DELTA,'K_Pa':K,'G_Pa':G,'rho_kg_m3':1000,'s00':S00,'alpha':0})
freeze={'selection':'One symmetric, uniform dead-traction loading of outer top and bottom material layers; no axial displacement constraint, lateral DOFs free. No further fixture variant allowed this turn.','physical_rationale':'Apply distributed normal force on the physical outer faces, avoiding the rigid common axial displacement and its known end-corner constraint. This is a candidate improvement, not an assumed remedy.','precheck_type':'One non-damaging independent elastic equilibrium solve, no topology changes or native fracture','traction_Pa':20.,'force_per_loaded_node_N':20*DX**2,'top_and_bottom_nodes_each':int(rt['high'].sum()),'face_area_m2':.005*PERIOD,'initial_coordinates':'Exact pristine reference after explicit A0 pre-cuts, zero initial displacement and velocity','rigid_motion_handling':'Balanced forces, exact symmetric starting state; no local artificial physical pins. Center translation removed only as a gauge in reporting, not internal deformation.','initial_guess':'Undeformed, no parameter fitting','energy_implementation':'Existing validated analysis/dual_track_v1/energy_y.so, force-compatible expanded LPS with fixed full-family m','geometry_sha256':sha(O/'reference_topology.npz'),'material':{'K_Pa':K,'G_Pa':G,'rho_kg_m3':1000,'s00':S00,'alpha':0},'optimizer':{'method':'L-BFGS-B','coordinate_scale_m':1e-6,'energy_scale_J':1e-12,'maxiter':5000,'maxfun':15000,'ftol':2.220446049250313e-15,'gtol':1e-10,'maxcor':30,'maxls':40},'budget':{'wall_s':300,'peak_RSS_MiB':4096,'output_MiB':200},'gates':{'force_max_abs_N':1e-13,'force_L2_relative_to_external':1e-5,'noncontact_equilibrium_requires_final_contact_U_J_below':1e-24,'max_bond_stretch':.03,'symmetry_load_relative':1e-12,'front_leads_every_competing_region':True},'regions':'Exactly prior fixed horizon-wide reference regions; notch front prioritized, then notch mouth, grip, other finite-X surface, interior','stop':'Do not launch a formal fracture if this one precheck is unaccepted, front is not the leading tensile-stretch region, contact is unresolved, or physical budget cannot be supported. Do not rerun another loading variant.','formal_candidate_if_admitted':{'dt_s':5e-7,'loading_duration_s':.024,'ramp':'half-cosine symmetric traction, amplitude must be frozen from the accepted elastic precheck before any formal run','damping':'none (0); measure actual kinetic content','maximum_wall_s':1800,'maximum_RSS_GiB':8,'maximum_output_GiB':2,'MPI_ranks':2},'created_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
dump(O/'PRECHECK_FREEZE.json',freeze)
dump(O/'PREPARATION_RESOURCES.json',{'wall_s':time.perf_counter()-start,'protected_files':len(protect),'protected_bytes':sum(v['bytes'] for v in protect.values())})
print('PRECHECK_FROZEN',O,flush=True)
