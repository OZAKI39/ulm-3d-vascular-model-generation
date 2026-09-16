"""Package frozen first RMBW outputs; compare production C++ log interpolation.
No modification, symmetrization, fitting or smoothing of reference data.
"""
from pathlib import Path
import sys,json,hashlib,ctypes,csv
import numpy as np,h5py
R=Path(__file__).resolve().parents[1];P=json.loads((R/'provenance/TASK_PATHS.json').read_text());W=Path(P['local_work'])
sys.path.insert(0,str(R/'src'))
from reference_wall_hydrodynamics_v0 import Lookup
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
plan=json.loads((R/'contracts/TABLE_GENERATION_PLAN.json').read_text());t=np.load(R/'raw/RMBW_table_FIRST.npz');eps=t['epsilon'];T=t['R_total_scaled'];B=t['M_scaled'];E=T-np.eye(6)
path=R/'tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5';assert not path.exists(),'No table overwrite'
classes=np.where(eps<=.2,1,np.where(eps<=5,2,3)).astype(np.int32)
lim=['finite-gap RR absolute accuracy unverified','TR/RT continuous-gap accuracy unverified','high-order far-field RR/TR not certified','RR_parallel literature constant 0.3817 versus 0.3709 unresolved']
with h5py.File(path,'w') as f:
 for name,d in [('epsilon',eps),('M_RMBW_SI',t['M_SI']),('R_total_RMBW_SI',t['R_total_SI']),('M_scaled',B),('R_total_scaled',T),('R_wall_excess_scaled',E),('M_scaled_eigenvalues',np.linalg.eigvalsh((B+B.transpose(0,2,1))/2)),('R_total_scaled_eigenvalues',np.linalg.eigvalsh((T+T.transpose(0,2,1))/2)),('reciprocity_error',np.max(abs(T-T.transpose(0,2,1)),axis=(1,2))/np.max(abs(T),axis=(1,2))),('condition_number',np.linalg.cond(T)),('validity_class',classes)]:f.create_dataset(name,data=d,compression='gzip')
 f.create_dataset('reference_mode_flags',data=np.array([json.dumps({'class':int(c),'normal':'Brenner supported','TT_parallel':'qualified','RR_TR':'REFERENCE_LIMITATIONS','far_fallback':bool(e>9.018296)}) for c,e in zip(classes,eps)],dtype=h5py.string_dtype()))
 f.attrs['format']='RMBW_WALL_RESISTANCE_TABLE_V0';f.attrs['reference_radius_m']=float(t['radius_m']);f.attrs['reference_mu_Pa_s']=float(t['mu_Pa_s']);f.attrs['convention_sha256']=plan['convention_sha256'];f.attrs['limitations']=json.dumps(lim);f.attrs['runtime_uses']='R_wall_excess_scaled only';f.attrs['generator_commit']=plan['reference_commit']
lib=ctypes.CDLL(str(W/'libwall_math.so'));lib.wall_lookup_open.argtypes=[ctypes.c_char_p];lib.wall_lookup_open.restype=ctypes.c_void_p;handle=lib.wall_lookup_open(str(path).encode());assert handle
ptr=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS');lib.wall_lookup_evaluate.argtypes=[ctypes.c_void_p,ctypes.c_int,ptr,ptr]
hold=np.array(plan['holdout_epsilon']);interp=np.empty((len(hold),6,6));assert lib.wall_lookup_evaluate(handle,len(hold),hold,interp)==0;interp+=np.eye(6)
independent=Lookup(path).excess(hold)+np.eye(6);assert np.allclose(interp,independent,rtol=1e-13,atol=1e-14)
rng=np.random.default_rng(plan['seed']+1);vectors=rng.normal(size=(128,6));vectors/=np.linalg.norm(vectors,axis=1)[:,None]
out=[];maxerr=0.;worst=None
for size in plan['radii_m']:
 raw=np.load(R/'raw'/('RMBW_'+size+'_FIRST.npz'));direct=raw['R_total_scaled'];mob=raw['M_scaled'];inv=np.linalg.solve(interp,np.broadcast_to(np.eye(6),interp.shape))
 # Work-conjugate normalized vectors cover both imposed velocities and loads.
 ref_actions=direct@vectors.T;act=interp@vectors.T;error=np.max(np.linalg.norm(act-ref_actions,axis=1)/np.maximum(np.linalg.norm(ref_actions,axis=1),1e-300),axis=1)
 ref_vel=mob@vectors.T;vel=inv@vectors.T;verror=np.max(np.linalg.norm(vel-ref_vel,axis=1)/np.maximum(np.linalg.norm(ref_vel,axis=1),1e-300),axis=1)
 entry=np.max(abs(interp-direct)/np.maximum(np.maximum(abs(interp),abs(direct)),1e-30),axis=(1,2))
 eig=np.linalg.eigvalsh((interp+interp.transpose(0,2,1))/2);eigref=np.linalg.eigvalsh((direct+direct.transpose(0,2,1))/2);eigerr=np.max(abs(eig-eigref)/eigref,axis=1)
 rec=np.max(abs(interp-interp.transpose(0,2,1)),axis=(1,2))/np.max(abs(interp),axis=(1,2))
 for i,e in enumerate(hold):
  row={'size':size,'epsilon':e,'action_relative_error':error[i],'mobility_action_relative_error':verror[i],'max_entry_relative_error':entry[i],'eigenvalue_relative_error':eigerr[i],'reciprocity_error':rec[i],'min_scaled_eigenvalue':eig[i,0],'qualified':bool(e<=5)}
  for name,k,l in [('normal_TT',2,2),('parallel_TT',0,0),('RR_parallel',3,3),('RR_normal',5,5),('TR',0,4)]:row[name+'_relative_error']=abs(interp[i,k,l]-direct[i,k,l])/max(abs(interp[i,k,l]),abs(direct[i,k,l]),1e-30)
  out.append(row)
  if e<=5 and max(error[i],verror[i])>maxerr:maxerr=max(error[i],verror[i]);worst=row
np.savez_compressed(R/'validation/TABLE_HOLDOUT_CPP.npz',epsilon=hold,R_interpolated_scaled=interp,action_vectors=vectors)
with (R/'validation/TABLE_HOLDOUT_ERRORS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
qualified=[x for x in out if x['qualified']];summary={'status':'PASS' if maxerr<=1e-3 and min(x['min_scaled_eigenvalue'] for x in qualified)>0 and max(x['reciprocity_error'] for x in qualified)<=1e-10 else 'FAIL',
 'grid_points':len(eps),'holdout_count_per_radius':len(hold),'radii_tested':list(plan['radii_m']),'vectors_per_epsilon':len(vectors),'max_qualified_action_error':maxerr,'worst':worst,
 'max_reciprocity_error':max(x['reciprocity_error'] for x in qualified),'min_scaled_eigenvalue':min(x['min_scaled_eigenvalue'] for x in qualified),'min_scaled_mobility_eigenvalue':float(np.min(np.linalg.eigvalsh(B))),
 'table_target':'fresh frozen RMBW outputs, not independently exact Stokes truth','production_interpolation':'C++ shared-library path','python_independent_interpolation':'PASS'}
(R/'validation/TABLE_INTERPOLATION_AUDIT.json').write_text(json.dumps(summary,indent=2)+'\n')
contract={'status':'FROZEN','table_file':str(path.relative_to(R)),'table_sha256':sha(path),'RMBW_commit':plan['reference_commit'],'source_provenance':json.loads((R/'provenance/reference_audit/RIGIDMULTIBLOBSWALL_PROVENANCE.json').read_text()),
 'generator_sha256':plan['generator_sha256'],'packager_sha256':sha(Path(__file__)),'reference_audit_manifest_sha256':sha(R/'provenance/reference_audit/SHA256SUMS'),
 'grid_points':len(eps),'epsilon_range':[.001,20.],'grid_dataset':'epsilon','matrix_convention':json.loads((R/'provenance/reference_audit/WALL_REFERENCE_CONVENTION.json').read_text()),
 'units':'D=diag(sqrt(MT0),sqrt(MR0)); R_SI=D^-1 R_scaled D^-1; excess=total-I','interpolation':'piecewise linear in log(epsilon), full 6x6, no coefficient compression',
 'reference_classes':{'NEAR_FIELD_HIGH_CONFIDENCE':[.001,.2],'QUALIFIED_REFERENCE':[.001,5.],'LIMITED_REFERENCE':[5.,20.]},'limitations':lim,
 'epsilon_below_min':'clamp resistance to epsilon=0.001; independent hard nonoverlap','epsilon_above_max':'zero excess; DEVELOPMENT BULK HANDOFF','runtime_RMBW_dependency':False,'bulk_double_counting':False}
(R/'contracts/RMBW_WALL_TABLE_CONTRACT.json').write_text(json.dumps(contract,indent=2)+'\n')
handoff={'epsilon':20.,'closure':'DEVELOPMENT BULK HANDOFF','jump_definition':'R_total_scaled(20+)-R_total_scaled(20); zero wall excess above20',
 'normal_TT_jump':float(1-T[-1,2,2]),'parallel_TT_jump':float(1-T[-1,0,0]),'RR_parallel_jump':float(1-T[-1,3,3]),'RR_normal_jump':float(1-T[-1,5,5]),'TR_jump':float(-T[-1,0,4]),
 'max_matrix_action_jump_upper_bound':float(np.linalg.norm(T[-1]-np.eye(6),2)),'trajectory_effect':'MEASURE_IN_SYNTHETIC_HANDOFF_CASE','not_exact_zero_wall_physics':True}
(R/'validation/FAR_FIELD_HANDOFF_AUDIT.json').write_text(json.dumps(handoff,indent=2)+'\n')
print(json.dumps(summary,indent=2));assert summary['status']=='PASS','Preserve failure; do not change reference values'
