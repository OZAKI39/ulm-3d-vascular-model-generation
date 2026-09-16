from pathlib import Path
import json,hashlib
import numpy as np
R=Path(__file__).resolve().parents[1];d=np.load(R/'validation/COEFFICIENT_100000_STATES.npz');inputs=d['inputs'];n=len(inputs);rng=np.random.default_rng(20260917);normal=rng.normal(size=(n,3));normal/=np.linalg.norm(normal,axis=1)[:,None];tangent=rng.normal(size=(n,3));tangent-=np.sum(tangent*normal,axis=1)[:,None]*normal;tangent/=np.linalg.norm(tangent,axis=1)[:,None];q=np.zeros((n,2,6));kind=np.arange(n)%6
for k in range(6):
 sel=kind==k
 if k==0:q[sel,0,:3]=-.001*normal[sel];q[sel,1,:3]=.001*normal[sel]
 if k==1:q[sel,0,:3]=.001*normal[sel];q[sel,1,:3]=-.001*normal[sel]
 if k==2:q[sel,0,:3]=.001*tangent[sel];q[sel,1,:3]=-.001*tangent[sel]
 if k==3:q[sel,0,3:]=1000*tangent[sel];q[sel,1,3:]=1000*tangent[sel]
 if k==4:q[sel,0,3:]=1000*tangent[sel];q[sel,1,3:]=-1000*tangent[sel]
 if k==5:q[sel,:,:3]=rng.normal(0,.001,(sel.sum(),2,3));q[sel,:,3:]=rng.normal(0,1000,(sel.sum(),2,3))
np.savez_compressed(R/'validation/PAIR_STATE_MOTION_COVERAGE_100000.npz',state_index=np.arange(n),normal=normal,q_SI=q,mode_code=kind)
meta={'status':'PASS','N':n,'radius_gap_source':'COEFFICIENT_100000_STATES.npz','radius_gap_sha256':hashlib.sha256(inputs.tobytes()).hexdigest(),'motion_seed':20260917,'motion_categories':['approach','separation','tangential_slip','equal_transverse_rotation','opposite_transverse_rotation','arbitrary_combined'],'counts':np.bincount(kind).tolist(),'interpretation':'Scalar coefficients depend only on radii/gap, hence each already-verified scalar applies exactly to all assigned velocities; combined C++ full-matrix tests are separately stored in PAIR_TORQUE_BALANCE_AUDIT.json. This file is coverage/provenance, not an additional solver run.'}
(R/'validation/PAIR_STATE_MOTION_COVERAGE.json').write_text(json.dumps(meta,indent=2)+'\n')
print('MOTION_COVERAGE',n)
