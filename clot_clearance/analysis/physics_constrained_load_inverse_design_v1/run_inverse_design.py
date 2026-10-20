"""One CLI for freeze/preflight and resumable real optimizer stages."""
from common import *
import argparse
def prepare():
 O.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True)
 from traction_library import build,BOUNDS
 from traction_to_v3 import build as mapping
 build();mapping()
 frozen=dict(status='FROZEN_BEFORE_FIRST_NATIVE_TRIAL',mechanism='STREAMING_SHEAR_UNSTEADY_PLANE_STOKES',parameters={k:dict(bounds=v,SI_unit='m/s' if k=='drive_speed_m_s' else 's',source='Prescribed moving upper plane, independent of unknown acoustic/bubble conversion',extrapolation=False) for k,v in BOUNDS.items()},scope='Low-Re physical plane-flow screen; curved V3 patch template has unverified geometry error',small_geometry_source=str(SPEC),material='Unchanged original validated variational LPS and precrack; all material nodes free to move, original bottom anchor layers fixed',load='Plane shear along +x on existing flat top; trapezoidal area quadrature; no displacement grip; no numerical damping or background drag',objective='history max(mean d on remaining z=0 bulk ligament - mean d on internal bonds of upper-half candidate + 0.01*mean clip(kappa/s0,0,1) on remaining ligament). Original precuts excluded. If no damage, this is only a tensile-strain precursor, not crack evidence.',objective_masks='candidate original material z>0; boundary original bulk bonds crossing z=0 with initial_cut false; internal original bulk bonds inside candidate',algorithm='3 seeded random initial native evaluations; 1 actual exploratory constrained BoTorch proposal/evaluation; another next-candidate suggestion only. Not statistically established Bayesian optimization.',numerical_checks='Existing VV forecast max damage increment .1, dt .25us to .03125us. Original-ID native connectivity, irreversibility, force port and energy ledger.',high_fidelity_admission='At most one new V3 only if meaningful localized new cohesive damage and valid physical traction; otherwise NOT_RUN_EVIDENCE_INSUFFICIENT',library_sha256=sha(LIB),reused_source_sha256={str(p.relative_to(R)):sha(p) for p in [R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',SPEC/'reference.npz',SPEC/'reference.native',SPEC/'body.data']},implementation_sha256={p.name:sha(p) for p in A.glob('*.py')})
 p=O/'DESIGN_FREEZE.json'
 if not p.exists():dump(p,frozen)
 print('PREPARED',sha(p),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--limit',type=int,default=4);p.add_argument('--next',action='store_true');args=p.parse_args()
 if args.prepare:prepare()
 else:
  assert (O/'DESIGN_FREEZE.json').exists()
  from inverse_optimizer import run
  run(args.limit,args.next)
