"""Read-only scientific audit adapters; never imported by production physics."""
from pathlib import Path
from copy import deepcopy
import csv, json, hashlib, time, socket, os
import numpy as np
from scipy.optimize import linprog
from .particle81_simulation import environment, dump, DT, HORIZON
from .particle8_replay import canonical_hash

ROLES = ['OUTLET_01','OUTLET_02','OUTLET_03','NO_EXIT']
FLOW_SHA = '129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d'

def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def csvwrite(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
        w.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)

def verify_events(candidates,events):
    raw={r['common_event']['particle_id']:r for r in candidates}
    assert {e['particle_id'] for e in events} == {i for i,r in raw.items() if r['admission']['accepted']}
    for e in events:
        r=raw[e['particle_id']]
        assert canonical_hash(r['common_event'])==e['common_event_sha256']
        assert canonical_hash(r['admission'])==e['admission_sha256']
        assert e['anchor_m']==e['birth_center_m']==r['common_event']['anchor_m']==r['admission']['birth_center_m']
        assert e['birth_time_s']==r['scheduled_time_s']
        assert e['radius_m']==r['admission']['radius_m']
    return dict(candidate_count=len(raw),accepted_count=len(events),same_initial_position=True,events_unchanged=True)

def p65_job(args):
    event,folder,identity=args
    from .particle82a_integration import integrate_admitted
    start=time.time();before=canonical_hash(event)
    m=integrate_admitted(deepcopy(event),output=folder)
    assert canonical_hash(event)==before==m['birth_metadata_sha256']
    meta=Path(folder)/'trajectories'/f"mb_{event['particle_id']:06d}.json"
    m['audit_identity']=identity;m['audit_role']='SAME_500_P65_CURRENT_HANDOFF_ONLY_P9_CORRECTION_OFF'
    m['receipt']=dict(host=socket.gethostname(),pid=os.getpid(),start_unix_s=start,end_unix_s=time.time())
    dump(meta,m)
    return dict(particle_id=event['particle_id'],outlet=m['exit_outlet'] or 'NO_EXIT',outcome=m['end_reason'],failure_detail=m['failure_detail'],radius_m=m['radius_m'],minimum_gap_m=m.get('minimum_original_wall_gap_m'),travel_time_s=m.get('last_elapsed_time_s'),trajectory_sha256=m['samples_sha256'],wall_seconds=m['wall_seconds'],birth_metadata_sha256=before)

def point_job(args):
    rows,folder,config=args
    from .particle82_tracers import trace_positions
    env=environment();result=[]
    for raw in rows:
        e=raw['common_event'];pid=e['particle_id'];p=e['anchor_m']
        t=trace_positions([p],**config)[0];path=t.pop('path')
        assert np.array_equal(path[0,1:],p)
        f=Path(folder)/f'point_{pid:06d}.npz';f.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(f,path=path)
        face=e['anchor_triangle'];ad=raw['admission']
        result.append(dict(candidate_id=pid,particle_id=pid,birth_time_s=raw['scheduled_time_s'],anchor_m=p,anchor_triangle=face,
            triangle_flux_weight_m3_s=float(env.sampler.weights[face]),triangle_sampling_probability=float(env.sampler.weights[face]/env.sampler.Q_m3_s),
            sampling_provenance='EXISTING_EXACT_POSITIVE_P1_FLUX_DIRICHLET_MIXTURE; RNG STATES IN ORIGINAL COMMON_EVENT',
            common_event_sha256=canonical_hash(e),admission_sha256=canonical_hash(ad),first_diameter_um=e['first_diameter_um'],first_radius_m=e['first_radius_m'],
            accepted=ad['accepted'],admission_status=ad['status'],rejection_causes=sorted(set(a['cause'] for a in ad['attempts'] if a['cause']!='ACCEPTED')) if not ad['accepted'] else [],
            accepted_radius_m=ad['radius_m'],point_outlet=t['outlet'] or 'NO_EXIT',point_end_reason=t['end_reason'],
            travel_time_s=float(path[-1,0]),sample_count=len(path),initial_position_exact=True,terminal_classifier_verified=bool(t['outlet']),path_sha256=sha(f)))
    return result

def downstream_feasibility(normals,velocity):
    """Max uhat.v subject to N v >= 0 and box |v_j| <= 1; save witness + LP dual."""
    N=np.asarray(normals,float);u=np.asarray(velocity,float);speed=np.linalg.norm(u)
    if not speed: return dict(status='ZERO_FEM_UNRESOLVED',feasible_downstream_direction=None)
    u=u/speed;lp=linprog(-u,A_ub=-N,b_ub=np.zeros(len(N)),bounds=[(-1,1)]*3,method='highs',options={'dual_feasibility_tolerance':1e-10,'primal_feasibility_tolerance':1e-10})
    if not lp.success:return dict(status=lp.message,feasible_downstream_direction=None)
    v=lp.x;tol=256*np.finfo(float).eps*max(1,np.linalg.cond(N))*max(1,len(N))
    primal=float(np.min(N@v));value=float(u@v)
    valid=primal>=-tol
    # A no-downstream certificate: -uhat is a nonnegative combination of normals.
    from scipy.optimize import nnls
    lam,res=nnls(N.T,-u);cert=bool(res<=tol)
    feasible=bool(value>tol and valid)
    resolved=feasible or cert
    return dict(status='RESOLVED' if resolved else 'OPTIMIZATION_TOLERANCE_UNRESOLVED',feasible_downstream_direction=feasible if resolved else None,
        direction=v.tolist(),unit_direction=(v/np.linalg.norm(v)).tolist() if np.linalg.norm(v) else [0.,0.,0.],
        max_downstream_projection=value,minimum_normal_projection=primal,verification_tolerance=tol,
        no_downstream_certificate=cert,certificate_multipliers=lam.tolist(),certificate_residual=float(res),
        lp_inequality_dual=lp.ineqlin.marginals.tolist(),objective='max dot(unit FEM velocity,v), Nv>=0, -1<=v_j<=1')
