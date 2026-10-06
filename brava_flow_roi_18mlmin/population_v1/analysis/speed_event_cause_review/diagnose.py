"""Original CUDA RHS on stored windows; original CPU solver for exact accounting."""
from pathlib import Path
import json
import signal
import sys
import time
import dataclasses
import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
POP = HERE.parents[1]
sys.path.insert(0, str(POP / 'source/particle_3d/src'))
import torch
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.gpu_field import GPUField
from particle_3d.population_v1.gpu_wall import GPUWall
from particle_3d.population_v1.gpu_resistance import assemble_isolated, solve_isolated
from particle_3d.resistance_solver import solve_resistance
from particle_3d.resistance_assembly import ResistanceSystem
from particle_3d.population_v1.components import assemble_components, solve_components
from particle_3d.nearfield_handoff import handoff_constraints
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.resistance_solver import contact_jacobian
from particle_3d.particle_shapes import Sphere
from diagnostic_core import held_alignment, increment_loads, segment_cells, strict_probe, direction_angle

torch.set_num_threads(1)


def save(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def npv(x):
    return x.detach().cpu().numpy()


class OriginalRHS:
    def __init__(self):
        identity = json.loads((POP / 'runs/realization_0000_contactfix/IDENTITY.json').read_text())
        self.env = load_environment(identity['flow_provenance'])
        # In-memory accelerators; no writes in the production data folder.
        self.field = GPUField(self.env.field)
        self.wall = GPUWall(self.env.wall)

    def batch(self, positions, radius, hints=None):
        x = torch.as_tensor(np.asarray(positions), dtype=torch.float64, device='cuda')
        r = torch.full((len(x),), radius, dtype=x.dtype, device=x.device)
        hints = None if hints is None else torch.as_tensor(hints, dtype=torch.int64, device=x.device)
        sampled = self.field.sample(x, hints)
        cell, data = sampled[:2]
        if bool((cell < 0).any()):
            raise ValueError('INVALID_ORIGINAL_FIELD_SAMPLE')
        wall = self.wall.query(x, r)
        if bool((wall[1][:, 0] < max(2e-9, .001*radius)-self.wall.roundoff).any()):
            raise ValueError('ORIGINAL_HANDOFF_ADMISSIBILITY_FAILED')
        free = torch.cat([data[:, :3], .5*data[:, 13:16]], 1)
        R, b, diag, assembly = assemble_isolated(r, self.env.mu, wall[1], data[:, 4:13].reshape(-1, 3, 3), free)
        U, hydro, audit = solve_isolated(R, b, diag, free, wall[3], wall[4])
        fallback = audit['fallback'] | (assembly['affine_solve_info'] != 0)
        for index in torch.nonzero(fallback).flatten().cpu().tolist():
            # Exactly the production small-component authority; complete singleton,
            # original GPU free packet plus authoritative CPU gradient/handoff.
            shapes = {0: Sphere(npv(x[index]), radius)}
            policy = NearFieldRegularizationV1()
            system = assemble_components(shapes, {0: npv(free[index])}, self.env.mu,
                self.env.wall, policy=policy,
                gradient_at=lambda p:self.env.field.sample(p).velocity_gradient_s_inv)
            cs,_ = handoff_constraints(shapes, self.env.wall, policy)
            sol = solve_components(system, particles=shapes, constraints=cs)
            local = system.systems[0]
            J, ordered = contact_jacobian(local, shapes, cs)
            row = sol.record['components'][0]
            for target, value in [(R,local.matrix.toarray()),(b,local.rhs),(U,sol.velocity),(hydro,sol.unconstrained)]:
                target[index] = torch.as_tensor(value,device=x.device,dtype=x.dtype)
            wall[4][index]=len(cs);wall[3][index]=0;wall[2][index]=-1
            audit['multipliers'][index]=0
            if cs:
                wall[3][index,:len(cs)]=torch.as_tensor(J[:,:3],device=x.device,dtype=x.dtype)
                wall[2][index,:len(cs)]=torch.tensor([c.canonical_id[3] for c in ordered],device=x.device)
                audit['multipliers'][index,:len(cs)]=torch.tensor(row['multipliers'],device=x.device,dtype=x.dtype)
        q = torch.zeros_like(b)
        q[:, :3] = (wall[3]*audit['multipliers'][:, :, None]).sum(1)
        normal = wall[1][:, 7:10]
        normalized = normal/torch.linalg.vector_norm(normal, dim=1)[:, None]
        Rnormal = torch.zeros_like(R)
        Rnormal[:, :3, :3] = assembly['normal_coefficient'][:, None, None]*(normalized[:, :, None]*normalized[:, None, :])
        result = {k: npv(v) for k, v in dict(x=x, tetra=cell, flow=data, triangle=wall[0], wall=wall[1],
            contact_triangles=wall[2], contact_normals=wall[3], contact_count=wall[4],
            R=R, b=b, diag=diag, free=free, U=U, hydro=hydro, q=q,
            residual=(R@U[:, :, None]).squeeze(2)-b-q,
            b_bulk=diag*free, b_planar=b-diag*free, R_normal=Rnormal,
            R_planar=R-torch.diag_embed(diag)-Rnormal,
            multipliers=audit['multipliers'], backward_residual=audit['backward_residual'],
            cpu_fallback=fallback,
            wall_weight=assembly['wall_weight'], normal_coefficient=assembly['normal_coefficient']).items()}
        result['feature'] = np.sum(result['wall'][:, 10:13] > 128*np.finfo(float).eps, axis=1)
        G = result['flow'][:, 4:13].reshape(-1, 3, 3)
        gn = np.einsum('nij,nj->ni', G, result['wall'][:, 7:10])
        result['shear'] = gn - result['wall'][:, 7:10]*np.einsum('ni,ni->n', gn, result['wall'][:, 7:10])[:, None]
        return result


def state_packet(z, k, t=None):
    count = int(z['contact_count'][k])
    active = np.flatnonzero(z['multipliers'][k, :count] > 0)
    return dict(time_s=None if t is None else float(t[k]), x_m=z['x'][k].tolist(),
        speed_mm_s=float(np.linalg.norm(z['U'][k, :3])*1e3),
        hydro_speed_mm_s=float(np.linalg.norm(z['hydro'][k, :3])*1e3),
        u_mm_s=(z['flow'][k, :3]*1e3).tolist(), u_speed_mm_s=float(np.linalg.norm(z['flow'][k, :3])*1e3),
        G_s_inv=z['flow'][k, 4:13].reshape(3, 3).tolist(), curl_s_inv=z['flow'][k, 13:16].tolist(),
        shear_s_inv=z['shear'][k].tolist(), tetra=int(z['tetra'][k]),
        gap_um=float(z['wall'][k, 0]*1e6), normal=z['wall'][k, 7:10].tolist(),
        nearest_point_m=z['wall'][k, 1:4].tolist(), triangle=int(z['triangle'][k]),
        feature={3:'FACE', 2:'EDGE', 1:'VERTEX', 0:'VERTEX'}[int(z['feature'][k])],
        barycentric=z['wall'][k, 10:13].tolist(), wall_weight=float(z['wall_weight'][k]),
        normal_coefficient_kg_s=float(z['normal_coefficient'][k]),
        candidate_contact_triangles=z['contact_triangles'][k, :count].tolist(),
        active_contact_triangles=z['contact_triangles'][k, active].tolist(),
        multipliers_kinematic_N=z['multipliers'][k, :count].tolist(),
        U_translation_m_s=z['U'][k, :3].tolist(), U_rotation_s_inv=z['U'][k, 3:].tolist())


def account(z, i, j, t=None):
    loads, target, r0, r1 = increment_loads(z['R'][i], z['R'][j], z['b'][i], z['b'][j],
        z['U'][i], z['U'][j], z['q'][i], z['q'][j])
    bulk_delta=z['b_bulk'][j]-z['b_bulk'][i]
    parts = dict(delta_bulk_translation_b=np.r_[bulk_delta[:3],np.zeros(3)],
        delta_bulk_rotation_b=np.r_[np.zeros(3),bulk_delta[3:]],
        delta_planar_shear_b=z['b_planar'][j]-z['b_planar'][i],
        delta_normal_resistance=-(z['R_normal'][j]-z['R_normal'][i])@z['U'][i],
        delta_planar_resistance=-(z['R_planar'][j]-z['R_planar'][i])@z['U'][i],
        delta_constraint_q=loads[2], delta_solver_residual=loads[3])
    system = ResistanceSystem((0,), sparse.csr_matrix(z['R'][j]), z['b'][j], z['diag'][j],
                              z['free'][j], [], [])
    response = {}
    for name, load in parts.items():
        solved = solve_resistance(dataclasses.replace(system, rhs=load))
        response[name] = dict(load_translation_N=load[:3].tolist(), load_torque_N_m=load[3:].tolist(),
            translation_mm_s=(solved.velocity[:3]*1e3).tolist(), rotation_s_inv=solved.velocity[3:].tolist())
    actual = z['U'][j]-z['U'][i]
    reconstructed = sum(np.r_[np.array(r['translation_mm_s'])/1e3, r['rotation_s_inv']] for r in response.values())
    error = reconstructed-actual
    return dict(from_state=state_packet(z, i, t), to_state=state_packet(z, j, t),
        components=response, actual_translation_delta_mm_s=(actual[:3]*1e3).tolist(),
        actual_rotation_delta_s_inv=actual[3:].tolist(),
        scalar_speed_delta_mm_s=float((np.linalg.norm(z['U'][j,:3])-np.linalg.norm(z['U'][i,:3]))*1e3),
        direction_change_deg=direction_angle(z['U'][i,:3], z['U'][j,:3]),
        reconstruction_error_translation_m_s=float(np.max(abs(error[:3]))),
        reconstruction_error_rotation_s_inv=float(np.max(abs(error[3:]))),
        equation_load_error_force_N=float(np.max(abs(loads.sum(0)[:3]-target[:3]))),
        equation_load_error_torque_N_m=float(np.max(abs(loads.sum(0)[3:]-target[3:]))),
        role='Exact algebra using R_at_later_state and original solve_resistance; not unique causal percentages')


def topology_key(z, k, wall):
    nodes = wall.global_node_ids[int(z['triangle'][k])]
    active = z['wall'][k, 10:13] > 128*np.finfo(float).eps
    return tuple(sorted(int(i) for i in nodes[active]))


def interfaces(rhs, z, t, radius):
    field = rhs.env.field
    crossings, partitions = [], []
    for k in range(len(t)-1):
        spans, hits = segment_cells(field, z['x'][k], z['x'][k+1])
        partitions.append(dict(segment=k, cells=[v[2] for v in spans], fractions=[[v[0], v[1]] for v in spans]))
        for hit in hits:
            hit.update(segment=k, time_s=float(t[k]+hit['fraction']*(t[k+1]-t[k])), kind='TETRA')
            G0, G1 = field.gradients_s_inv[[hit['left'], hit['right']]]
            hit['G_jump_s_inv'] = float(np.linalg.norm(G1-G0))
            hit['priority'] = hit['G_jump_s_inv']*(radius*max(z['wall_weight'][k:k+2]) +
                np.linalg.norm(z['x'][k+1]-z['x'][k]))
            crossings.append(hit)
    wall_switches = []
    for k in range(len(t)-1):
        key0, key1 = topology_key(z, k, rhs.env.wall), topology_key(z, k+1, rhs.env.wall)
        active0 = tuple(z['contact_triangles'][k, :int(z['contact_count'][k])][z['multipliers'][k,:int(z['contact_count'][k])]>0])
        active1 = tuple(z['contact_triangles'][k+1,:int(z['contact_count'][k+1])][z['multipliers'][k+1,:int(z['contact_count'][k+1])]>0])
        if key0 != key1 or active0 != active1:
            wall_switches.append(dict(kind='WALL_FEATURE_OR_ACTIVE_CONSTRAINT', segment=k,
                bracket_time_s=[float(t[k]), float(t[k+1])], topology_before=list(key0), topology_after=list(key1),
                active_before=[int(i) for i in active0], active_after=[int(i) for i in active1],
                priority=float(np.linalg.norm(z['U'][k+1,:3]-z['U'][k,:3])),
                near_wall=bool(max(z['wall_weight'][k:k+2])>0 or active0 or active1)))
    # All crossings retained. Include near-wall derivative interfaces explicitly:
    # bulk G jumps do not necessarily enter translation when wall weights are zero.
    peak = int(np.argmax(np.diff(np.linalg.norm(z['U'][:,:3],axis=1))))
    chosen = [v for v in crossings if v['segment'] == peak]
    near_crossings = [v for v in crossings if max(z['wall_weight'][v['segment']:v['segment']+2])>0]
    for hit in sorted(near_crossings, key=lambda h:-h['G_jump_s_inv'])[:3]:
        if hit not in chosen:
            chosen.append(hit)
    for hit in sorted(crossings, key=lambda h:-h['priority']):
        if hit not in chosen and len(chosen) < 8:
            chosen.append(hit)
    near_switches=[v for v in wall_switches if v['near_wall']]
    releases=[v for v in near_switches if len(v['topology_before'])==3 and len(v['topology_after'])<3]
    ordered_switches=releases[:1]+sorted(near_switches,key=lambda h:-h['priority'])
    for hit in ordered_switches:
        if hit in chosen:continue
        if len(chosen) < 12:
            chosen.append(hit)
    measurements = []
    for hit in chosen:
        k = hit['segment']; a,b=z['x'][k:k+2]; d=b-a; length=np.linalg.norm(d)
        if hit['kind'] == 'TETRA':
            f=hit['fraction']; width=min(f-hit['left_fraction'], hit['right_fraction']-f)
        else:
            lo,hi=0.,1.
            topology_changes=hit['topology_before']!=hit['topology_after']
            def switch_key(packet,index):
                if topology_changes:return topology_key(packet,index,rhs.env.wall)
                count=int(packet['contact_count'][index])
                return tuple(packet['contact_triangles'][index,:count][packet['multipliers'][index,:count]>0])
            startkey=switch_key(z,k)
            for _ in range(40):
                mid=(lo+hi)/2
                sample=rhs.batch([a+mid*d],radius)
                if switch_key(sample,0)==startkey:
                    lo=mid
                else:
                    hi=mid
            f=(lo+hi)/2; width=min(f,1-f)
        entry=dict(interface=hit, located_time_s=float(t[k]+f*(t[k+1]-t[k])), approaches=[])
        for fraction in [.01,.001,.0001]:
            delta=width*fraction
            if delta*length <= 16*field.geometry.candidate_padding_m:
                entry['approaches'].append(dict(fraction=fraction,status='NOT_EVALUABLE',reason='Probe spacing below production locator roundoff margin'))
                continue
            positions=np.array([a+(f-delta)*d,a+(f+delta)*d])
            expected=[hit['left'],hit['right']] if hit['kind']=='TETRA' else [None,None]
            try:
                checks=[strict_probe(field,x,c) for x,c in zip(positions,expected)]
                probe=rhs.batch(positions,radius)
            except ValueError as error:
                entry['approaches'].append(dict(fraction=fraction,status='NOT_EVALUABLE',reason=str(error)))
                continue
            rec=account(probe,0,1)
            rec.update(fraction=fraction,status='MEASURED',separation_m=float(np.linalg.norm(positions[1]-positions[0])),
                barycentric=[v[1].tolist() for v in checks],
                delta_u_m_s=float(np.linalg.norm(probe['flow'][1,:3]-probe['flow'][0,:3])),
                delta_G_s_inv=float(np.linalg.norm(probe['flow'][1,4:13]-probe['flow'][0,4:13])),
                delta_gap_m=float(probe['wall'][1,0]-probe['wall'][0,0]),
                delta_normal=float(np.linalg.norm(probe['wall'][1,7:10]-probe['wall'][0,7:10])),
                delta_RHS_m_s=float(np.linalg.norm(probe['U'][1,:3]-probe['U'][0,:3])))
            entry['approaches'].append(rec)
        measurements.append(entry)
    weights=[]
    for threshold in [.01,.05,.1,1.]:
        ratios=z['wall'][:,0]/radius
        for k in np.flatnonzero((ratios[:-1]-threshold)*(ratios[1:]-threshold)<0):
            weights.append(dict(gap_over_radius=threshold,bracket_time_s=[float(t[k]),float(t[k+1])],
                                role='Accepted-segment bracket, not interpolated dynamics'))
    return dict(all_segment_partitions=partitions,all_tetra_crossings=crossings,
                wall_switches=wall_switches,weight_transition_brackets=weights,measurements=measurements)


def main():
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('WINDOW_600S_LIMIT')))
    started=time.perf_counter(); rhs=OriginalRHS()
    print('ENVIRONMENT_AND_ORIGINAL_CUDA_READY',time.perf_counter()-started,flush=True)
    for event in json.loads((HERE/'FROZEN_PLAN.json').read_text())['events']:
        pid=event['particle_id']; clock=time.perf_counter();signal.alarm(600)
        try:
            ledger=json.loads((HERE/f'ledger_{pid}.json').read_text())
            assert ledger['isolation_confirmed'] and not ledger['cpu_component_rows']
            with np.load(HERE/f'original_{pid}.npz') as raw:
                t=raw['physical_time_s']; x=raw['x_m']; vel=raw['velocity_m_s']; hints=raw['tetra_id']
                dt,held,alignment_error=held_alignment(t,x,vel)
                z=rhs.batch(x[:-1],ledger['radius_m'],hints[:-1])
            t=t[:-1]
            reproduction=np.max(abs(z['U'][:,:3]-held))
            if reproduction>1e-9:
                raise ValueError(('RHS_REPRODUCTION_FAILED',float(reproduction)))
            np.savez_compressed(HERE/f'measured_{pid}.npz',**z,time_s=t,dt_s=dt,
                                held_m_s=held,displacement_error_m=alignment_error)
            speeds=np.linalg.norm(z['U'][:,:3],axis=1)
            peak=int(np.argmax(np.diff(speeds)))
            if pid==822:
                before,after=np.searchsorted(t,[.6505,.651])
            else:
                before,after=np.searchsorted(t,[event['center_s']-.8*.069082516488968,event['center_s']+.8*.069082516488968])
            accounting=[dict(label='event_before_after',**account(z,int(before),int(after),t)),
                        dict(label='largest_adjacent_scalar_increase_in_window',**account(z,peak,peak+1,t))]
            inter=interfaces(rhs,z,t,ledger['radius_m'])
            result=dict(particle_id=pid,window_s=event['window_s'],actual_range_s=[float(t[0]),float(t[-1])],
                samples=len(t),radius_m=ledger['radius_m'],maximum_alignment_error_m=float(np.max(abs(alignment_error))),
                maximum_rhs_vs_held_error_m_s=float(reproduction),maximum_backward_residual=float(z['backward_residual'].max()),
                actual_dt_s=dict(min=float(dt.min()),median=float(np.median(dt)),max=float(dt.max()),
                                count_below_nominal=int(np.sum(dt<.0005-1e-12))),
                birth_row_has_initial_free_velocity=pid==822,accounting=accounting,interfaces=inter,
                original_rejected_trials_in_window=len(ledger['rejected_rows']),
                original_rejects_responsible_target=sum(pid in r.get('responsible_particle_ids',[]) for r in ledger['rejected_rows']),
                original_handoff_events=ledger['handoff_rows'],original_cpu_fallback_count=0,
                minimum_swept_neighbor_margin_m=ledger['minimum_swept_margin_above_pair_shell_m'],
                restart_limitation=ledger['state_restart_limit'],elapsed_s=time.perf_counter()-clock)
            save(f'result_{pid}.json',result)
            print('DONE',pid,'RHS_error',reproduction,'crossings',len(inter['all_tetra_crossings']),
                  'interface_tests',len(inter['measurements']),'seconds',result['elapsed_s'],flush=True)
        except Exception as error:
            import traceback
            save(f'error_{pid}.json',dict(error=repr(error),traceback=traceback.format_exc(),elapsed_s=time.perf_counter()-clock))
            print('STOP_BRANCH',pid,repr(error),flush=True)
        finally:
            signal.alarm(0)


if __name__=='__main__':
    main()
