"""Three instantaneous contact checks; no state advancement or changed model."""
import csv,json
import numpy as np
from analyze_stops import (OUT,RUN,VIS,read,dump,csvout,FrozenFEMField,WallGeometry,Sphere,
    NearFieldRegularizationV1,handoff_constraints,assemble_v1,augment_planar_system,wall_gap,
    solve_resistance,contact_jacobian,lp_direction,sha)

def main():
    rows=list(csv.DictReader((OUT/'data/stopped_tracks.csv').open()))
    rows=[r for r in rows if r['downstream_direction_test']=='DOWNSTREAM_DIRECTION_EXISTS']
    a=np.load(RUN/'data/gpu_mesh_input.npz');field=FrozenFEMField(a['points'],a['tetra'],a['velocity'],a['pressure'])
    wall=WallGeometry(a['wall_triangles']);policy=NearFieldRegularizationV1();result=[]
    for r in rows:
        pid=int(r['particle_id']);folder=RUN/'tracks'/f'mb_{pid:06d}'
        last=read(folder/'support.json')['last_32_accepted_contact_solves'][-1]
        position=np.array(last['position']);shape=Sphere(position,float(r['diameter_um'])*.5e-6)
        # Use the original binary radius rather than rounded CSV algebra.
        shape=Sphere(position,read(folder/'metrics.json')['radius_m']);shapes={pid:shape};s=field.sample(position)
        free=np.r_[s.velocity_m_s,.5*s.vorticity_s_inv]
        system=assemble_v1(shapes,{pid:free},.00345312,wall,policy=policy)
        system=augment_planar_system(system,shapes,.00345312,wall,lambda x:field.sample(x).velocity_gradient_s_inv,wall_gap)
        # The saved nearest feature is interior; no open-rim fallback applies.
        assert all(v['applicability']=='NEAREST_SINGLE_PLANAR_WALL' for v in system.planar_diagnostics)
        assert all(v['applicability']=='NEAREST_SINGLE_PLANAR_WALL' for v in last['planar'])
        contacts,_=handoff_constraints(shapes,wall,policy);sol=solve_resistance(system,particles=shapes,constraints=contacts)
        j,ordered=contact_jacobian(system,shapes,contacts);n=j[:,:3]
        R=system.matrix.toarray();b=system.rhs
        force=b[:3]-R[:3,3:]@np.linalg.solve(R[3:,3:],b[3:])
        test=lp_direction(n,force);assert test['status']=='NO_DOWNSTREAM_DIRECTION_CERTIFIED'
        fluid=lp_direction(n,s.velocity_m_s);assert fluid['status']=='DOWNSTREAM_DIRECTION_EXISTS'
        unitF=force/np.linalg.norm(force);unitU=s.velocity_m_s/np.linalg.norm(s.velocity_m_s)
        balance=force+n.T@np.array(sol.record['multipliers'])
        result.append(dict(particle_id=pid,diameter_um=float(r['diameter_um']),
            fluid_forward_max_projection=fluid['maximum_forward_projection'],
            reduced_force_forward_max_projection=test['maximum_forward_projection'],
            reduced_force_no_forward_certificate_residual=test['dual_certificate_residual'],
            reduced_force_vs_fluid_angle_deg=float(np.degrees(np.arccos(np.clip(unitF@unitU,-1,1)))),
            force_balance_residual_N=float(np.linalg.norm(balance)),
            normalized_force_dot_fluid_forward_witness=float(unitF@fluid['direction']),
            replay_speed_m_s=float(np.linalg.norm(sol.velocity[:3]))))
    csvout('three_forward_exceptions.csv',result)
    dump('data/three_forward_exceptions.json',dict(all_pass=True,count=len(rows),results=result,
        definition='F_eff=b_v-R_vw solve(R_ww,b_w), instantaneous force after eliminating free rotation',
        scope='The current resistance-weighted stationary optimum need not maximize dot(FEM_velocity,v). No trajectory or pressure experiment.',
        original_data_changed=False,new_trajectory_integrations=0))
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
