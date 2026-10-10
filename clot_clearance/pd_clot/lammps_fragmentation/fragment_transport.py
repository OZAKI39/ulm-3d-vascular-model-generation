"""One-way native P1 CFD, exponential drag and native wall constrained translations."""
import csv,json,time,resource
import numpy as np
from .config import dump
from .y_clot_morphology_v3 import OUT
from .native_y_wall_query import NativeYWallQuery
from .failure_map import csv_out

CONFIG=dict(model='ONE_WAY_UNOBSTRUCTED_BACKGROUND_CFD_DEMONSTRATION',rho_f_kg_m3=1000.,mu_f_Pa_s=.001,end_relative_time_s=5.,maximum_step_s=.005,maximum_displacement_m=20e-6,output_interval_s=.02,wall_center_proxy_m=5e-6,unresolved_packet_relaxation_s=.0001,maximum_constraint_iterations=12,drag='Schiller-Naumann: factor=1+.15Re^.687 for Re<1000; .44Re/24 otherwise',drag_source='https://cpp.openfoam.org/v12/SchillerNaumann_8C_source.html',constraint='Frictionless inelastic unilateral constraint within every integrator step; rigid member offsets retained; no rotation, no renderer correction',near_wall_uncertain_clearance_m=100e-6)

def freeze():
    assert not (OUT/'TRANSPORT_CONFIG_FROZEN.json').exists();dump(OUT/'TRANSPORT_CONFIG_FROZEN.json',{**CONFIG,'freeze_epoch_s':time.time()})
    (OUT/'BACKGROUND_CFD_LIMITATION.md').write_text('''# ONE_WAY_UNOBSTRUCTED_BACKGROUND_CFD_DEMONSTRATION
Only real C1 detached graph components are handed off at4660 steps (2.33ms), with exact members, reference mass, and actual COM velocity. Attached bodies remain fixed at their C1 final positions. Native result_121.vtu is unmodified, unobstructed water flow; it does not supply the artificial fracture traction.
Every motion step queries strict P1 tetra velocity at current COM. Frozen-field exponential drag update uses adaptive dt<=5ms and displacement<=20um; output times are reached by actual integration, not interpolated LAMMPS frames. Resolved graph components (>=3 nodes) translate rigidly, with Stokes inertia time rho*d_eq^2/(18mu) and Schiller-Naumann finite-Re factor. Source: https://cpp.openfoam.org/v12/SchillerNaumann_8C_source.html . This empirical sphere correction is an engineering choice for irregular fragments.
1-2node debris are numerical mass packets with declared100us relaxation, inherited velocity and retained mass; no clinical diameter or biological particle inertia is inferred from50um grid spacing. Their d_eq remains only a reference-volume graph statistic.
Native WALL constraints are solved inside the ODE step on every rigid member center;5um engineering clearance proxy, reduced only at initialization to existing positive clearance if necessary (no initial position change). Candidate constrained translation and normal velocity removal are recorded. Sweep intersections are checked. Actual port crossings bypass sidewall constraints; a component exits when its COM crosses an actual native port.10mm daughter classification differs from12mm outlet.
No fluid feedback, obstruction update, interfragment collision, collision with the retained clot, fragment rotation/deformation, lift, added mass, turbulence, adhesion, acoustic or bubble force. Near-wall CFD gradients and drag are not independently validated. Paths near/through the retained clot are flagged uncertain. No transport coordinate change occurs in the renderer.
''')

def run():
    started=time.perf_counter();assert not (OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz').exists()
    c=json.loads((OUT/'TRANSPORT_CONFIG_FROZEN.json').read_text());hand=json.loads((OUT/'FRAGMENT_HANDOFF.json').read_text());rows=hand['components'];assert rows
    q=NativeYWallQuery();flow=q.flow;src=flow.sources;tr=np.load(OUT/'cases/C1_GAUSSIAN_A200/accepted_trajectory.npz');X0=tr['position_m'][-1];attached=X0[:int(tr['clot_count'])][tr['attached'][-1]]
    from scipy.spatial import cKDTree
    retained_tree=cKDTree(attached)
    # Native interpolation roundtrip: exact tetra centroids have all weights1/4.
    check_ids=np.arange(0,len(src.tet),max(1,len(src.tet)//120));centers=np.asarray(src.flow.points[src.tet[check_ids]],float).mean(1);sample=flow.sample(centers);expected=src.u[src.tet[check_ids]].mean(1);err=float(abs(sample['velocity_m_s']-expected).max());assert sample['valid'].all() and err<1e-12
    outputs=np.round(np.arange(0,c['end_relative_time_s']+c['output_interval_s']/2,c['output_interval_s']),8);C=len(rows);P=np.empty((len(outputs),C,3));V=np.empty_like(P);active=np.zeros((len(outputs),C),bool);states=[];events=[];step_rows=[];summ=[]
    for j,r in enumerate(rows):
        pos=np.array(r['COM_xyz_m']);v=np.array(r['COM_velocity_mps']);ids=np.array(r['member_atom_ids'])-1;offsets=X0[ids]-pos;initial=q.query(pos+offsets);minimum=float(initial['clearance_m'].min());assert minimum>=0,'Handoff contains actual wall penetration'
        radius=min(c['wall_center_proxy_m'],max(0,minimum-1e-9));t=0.;steps=0;contacts=0;maxcorrection=0.;maxre=0.;maxsweep=0;near=minimum<c['near_wall_uncertain_clearance_m'];overlap=False;status='TRACKING';classified=[];minaccepted=minimum
        member_ports=['']*len(ids)
        def legal_member_exit(ii,old,new):
            remembered=member_ports[ii]
            if remembered:
                port=q.ports[remembered]
                if np.dot(new-port['center'],port['outward'])>0:return True
            crossed=q.port_crossing(old,new)
            if crossed:member_ports[ii]=crossed;return True
            return False
        P[0,j]=pos;V[0,j]=v;active[0,j]=True
        for frame,target in enumerate(outputs[1:],1):
            while t<target-1e-12 and status=='TRACKING':
                fld=flow.sample([pos]);assert fld['valid'][0],f'Invalid CFD COM at component{r["component_id"]}'
                u=fld['velocity_m_s'][0];speed=max(np.linalg.norm(v),np.linalg.norm(u),1e-9);dt=min(c['maximum_step_s'],c['maximum_displacement_m']/speed,target-t)
                d=r['equivalent_diameter_m'];re=c['rho_f_kg_m3']*d*np.linalg.norm(v-u)/c['mu_f_Pa_s'];maxre=max(maxre,re)
                if r['node_count']>=3:
                    factor=1+.15*re**.687 if re<1000 else .44*re/24;tau=1000*d*d/(18*c['mu_f_Pa_s']*factor)
                else:tau=c['unresolved_packet_relaxation_s']
                fac=-np.expm1(-dt/tau);vn=u+(v-u)*(1-fac);dis=u*dt+(v-u)*tau*fac;trial=pos+dis;port=q.port_crossing(pos,trial)
                if port:
                    normal=q.ports[port]['outward'];center=q.ports[port]['center'];fraction=-np.dot(pos-center,normal)/np.dot(trial-pos,normal);trial=pos+(trial-pos)*fraction;dt*=fraction;vn=u+(v-u)*np.exp(-dt/tau);status='EXIT_'+port;events.append(dict(component_id=r['component_id'],time_s=t+dt,event=status,xyz_m=trial.tolist(),mass_kg=r['mass_kg']))
                else:
                    # Projection is the explicit unilateral constraint substep, never postprocessing.
                    correction=np.zeros(3);normals=[]
                    for attempt in range(c['maximum_constraint_iterations']):
                        test=q.query(trial+offsets);clear=test['clearance_m'].copy()
                        for ii in np.flatnonzero(clear<0):
                            if legal_member_exit(ii,pos+offsets[ii],trial+offsets[ii]):clear[ii]=np.inf
                        deficit=radius-clear;worst=int(np.argmax(deficit))
                        if deficit[worst]<=1e-10:break
                        nn=test['inward_normal'][worst];shift=(deficit[worst]+1e-9)*nn;trial+=shift;correction+=shift;normals.append(nn)
                    else:
                        status='STALLED_CONSTRAINT';trial=pos.copy();vn=np.zeros(3);events.append(dict(component_id=r['component_id'],time_s=t,event=status,xyz_m=pos.tolist(),mass_kg=r['mass_kg']))
                    for nn in normals:vn-=min(np.dot(vn,nn),0)*nn
                    if normals:
                        contacts+=1;corr=float(np.linalg.norm(correction));maxcorrection=max(maxcorrection,corr)
                        events.append(dict(component_id=r['component_id'],time_s=t+dt,event='WALL_CONSTRAINT',xyz_m=trial.tolist(),mass_kg=r['mass_kg'],constraint_correction_m=corr))
                    final=q.query(trial+offsets);validclear=final['clearance_m'].copy()
                    for ii in np.flatnonzero(validclear<0):
                        if legal_member_exit(ii,pos+offsets[ii],trial+offsets[ii]):validclear[ii]=np.inf
                    assert validclear.min()>=-1e-10,dict(component=r['component_id'],t=t,pos=pos.tolist(),trial=trial.tolist(),clearance=validclear.tolist(),member_ports=member_ports)
                    length=np.linalg.norm(trial-pos);old=q.query(pos+offsets);hits=0
                    for ii in np.flatnonzero(old['clearance_m']<=length+1e-12):
                        if q.swept_wall_crossing(pos+offsets[ii],trial+offsets[ii])[0]:hits+=1
                    if hits:
                        # Reject a swept crossing within the integrator. Keep the last valid state.
                        events.append(dict(component_id=r['component_id'],time_s=t,event='STALLED_SWEPT_WALL_GUARD',xyz_m=pos.tolist(),mass_kg=r['mass_kg']));status='STALLED_SWEPT_WALL_GUARD';trial=pos.copy();vn=np.zeros(3)
                    maxsweep+=hits;minaccepted=min(minaccepted,float(validclear.min()));near|=validclear.min()<c['near_wall_uncertain_clearance_m']
                for section in q.classification_crossing(pos,trial):
                    if section not in classified:classified.append(section);events.append(dict(component_id=r['component_id'],time_s=t+dt,event=section,xyz_m=trial.tolist(),mass_kg=r['mass_kg']))
                overlap|=bool(retained_tree.query(trial+offsets)[0].min()<50e-6)
                pos=trial;v=vn;t+=dt;steps+=1
                step_rows.append(dict(component_id=r['component_id'],relative_time_s=t,absolute_time_s=t+.00233,dt_s=dt,x_m=pos[0],y_m=pos[1],z_m=pos[2],vx_m_s=v[0],vy_m_s=v[1],vz_m_s=v[2],flow_tetra_id=int(fld['tetra_id'][0]),Re_proxy=re,status=status))
            P[frame,j]=pos;V[frame,j]=v;active[frame,j]=status=='TRACKING'
        if status=='TRACKING':status='NOT_REACHED_CLASSIFICATION_PLANE' if not classified else 'CLASSIFIED_NOT_EXITED'
        summ.append(dict(component_id=r['component_id'],category=r['category'],node_count=r['node_count'],mass_kg=r['mass_kg'],final_status=status,classification_sections=classified,end_time_s=t,steps=steps,wall_constraint_steps=contacts,maximum_constraint_correction_m=maxcorrection,minimum_accepted_center_clearance_m=minaccepted,maximum_Re_proxy=maxre,Re_proxy_physical_meaning='hydrodynamic approximation' if r['node_count']>=3 else 'grid-volume diagnostic only; packet law uses declared tau',near_wall_uncertain=bool(near),retained_clot_proximity_uncertain=bool(overlap),swept_guard_rejections=maxsweep,proxy_radius_m=radius,displacement_m=float(np.linalg.norm(pos-r['COM_xyz_m']))))
        print('TRANSPORT',j+1,'/',C,r['component_id'],r['node_count'],status,t,steps,contacts,flush=True)
    np.savez_compressed(OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz',relative_time_s=outputs,absolute_time_s=outputs+.00233,component_id=[r['component_id'] for r in rows],position_m=P,velocity_m_s=V,active=active)
    csv_out(OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv',step_rows);csv_out(OUT/'FRAGMENT_TRANSPORT_SUMMARY.csv',summ)
    fields=['component_id','time_s','event','xyz_m','mass_kg','constraint_correction_m']
    with (OUT/'FRAGMENT_TRANSPORT_EVENTS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(events)
    dump(OUT/'FRAGMENT_TRANSPORT_EVENTS.json',dict(events=events,components=summ,transported_mass_kg=sum(r['mass_kg'] for r in rows),all_mass_retained_or_exit_accounted=True,no_accepted_center_wall_penetration=all(r['minimum_accepted_center_clearance_m']>=-1e-10 for r in summ),model=c['model'],configuration=c,wall_seconds=time.perf_counter()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,wall_query_seconds=q.query_seconds,flow_queries=flow.queries))
    dump(OUT/'FLOW_SAMPLER_AND_DOMAIN_VALIDITY.json',dict(status='PASS',source_sha256=src.flow_id,mesh_sha256=src.mesh_id,field='velocity',units='m/s',method='P1 native tetra barycentric interpolation, no extrapolation',centroid_checks=len(check_ids),max_interpolation_error_m_s=err,maximum_source_velocity_m_s=float(np.linalg.norm(src.u,axis=1).max()),out_of_domain_evaluations=flow.rejected,note='Includes rejected candidate geometry probes; no accepted velocity extrapolation',components=summ))

if __name__=='__main__':
    import sys
    {'freeze':freeze,'run':run}[sys.argv[1]]()
