"""Read-only native handoff / timewarp audit. No solver launch entrypoints."""
import csv,json,time
from pathlib import Path
import numpy as np
from .config import ROOT,dump
from .y_clot_morphology_v3 import OUT as PRIOR,RESULTS as PRIOR_RESULTS,load_geometry
from .y_bifurcation_sources import YSources,digest
from .background_flow_sampler import BackgroundFlowSampler
from .failure_map import csv_out
from .fragment_analysis import read_dump

TASK='y_clot_continuous_handoff_and_honest_render_v1'
OUT=ROOT/'verification'/TASK
RESULTS=ROOT/'results'/TASK
CASE=PRIOR/'cases/C1_GAUSSIAN_A200'

def quant(x):
    a=np.asarray(x,float);return dict(min=float(a.min()),median=float(np.median(a)),max=float(a.max())) if len(a) else None

def angle(a,b):
    d=np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(np.dot(a,b)/d,-1,1)))) if d>1e-30 else None

def load_transport_steps():
    data={}
    for r in csv.DictReader((PRIOR/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv').open()):
        cid=int(r['component_id']);data.setdefault(cid,[]).append([float(r[k]) for k in ['relative_time_s','dt_s','x_m','y_m','z_m','vx_m_s','vy_m_s','vz_m_s','flow_tetra_id','Re_proxy']])
    return {cid:np.asarray(rows) for cid,rows in data.items()}

def register():
    s=YSources();paths=[PRIOR/'FINAL_STATUS.json',PRIOR/'REVIEW_Y_CLOT_REALISTIC_MORPHOLOGY_WALL_CONTACT_AND_FLOW_TRANSPORT_V1.md',CASE/'trajectory.lammpstrj',CASE/'surviving_bonds.local',CASE/'broken_bonds.local',CASE/'accepted_trajectory.npz',CASE/'broken_bonds.npz',CASE/'lineage/membership_and_shapes.json',CASE/'lineage/lifetimes.csv',CASE/'in.lammps',PRIOR/'FRAGMENT_HANDOFF.json',PRIOR/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz',PRIOR/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv',PRIOR/'FRAGMENT_TRANSPORT_EVENTS.json',PRIOR/'TRANSPORT_CONFIG_FROZEN.json',PRIOR_RESULTS/'movie_frame_manifest/frames.json',PRIOR_RESULTS/'animations/y_clot_realistic_fracture_and_transport.mp4',PRIOR_RESULTS/'figures/solver_vs_medical_view_same_frames.png']
    paths+=list((OUT/'SOURCE_BEFORE').glob('*.py'))
    dump(OUT/'SOURCE_AND_PRIOR_RESULT_REGISTER.json',dict(prior_verification=str(PRIOR),prior_results=str(PRIOR_RESULTS),files=[dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,sha256=digest(p)) for p in paths],native_mesh_sha256=s.mesh_id,native_flow_sha256=s.flow_id,native_geometry_sha256=s.geometry_id,SI_units=True,original_asset_register=str(PRIOR/'SOURCE_ASSET_REGISTER.json'),prior_wall_scope=str(PRIOR/'INDEPENDENT_WALL_AND_TRANSPORT_AUDIT.json'),policy='Read-only reuse. Liu K/G/s00,1.0dx,three anchor layers,geometry,CFD fixed; no B1/B2,PD2,CFD,FSI or bubble solves'))

def audit():
    started=time.perf_counter();register();g=load_geometry();n=g['clot_count'];tr=np.load(CASE/'accepted_trajectory.npz');hand=json.loads((PRIOR/'FRAGMENT_HANDOFF.json').read_text());h=hand['components'];trans=np.load(PRIOR/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz');data=load_transport_steps();flow=BackgroundFlowSampler();summary=[];detail=[];hist=[];cons=json.loads((PRIOR/'FRAGMENT_TRANSPORT_EVENTS.json').read_text())
    assert len(h)==208 and sum(r['node_count'] for r in h)==225
    # Read the native file, not just a derived JSON, at the last five actual saved times.
    native={}
    for step,cols,a in read_dump(CASE/'trajectory.lammpstrj'):
        if step>=tr['step'][-5]:native[step]=(a[:,[cols.index(k) for k in ['x','y','z']]],a[:,[cols.index(k) for k in ['vx','vy','vz']]])
    for j,r in enumerate(h):
        ids=np.array(r['member_atom_ids'])-1;cid=r['component_id'];x=tr['position_m'][-1,ids].mean(0);v=tr['velocity_m_s'][-1,ids].mean(0);m=r['mass_kg'];assert int(trans['component_id'][j])==cid
        for k in range(len(tr['step'])-5,len(tr['step'])):
            step=int(tr['step'][k]);X,V=native[step];assert np.array_equal(X,tr['position_m'][k]) and np.array_equal(V,tr['velocity_m_s'][k]);labels=tr['component_id'][k];same=bool(len(np.unique(labels[ids]))==1 and (labels==labels[ids[0]]).sum()==len(ids));independent=same and not tr['attached'][k,ids].any()
            hist.append(dict(component_id=cid,node_count=len(ids),member_atom_ids=json.dumps(r['member_atom_ids']),source_step=step,physical_time_s=step*.5e-6,exact_component_at_this_frame=same,detached_at_this_frame=bool(independent),COM_x_m=float(X[ids,0].mean()),COM_y_m=float(X[ids,1].mean()),COM_z_m=float(X[ids,2].mean()),COM_vx_m_s=float(V[ids,0].mean()),COM_vy_m_s=float(V[ids,1].mean()),COM_vz_m_s=float(V[ids,2].mean()),history_semantics='Exact component if flags true; otherwise only final member cohort, not an independent prior fragment'))
        posjump=float(np.linalg.norm(trans['position_m'][0,j]-x));vjump=float(np.linalg.norm(trans['velocity_m_s'][0,j]-v));handjump=float(np.linalg.norm(np.array(r['COM_velocity_mps'])-v));previous=x.copy();previous_v=v.copy();rowdata=data[cid]
        # Follow at least first five accepted steps, and all saved early steps through20ms.
        count=max(5,int(np.searchsorted(rowdata[:,0],.02,side='right')));u0=flow.sample([x])['velocity_m_s'][0]
        for k,z in enumerate(rowdata[:count]):
            t,dt=z[:2];newx=z[2:5];newv=z[5:8];sample=flow.sample([previous,newx]);assert sample['valid'].all();u=sample['velocity_m_s'][0];uend=sample['velocity_m_s'][1];re=1000*r['equivalent_diameter_m']*np.linalg.norm(previous_v-u)/.001
            tau=.0001 if len(ids)<=2 else 1000*r['equivalent_diameter_m']**2/(18*.001*(1+.15*re**.687 if re<1000 else .44*re/24))
            expected=u+(previous_v-u)*np.exp(-dt/tau);force=m*(u-previous_v)/tau;avgforce=m*(expected-previous_v)/dt
            contact=[e for e in cons['events'] if e['component_id']==cid and e['event']=='WALL_CONSTRAINT' and abs(e['time_s']-t)<1e-11]
            entry=dict(component_id=cid,node_count=len(ids),mass_kg=m,member_atom_ids=json.dumps(r['member_atom_ids']),category=r['category'],accepted_step_index=k+1,transport_elapsed_s=t,physical_time_s=.00233+t,dt_s=dt,relaxation_s=tau,Re_proxy=re,CFD_tetra_start=int(sample['tetra_id'][0]),wall_constraint=bool(contact),speed_m_s=float(np.linalg.norm(newv)),turn_from_handoff_deg=angle(v,newv),pure_drag_velocity_residual_m_s=float(np.linalg.norm(expected-newv)),drag_scope='Reconstructed applied model at actual pre-step state; not a newly solved trajectory; numerical proxy for1-2nodes')
            for prefix,a in [('COM_m',newx),('velocity_m_s',newv),('u_CFD_start_m_s',u),('u_CFD_end_m_s',uend),('drag_start_N',force),('drag_step_average_N',avgforce)]:
                for axis,value in zip('xyz',a):entry[prefix+'_'+axis]=float(value)
            detail.append(entry);previous=newx;previous_v=newv
        # Exact first saved accepted state at or after five100us packet relaxation times.
        k500=min(int(np.searchsorted(rowdata[:,0],.0005)),len(rowdata)-1);z500=rowdata[k500];z20=rowdata[min(int(np.searchsorted(rowdata[:,0],.02)),len(rowdata)-1)]
        native_prev=tr['position_m'][-2,ids].mean(0);laststep_disp=float(np.linalg.norm(x-native_prev));speed0=float(np.linalg.norm(v));u0speed=float(np.linalg.norm(u0))
        item=dict(component_id=cid,node_count=len(ids),category=r['category'],mass_kg=m,member_atom_ids=json.dumps(r['member_atom_ids']),handoff_step=4660,handoff_physical_time_s=.00233,position_jump_m=posjump,velocity_jump_m_s=vjump,handoff_JSON_velocity_error_m_s=handjump,initial_speed_m_s=speed0,initial_CFD_speed_m_s=u0speed,initial_angle_to_CFD_deg=angle(v,u0),first_step_s=float(rowdata[0,0]),first_step_speed_m_s=float(np.linalg.norm(rowdata[0,5:8])),first_step_direction_change_deg=angle(v,rowdata[0,5:8]),sample_after_500us_s=float(z500[0]),sample_after_500us_speed_m_s=float(np.linalg.norm(z500[5:8])),speed_ratio_after_500us=float(np.linalg.norm(z500[5:8])/max(speed0,1e-30)),sample_after_500us_angle_deg=angle(v,z500[5:8]),sample_near20ms_s=float(z20[0]),sample_near20ms_speed_m_s=float(np.linalg.norm(z20[5:8])),last_native_saved_step_displacement_m=laststep_disp,initial_CFD_displacement_over100us_m=u0speed*1e-4,initial_CFD_100us_pixels_old_full=u0speed*1e-4/.026*512,initial_CFD_100us_pixels_old_local=u0speed*1e-4/.0021*512,early_steps_audited=count)
        for prefix,a in [('native_handoff_COM_m',x),('native_handoff_velocity_m_s',v),('transport_t0_COM_m',trans['position_m'][0,j]),('transport_t0_velocity_m_s',trans['velocity_m_s'][0,j])]:
            for axis,value in zip('xyz',a):item[prefix+'_'+axis]=float(value)
        item['handoff_JSON_position_error_m']=float(np.linalg.norm(np.asarray(r['COM_xyz_m'])-x))
        summary.append(item)
    csv_out(OUT/'handoff_velocity_audit.csv',summary);csv_out(OUT/'last_five_LAMMPS_COM_frames.csv',hist);csv_out(OUT/'early_transport_native_steps_with_CFD_and_drag.csv',detail)
    groups={}
    for name,select in [('resolved_ge3',[r for r in summary if r['node_count']>=3]),('unresolved_1_2',[r for r in summary if r['node_count']<=2])]:
        groups[name]=dict(component_count=len(select),max_position_jump_m=max(r['position_jump_m'] for r in select),max_velocity_jump_m_s=max(r['velocity_jump_m_s'] for r in select),handoff_speed_m_s=quant([r['initial_speed_m_s'] for r in select]),CFD_speed_m_s=quant([r['initial_CFD_speed_m_s'] for r in select]),ratio_after_500us=quant([r['speed_ratio_after_500us'] for r in select]),direction_change_after_500us_deg=quant([r['sample_after_500us_angle_deg'] for r in select]))
    core=tr['attached'][-1];corev=tr['velocity_m_s'][-1,:n][core];report=dict(status='PASS',max_handoff_position_jump_m=max(r['position_jump_m'] for r in summary),max_handoff_velocity_jump_m_s=max(r['velocity_jump_m_s'] for r in summary),groups=groups,handoff_component_count=len(summary),total_handoff_mass_kg=sum(r['mass_kg'] for r in summary),native_last_five_frames_reread=True,actual_transport_steps_audited=len(detail),max_unconstrained_velocity_reproduction_error_m_s=max(r['pure_drag_velocity_residual_m_s'] for r in detail if not r['wall_constraint']),first_steps_with_wall_constraint=sum(r['wall_constraint'] for r in detail),all_detached_speeds_m_s=quant([r['initial_speed_m_s'] for r in summary]),old_full_view_pixel_m=.026/512,old_local_view_pixel_m=.0021/512,components_CFD_100us_motion_subpixel_old_local=sum(r['initial_CFD_100us_pixels_old_local']<1 for r in summary),components_CFD_speed_below_1mm_s=sum(r['initial_CFD_speed_m_s']<.001 for r in summary),attached_core=dict(nodes=int(core.sum()),mass_fraction=float(core.mean()),native_final_COM_speed_m_s=float(np.linalg.norm(corev.mean(0))),native_final_node_speed_m_s=quant(np.linalg.norm(corev,axis=1)),transport_role='Held at native handoff position; attached core is intentionally not advected. Its internal motion is outside the rigid detached-component transport model.'),cause_classification=dict(TRUE_VELOCITY_RESET_BUG='NOT_OBSERVED',PHYSICAL_MODEL_FORCE_DISCONTINUITY='ESTABLISHED: background fluid drag absent in fracture and enabled at handoff',PRESENTATION_TIMEWARP='See VIDEO_TIMEWARP_AUDIT.json',NEAR_WALL_SLOW_FLOW='LOCALLY_PRESENT: see CFD speeds and subpixel bounds',UNRESOLVED_DEBRIS_PROXY_RELAXATION='ESTABLISHED:100us relaxation, inherited initial velocity'),BACKGROUND_FLOW_NOT_APPLIED_DURING_FRACTURE=True,conditional_gate_D=dict(new_transport_solve='NOT_TRIGGERED: no reset, units or initialization error found',comparison_LAMMPS_plus_transport='NOT_RUN: continuously recomputed native component detachment/drag coupling has not been implemented or verified',future_contract='An explicit proposed coupling contract, not an implemented force'),wall_seconds=time.perf_counter()-started)
    report['cause_classification']['ATTACHED_BODY_SNAPSHOT_TRUNCATION']='ESTABLISHED: attached1702node cohort has native final COM speed0.0547157427m/s; Phase2 does not integrate this body and displays its saved position. This is a declared model truncation, not a numerical reset of transported components or evidence of physical stopping.'
    report['sampling_after_500us_elapsed_s']=quant([r['sample_after_500us_s'] for r in summary])
    report['speed_ratio_sampling_caveat']='First real accepted state at or after500us, not a common-time500us reconstruction; inspect each elapsed time in CSV.'
    assert report['max_handoff_position_jump_m']<1e-15 and report['max_handoff_velocity_jump_m_s']<1e-12
    assert report['max_unconstrained_velocity_reproduction_error_m_s']<1e-10
    dump(OUT/'HANDOFF_CONTINUITY_AUDIT.json',report)
    movie=json.loads((PRIOR_RESULTS/'movie_frame_manifest/frames.json').read_text());times=np.array([a['LAMMPS_time_s'] if a['chapter']=='LAMMPS_FRACTURE' else .00233+a['transport_relative_time_s'] for a in movie]);delta=np.diff(times);duplicates=np.flatnonzero(abs(delta)<1e-14)
    timewarp=dict(prior_fps=20,prior_total_frames=len(movie),old_handoff_frame_pairs=[[int(i),int(i+1)] for i in duplicates],old_extra_duplicate_duration_s=len(duplicates)/20,old_same_pose_display_span_s=(len(duplicates)+1)/20 if len(duplicates) else 0,old_multi_frame_hold_beyond_duplicate_endpoint=False,phase1_physical_seconds_per_playback_second=quant(delta[:117]*20),phase2_physical_seconds_per_playback_second=quant(delta[118:]*20),typical_physical_rate_ratio=2000.,phase1_typical_slowdown=2500.,phase2_typical_slowdown=1.25,clock_reset='Old panel switches from2.330ms to0.00s after handoff (relative clock, not data reset)',old_transport_first_nonzero_render_time_s=.04,missing_fast_relaxation_on_old_20ms_output_grid=True,planned_correction='Single endpoint; continuous absolute physical time; real accepted transport steps bracketing display-only same-member interpolation; report per-frame playback factor.')
    dump(OUT/'VIDEO_TIMEWARP_AUDIT.json',timewarp);dump(OUT/'handoff_timeline_audit.json',dict(prior=timewarp,audit=report,old_frame_absolute_times_s=times))
    print('HANDOFF_AUDIT',report['max_handoff_position_jump_m'],report['max_handoff_velocity_jump_m_s'],groups,flush=True)

if __name__=='__main__':audit()
