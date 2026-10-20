"""Read actual MEVO assets and classify a mechanical demand without inventing a transfer law."""
from common import *
import h5py,argparse
M=Path('/workspace/mevo_acoustic_microbubbles')
EXEC=M/'f10s_control_validation/.execution/14731fd56f0477c89cbbdea5ae1f13e8e38e3f3c'
DATA=EXEC/'results/straight_tube_exp1'
def inspect_assets():
 response=DATA/'f10n/inputs/response.h5';trajectory=DATA/'f10s/a_navigation_v1/runs/X_PAIR_3MS_DT1/trajectory.h5';cache=DATA/'f10n/cache/coarse.h5';kwave=M/'project/results/fast_pzt/f1/inputs/basis.h5'
 with h5py.File(response) as h:lo=h['lower'][:];hi=h['upper'][:];valid=h['valid'][:].astype(bool);harm=h['harmonics'][:];radial=json.loads(h.attrs['metadata_json'])
 with h5py.File(trajectory) as h:
  t=h['time_s'][:];x=h['positions_m'][:];ids=h['particle_id'][:];v=np.diff(x,axis=0)/np.diff(t)[:,None,None];sample={k:h['force_samples/'+k][:] for k in h['force_samples']};flow=str(h.attrs['flow']);schedule=json.loads(h.attrs['schedule_json'])
 with h5py.File(cache) as h:meta=json.loads(h.attrs['metadata_json']);central=h['pressure_basis'][:,20,20,20];axes=[h['axes/'+a+'_m'][:] for a in 'xyz']
 with h5py.File(kwave) as h:
  km=json.loads(h.attrs['metadata_json']);sh=h['P_complex'].shape;mid=tuple(n//2 for n in sh[1:]);kp=h['P_complex'][(slice(None),)+mid];kg=np.array([h['gradP_complex_'+a][(slice(None),)+mid] for a in 'xyz']).T
 np.savez_compressed(O/'MEVO_ACTUAL_EVIDENCE.npz',trajectory_time_s=t,particle_id=ids,positions_m=x,finite_difference_velocity_m_s=v,**{'force_sample_'+k:v for k,v in sample.items()},response_lower=lo,response_upper=hi,response_valid=valid,response_harmonics_m3=harm,current_proxy_center_basis=central,historical_kwave_center_basis=kp,historical_kwave_center_grad_basis_m_inv=kg)
 files=[response,trajectory,cache,kwave]+[EXEC/p for p in ['src/mevo_mb/fast_sim/models.py','src/mevo_mb/fast_sim/water_b2.py','src/mevo_mb/fast_sim/water_sequential.py','src/mevo_mb/bubbles/apecss.py','src/mevo_mb/bubbles/adaptive_table.py','src/mevo_mb/control/scheduled_pzt.py']]
 z=dict(status='ACTUAL_ASSETS_READ',execution_commit='14731fd56f0477c89cbbdea5ae1f13e8e38e3f3c',assets={str(p):dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},current_acoustics=dict(model=meta['model'],frequency_Hz=meta['frequency_Hz'],source_amplitude_meaning=meta['source_amplitude_meaning'],harmonic_convention=meta['harmonic_convention'],source_count=len(central),spatial_bounds_m=[[a[0],a[-1]] for a in axes],sampled_center_amplitude_per_source=np.abs(central),sampled_center_phase_rad=np.angle(central)),historical_kwave=dict(status='ACTUAL_SOLVE',scope='Historical F1 geometry, NOT latest F10S geometry or new streaming validation',frequency_Hz=km['frequency_Hz'],source_ids=km['source_ids'],source_reference_Pa=km['basis_source_amplitude_Pa'],grid_shape=sh,harmonic_convention='Re[P exp(+i omega t)]',simulation_versions=[q['k_wave_python_version'] for q in km['simulations']],sampled_center_amplitudes=np.abs(kp),sampled_center_phase_rad=np.angle(kp)),apecss=dict(metadata=radial,valid_leaf_count=int(valid.sum()),all_leaf_count=len(valid),outer_bounds=np.c_[lo.min(0),hi.max(0)],cell_policy='Point must lie in a valid adaptive leaf; bounding box alone is insufficient',parameter_order=['radius_m','frequency_Hz','pressure_Pa'],calibrated=False),trajectory=dict(duration_s=float(t[-1]),particles=len(ids),radius_m=1.78e-6,flow=flow,max_particle_speed_m_s=float(np.linalg.norm(v,axis=2).max()),center_of_mass_displacement_m=x[-1].mean(0)-x[0].mean(0),initial_center_span_m=np.ptp(x[0],axis=0),sampled_pressure_range_Pa=[float(np.abs(sample['pressure_Pa']).min()),float(np.abs(sample['pressure_Pa']).max())],sampled_primary_force_max_N=float(np.linalg.norm(sample['primary_force_N'],axis=2).max()),sampled_secondary_force_max_N=float(np.linalg.norm(sample['secondary_force_N'],axis=2).max()),schedule=schedule),streaming_transfer=dict(status='NOT_YET_EVALUABLE',reason='No validated mapping from acoustic phasor/radial harmonics/Bjerknes center motion to mean fluid velocity gradients and clot traction; F10S flow is STRICT_ZERO',background_flow_interface='Older FlowField interfaces support prescribed zero/affine/Y CFD, but F10S water_b2 explicitly rejects flow_enabled; none solves induced streaming',requires=['near-wall flow solver with bubble/cluster forcing','local geometry/gap and hydrodynamic boundary conditions','startup and sustained mean-flow validation','acoustic-pressure-to-traction inverse relation']))
 dump(O/'MEVO_ASSET_AUDIT.json',z);return z
def classify(req,audit=None):
 audit=audit or json.loads((O/'MEVO_ASSET_AUDIT.json').read_text());data=load_npz(O/'MEVO_ACTUAL_EVIDENCE.npz');details={}
 re=1000*abs(req['drive_speed_m_s'])*req.get('gap_m',25e-6)/.001
 details['fluid']=dict(status='CONDITIONALLY_PLAUSIBLE' if re<=.1 else 'OUTSIDE_VALIDATED_RANGE',Re=re,reason='Prescribed planar no-slip mean flow only; curved V3 patch error unvalidated')
 end=req['hold_until_s']+.4e-3;details['time']=dict(status='CONDITIONALLY_PLAUSIBLE' if end<=audit['trajectory']['duration_s'] else 'OUTSIDE_VALIDATED_RANGE',required_end_s=end,available_center_trajectory_s=audit['trajectory']['duration_s'],reason='Time coverage only; not proof of mean-flow formation')
 point=[req.get(k) for k in ['bubble_radius_m','frequency_Hz','acoustic_pressure_Pa']]
 if any(v is None for v in point):details['radial']=dict(status='NOT_YET_EVALUABLE',reason='Inverse traction-to-acoustic pressure/radius/frequency relation absent; required acoustic state not inferred from mechanical shear')
 else:
  point=np.array(point);ok=np.any(data['response_valid']&np.all(data['response_lower']<=point,axis=1)&np.all(point<=data['response_upper'],axis=1));details['radial']=dict(status='SUPPORTED_BY_AVAILABLE_MODEL' if ok else 'OUTSIDE_VALIDATED_RANGE',requested_point=point,scope='Radial table interpolation support only, not experimental material or streaming calibration')
 details['scale']=dict(status='NOT_YET_EVALUABLE',required_gap_m=req.get('gap_m',25e-6),existing_cluster_center_span_m=audit['trajectory']['initial_center_span_m'],reason='Cluster dimensions cannot establish a no-slip upper plane or mean flow layer')
 details['traction_realization']=dict(status='NOT_YET_EVALUABLE',reason=audit['streaming_transfer']['reason'])
 status='OUTSIDE_VALIDATED_RANGE' if any(q['status']=='OUTSIDE_VALIDATED_RANGE' for q in details.values()) else 'NOT_YET_EVALUABLE'
 return dict(status=status,requirement=req,dimensions=details)
def main():
 a=inspect_assets();results=[]
 for p in sorted((O/'trials').glob('*/SUMMARY.json')):
  z=json.loads(p.read_text())
  if z['status']=='COMPLETED':results.append(dict(trial=p.parent.name,**classify(z['params'],a)))
 # Interface branches use measured/table points, never fake optimization data.
 p=a['trajectory']['sampled_pressure_range_Pa'][0];demo=dict(drive_speed_m_s=.004,hold_until_s=.0008,bubble_radius_m=1.78e-6,frequency_Hz=2.5e6,acoustic_pressure_Pa=p)
 interface=dict(in_range_radial=classify(demo,a),out_of_range=classify(dict(demo,bubble_radius_m=3e-6),a))
 assert interface['in_range_radial']['dimensions']['radial']['status']=='SUPPORTED_BY_AVAILABLE_MODEL';assert interface['out_of_range']['status']=='OUTSIDE_VALIDATED_RANGE'
 dump(O/'MEVO_FEASIBILITY.json',dict(status='NOT_YET_EVALUABLE',trials=results,interface_branch_checks=interface,scope='No streaming achievement asserted. Available acoustic/radial/center-motion states read and bounded; traction conversion absent.'))
if __name__=='__main__':main()
