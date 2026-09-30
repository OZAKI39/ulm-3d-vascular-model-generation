"""Compare common exposure, not macro indices or raw connected-component count."""
import argparse,csv,json
from pathlib import Path
import numpy as np


def relative(value, reference):
    if reference==0:return 0. if value==0 else None
    return (value-reference)/abs(reference)


def collect(root, output):
    root=Path(root);out=Path(output);out.mkdir(parents=True,exist_ok=True)
    rows=[];runs={}
    for p in sorted(root.glob('*/SUMMARY.json')):
        folder=p.parent;s=json.loads(p.read_text());c=json.loads((folder/'CONFIG.json').read_text())
        cal=json.loads((folder/'CALIBRATION.json').read_text());h=json.loads((folder/'history.json').read_text());f=s['final']
        identity=json.loads((folder/'IDENTITY.json').read_text());z=np.load(folder/'states.npz')
        row=dict(run=folder.name,path=str(folder.resolve()),spacing_mm=c['clot']['particle_spacing_m']*1e3,
                 particle_count=f['particle_count'],initial_bonds=identity['bond_count'],fixed_particle_count=int(z['fixed'].sum()),
                 fixed_volume_fraction=float(z['volume'][z['fixed']].sum()/z['volume'].sum()),
                 calibrated_critical_energy_density_Pa=cal['critical_energy_density_Pa'],Gc_demo_J_m2=cal['Gc_demo_J_m2'],
                 Gc_measured_J_m2=cal['Gc_measured_J_m2'],calibration_relative_error=cal['relative_calibration_error'],
                 first_failure_N=None if s['events']['first_bond_failure'] is None else s['events']['first_bond_failure']['cycles'],
                 first_detachment_N=None if s['events']['first_detachment'] is None else s['events']['first_detachment']['cycles'],
                 first_resolved_detachment_N=None if s['events']['first_resolved_detachment'] is None else s['events']['first_resolved_detachment']['cycles'],
                 DeltaN=c['damage']['DeltaN'],dt_s=c['simulation']['dt_s'],C_E=c['regularized_damage']['C_E'],m_E=c['regularized_damage']['m_E'],
                 streaming_amplitude_m_s=c['streaming']['streaming_velocity_scale_m_s'],configuration_sha256=identity['config_sha256'],
                 elapsed_s=s['elapsed_wall_s'],**{k:f[k] for k in ['represented_cycles','total_volume_m3','attached_volume_fraction','detached_volume_fraction',
                 'resolved_fragment_volume_m3','resolved_fragment_volume_fraction','under_resolved_debris_volume_m3','singleton_volume_m3',
                 'singleton_detached_volume_fraction','lowrank_detached_volume_m3','P_lowrank','largest_resolved_fragment_volume_m3',
                 'damage_dissipation_estimate_J','relative_numerical_energy_residual','mean_damage','maximum_damage','broken_bond_fraction',
                 'maximum_displacement_m','fraction_of_damage_inside_influence_region','fraction_of_broken_bonds_inside_influence_region',
                 'fraction_of_detached_volume_originating_in_influence_region','damage_weighted_distance_from_source_m',
                 'total_clearance_volume_fraction','resolved_fragment_clearance_fraction','under_resolved_debris_clearance_fraction','fragment_resolution_status']})
        rows.append(row);runs[folder.name]=(row,h)
    if not rows:raise ValueError('No completed runs')
    with (out/'study_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    comparisons={};metrics=['attached_volume_fraction','detached_volume_fraction','mean_damage','maximum_damage',
                           'damage_dissipation_estimate_J','maximum_displacement_m','fraction_of_damage_inside_influence_region']
    if 'coarse' in runs:
        ref,rh=runs['coarse'];refmap={h['represented_cycles']:h for h in rh}
        for name in ['cycle_500','cycle_250','medium','fine','dt_half']:
            if name not in runs:continue
            row,hist=runs[name];common=[]
            for h in hist:
                N=h['represented_cycles']
                if N not in refmap:continue
                item={'N':N}
                for metric in metrics:
                    a,b=h[metric],refmap[N][metric]
                    item[metric+'_absolute_difference']=None if a is None or b is None else a-b
                    item[metric+'_relative_difference']=None if a is None or b is None else relative(a,b)
                common.append(item)
            comparisons[name]=dict(reference='coarse',common_exposure=common,
                first_failure_shift_N=None if row['first_failure_N'] is None or ref['first_failure_N'] is None else row['first_failure_N']-ref['first_failure_N'],
                first_detachment_shift_N=None if row['first_detachment_N'] is None or ref['first_detachment_N'] is None else row['first_detachment_N']-ref['first_detachment_N'],
                event_comparison='NOT_OBSERVED_IN_EXPOSURE' if row['first_failure_N'] is None or ref['first_failure_N'] is None else 'OBSERVED')
    reasons=[]
    if not all(k in runs for k in ['coarse','medium','fine']):reasons.append('Incomplete spacing study')
    if all(runs[n][0]['first_detachment_N'] is None for n in ['coarse','medium','fine'] if n in runs):
        reasons.append('No detachment in tested exposure; failure/fragmentation mesh convergence cannot be established')
    for name in ['medium','fine']:
        if name in comparisons:
            for metric in ['mean_damage','damage_dissipation_estimate_J']:
                delta=comparisons[name]['common_exposure'][-1][metric+'_relative_difference']
                if delta is not None and abs(delta)>.1:reasons.append(f'{name}/coarse {metric} differs by {delta:.2%}, exceeding descriptive 10% comparison tolerance')
    result=dict(runs=rows,comparisons=comparisons,mesh_status='NOT_CONVERGED' if reasons else 'CONVERGENCE_CRITERIA_MET_FOR_TESTED_METRICS',
                mesh_reasons=reasons,zero_event_convention='Null = not observed by final N, not cycle zero. Zero detached volumes do not validate fragment size or resolution.',
                energy_scope='Representative mechanical path; no multiplication of external work by cycle jump',
                fixed_slab_rule='Same physical z < origin_z + 0.125 mm; centers at boundary excluded with 1e-10 h roundoff tolerance')
    (out/'COMPARISONS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,default=Path('results/streaming_regularized_demo'))
    p.add_argument('--output',type=Path,default=Path('verification/regularization/comparisons'));a=p.parse_args();r=collect(a.runs,a.output)
    print(json.dumps({'runs':len(r['runs']),'mesh_status':r['mesh_status'],'reasons':r['mesh_reasons']},indent=2))


if __name__=='__main__':main()
