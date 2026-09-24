#!/usr/bin/env python3
"""Independent source/clock/outlet audit and compact full-sample derived export."""
from pathlib import Path
import sys,csv,gzip
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle81_simulation import OUTPUT,environment,dump,SEED,SONOVUE
from particle_3d.particle81_replay import Scene
from particle_3d.particle8_replay import read,digest
from particle_3d.injection_population import PopulationSource,C_MB,FluxClock,LinearProfile,ConstantMBConcentrationV0


def audit(root=OUTPUT):
    scene=Scene(root);root=Path(root);births=read(root/'data/birth_ledger.json');env=environment()
    if [e['particle_id'] for e in births['events']]!=list(range(1,len(births['events'])+1)):
        raise ValueError('Formal population audit requires the complete contiguous birth ledger')
    source=PopulationSource(SONOVUE,SEED);clock=FluxClock(LinearProfile([0],[env.sampler.Q_m3_s]),ConstantMBConcentrationV0())
    checks=[];field_probes=[];outlet_probes=[];max_clock_error=0.
    provenance=read(root/'data/frozen_provenance.json')['provenance']
    frozen_checks={key:digest(Path(provenance[key+'_path']))==provenance[key+'_sha256'] for key in ['mesh','flow']}
    for role,cap in provenance['boundary_manifest']['boundaries'].items():
        frozen_checks[role]=digest(Path(provenance['fem_root'])/cap['path'])==cap['sha256']
    with gzip.open(root/'data/trajectory_samples.csv.gz','wt',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(['particle_id','physical_time_s',*scene.catalog['original_sample_columns'],'speed_m_s','active_until_terminal'])
        for event in births['events']:
            pid=event['particle_id'];e=scene.entries[pid];a=scene.arrays[pid];meta=read(root/e['metadata_path'])
            original=source.next_mb();q=source.orientation()[0]
            same_size=original['diameter_um']==event['diameter_um'] and original['radius_m']==event['radius_m']
            same_q=np.array_equal(q,event['q']);same_clock=event['birth_time_s']==clock.time_at(pid)
            max_clock_error=max(max_clock_error,abs(C_MB*env.sampler.Q_m3_s*event['birth_time_s']-pid))
            rng=np.random.default_rng(event['position_seed']);proposal,faces=env.sampler.sample(rng)
            # Reproduce the original P7 sampler before admission conditioning.
            candidate=meta['admission_candidates'][0]
            same_proposal=np.array_equal(proposal[0],candidate['position_m']) and int(faces[0])==candidate['triangle_id']
            checks.append(dict(particle_id=pid,original_size_exact=same_size,original_isotropic_q_exact=bool(same_q),
                               deterministic_birth_exact=same_clock,first_flux_proposal_exact=bool(same_proposal)))
            if len(a):
                for k in np.unique(np.linspace(0,max(0,len(a)-2),min(9,len(a)),dtype=int)):
                    inside=env.field.sample(a[k,1:4]).inside_lumen
                    field_probes.append(dict(particle_id=pid,sample_index=int(k),inside_original_fem=bool(inside)))
                speed=np.linalg.norm(a[:,4:7],axis=1)
                for k,row in enumerate(a):writer.writerow([pid,event['birth_time_s']+row[0],*row,speed[k],int(k<len(a)-1)])
            if e['completed']:
                hit=env.classifier.first_event(a[-2,1:4],a[-1,1:4])
                outlet_probes.append(dict(particle_id=pid,outlet=hit.role if hit else None,
                    official_final_segment_match=hit is not None and hit.role==e['exit_outlet']))
    events=read(root/'data/particle8_1_events.json')['events']
    with (root/'data/events.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['particle_id','time_s','kind','order','reason','outlet'],lineterminator='\n')
        writer.writeheader();writer.writerows(events)
    record=dict(all_pass=all(all(v for k,v in c.items() if k!='particle_id') for c in checks) and
                all(p['inside_original_fem'] for p in field_probes) and all(p['official_final_segment_match'] for p in outlet_probes) and all(frozen_checks.values()),
        frozen_original_input_sha256_checks=frozen_checks,
        source_scene_sha256=scene.sha256,original_population_checks=checks,all_population_exact=True,
        maximum_cumulative_count_roundoff=max_clock_error,field_probe_scope='UP_TO_NINE_SAVED_PREEXIT_POINTS_PER_ADMITTED_MB; NOT_ALL_TIMESTEPS',
        field_probes=field_probes,completed_final_outlet_checks=outlet_probes,
        derived_sample_export='data/trajectory_samples.csv.gz',derived_sample_count=scene.catalog['physical_samples'],
        derived_export_sha256=digest(root/'data/trajectory_samples.csv.gz'),
        outlet_semantics='FIRST_CENTER_CROSSING_OFFICIAL_OPEN_CAP; NOT_COMPLETE_FINITE_SPHERE_CLEARANCE',
        failures_are_physiological_capture_probabilities=False,
        explicit_time_representation='Absolute acquisition time plus original per-MB elapsed time; elapsed time preserves substep precision')
    record['all_population_exact']=all(all(v for k,v in c.items() if k!='particle_id') for c in checks)
    dump(root/'data/physics_audit.json',record);print({k:record[k] for k in ['all_pass','maximum_cumulative_count_roundoff','derived_sample_count']},flush=True)
    return record


if __name__=='__main__':audit()
