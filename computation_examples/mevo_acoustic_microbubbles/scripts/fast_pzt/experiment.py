"""One bounded F1 scene; common interface for standalone and batched datasets."""
import argparse
import json
from dataclasses import asdict
from time import perf_counter
import numpy as np
from mevo_mb.fast_sim.models import Physics,ResponseSection,Basis,random_cloud
from mevo_mb.fast_sim.engine import FastSimulator
from mevo_mb.fast_sim.storage import save_run
from .common import BASE,ROOT,inputs,write,source_identity


def run_case(out,config):
    out.mkdir(parents=True,exist_ok=False); process_start=perf_counter()
    physics=Physics(**config.get('physics',{})); response=ResponseSection(inputs()[0],physics)
    basis=Basis(inputs()[1],physics.frequency_Hz) if config.get('basis',False) else None
    positions=np.asarray(config['positions_m']) if 'positions_m' in config else np.array([
        random_cloud(config['count'],physics.radius_m,config['cloud_radius_m'],seed,config.get('center_m',[0,0,0]))
        for seed in config.get('seeds',[config['seed']])])
    identity=dict(source_sha256=source_identity(),execution_commit=(ROOT/'SOURCE_F1_COMMIT.txt').read_text().strip() if (ROOT/'SOURCE_F1_COMMIT.txt').exists() else None,
                  response=response.metadata,basis=None if basis is None else basis.metadata)
    write(out/'CONFIG.json',config); write(out/'IDENTITY.json',identity)
    simulator=FastSimulator(positions,physics,response,basis)
    result=simulator.run(config['duration_s'],config['dt_s'],config['events'],config.get('sample_interval_s',0.0005),
                         config.get('wall_budget_s',180),graphs=config.get('graphs',True),substeps=config.get('substeps',1))
    metrics,events=save_run(out/'trajectory.h5',result,physics,identity,config)
    summary=result['summary']; summary['process_body_wall_s']=perf_counter()-process_start
    summary['configuration']=config; summary['physics_parameters']=asdict(physics)
    summary['formation_time_s']=[]
    for b in range(positions.shape[0] if positions.ndim==3 else 1):
        hits=np.flatnonzero(metrics['largest_size'][:,b]>=config.get('formation_threshold_size',16))
        summary['formation_time_s'].append(float(result['time_s'][hits[0]]) if len(hits) else None)
    protocol=json.loads((ROOT/'configs/fast_pzt/f1_protocol.json').read_text())
    q=summary['quality']; bounds=protocol['quality']
    summary['numerical_quality']='PASS' if summary['status']=='COMPLETE' and q['error_bits']==0 and q['max_accepted_penetration_m']<=bounds['max_accepted_penetration_over_radius']*physics.radius_m and q['max_contact_pass_correction_m']<=bounds['max_contact_pass_correction_over_radius']*physics.radius_m else 'FAIL'
    summary['edge_events_recorded']=len(result.get('edge_event_time_s',[]))
    write(out/'SUMMARY.json',summary)
    print(json.dumps({k:summary[k] for k in ['status','wall_s','process_body_wall_s','numerical_quality','final_largest_size','final_cluster_count','formation_time_s']}),flush=True)
    return summary


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--label',required=True); ap.add_argument('--config',required=True)
    a=ap.parse_args(); run_case(BASE/a.label,json.loads((ROOT/a.config).read_text()))


if __name__=='__main__': main()
