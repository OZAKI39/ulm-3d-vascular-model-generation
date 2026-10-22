import json,shutil
from pathlib import Path
import h5py,numpy as np
from scripts.stage4.common import ROOT,new_run,finish,sha
out=new_run('performance/secondary_numerical_diagnostic',dict(scope='Closed short trials and copied partial frames; not full secondary branch outcomes'))
shutil.copyfile(__file__,out/'diagnostic_script.py')
trial=json.loads((ROOT/'results/stage4/performance/component_contact_shadow/STATE.json').read_text())
coordinates=[];aggregates=[]
for key in ['baseline','candidate']:
    path=ROOT/trial[key]['run']
    with h5py.File(path/'particles/trajectories.h5') as h:
        assert h.attrs['complete'];coordinates.append(h['positions_m'][-1])
    m=json.loads((path/'metrics.json').read_text())['rows'][-1]
    aggregates.append(dict(run=trial[key]['run'],metrics_sha256=sha(path/'metrics.json'),values={k:m[k] for k in ['time_s','mean_nearest_neighbor_m','cluster_count','largest_cluster_fraction','near_contact_fraction']}))
distance=np.linalg.norm(coordinates[0]-coordinates[1],axis=1)
snapshots=json.loads((ROOT/'results/stage4/performance/secondary_contact_snapshot/STATE.json').read_text())
reynolds=[]
for row in snapshots['rows']:
    identity=json.loads(((ROOT/row['source']).parents[1]/'RUN_IDENTITY.json').read_text())
    medium=identity['config']['medium']
    with h5py.File(ROOT/row['snapshot']) as h:
        re=2*medium['density_kg_m3']*h['radii_m'][:][None,:]*np.linalg.norm(h['total_velocity_m_s'][:]-h['flow_velocity_m_s'][:],axis=-1)/medium['viscosity_Pa_s']
        index=np.unravel_index(np.argmax(re),re.shape)
        reynolds.append(dict(snapshot=row['snapshot'],snapshot_sha256=sha(ROOT/row['snapshot']),source_completed=False,
            maximum_saved_constrained_slip_Re=float(re[index]),time_s=float(h['time_s'][index[0]]),particle_id=int(h['particle_id'][index[1]]),
            samples_at_least_1=int((re>=1).sum()),samples=int(re.size)))
finish(out,dict(status='WARNING',shadow_source=trial['run'],shadow_queries=trial['shadow_global'],
    last_common_time_s=.04,position_error=dict(maximum_m=float(distance.max()),rms_m=float(np.sqrt(np.mean(distance*distance))),median_m=float(np.median(distance))),
    aggregation_comparison=aggregates,partial_snapshot_Re=reynolds,component_solver_deployed=False,
    scope='A tested secondary-ON trajectory is sensitive to roundoff-level contact-solver changes. This short test does not measure final branch-fraction sensitivity. Saved partial-frame Re violates small-Re Stokes applicability; sample counts are not time fractions or unique-particle counts. No new physics or corrected force was introduced.'))
