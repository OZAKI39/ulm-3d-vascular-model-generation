"""Read only a stable byte copy, never open a writer's HDF5 through h5py."""
import json,shutil,time
from pathlib import Path
import h5py,numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from mevo_mb.bubbles.contact import contact_candidates,constrain_pair_velocities
from scripts.stage4.profile_contact_components import component_projection
from scripts.stage4.common import ROOT,new_run,write,finish,sha
out=new_run('performance/secondary_contact_snapshot',dict(scope='Partial copied frames for runtime diagnosis only; not a completed run or restart checkpoint'),['statistics'])
shutil.copyfile(__file__,out/'diagnostic_script.py')
base=ROOT/'results/stage4/primary_secondary/20261003T115216806643Z'
rows=[]
for child in sorted(base.glob('*Z_B*_seed*')):
    source=child/'particles/trajectories.h5'
    destination=out/(child.name+'.h5')
    def identity(stat):return [stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns]
    for attempt in range(5):
        before=identity(source.stat())
        shutil.copyfile(source,destination)
        after=identity(source.stat())
        if before!=after:continue
        try:
            with h5py.File(destination) as h:
                count=len(h['time_s'])
                assert all(len(h[k])==count for k in h if k not in ['particle_id','radii_m'])
                assert np.all(np.diff(h['time_s'][:])>0)
                x=h['positions_m'][-1];r=h['radii_m'][:];v=h['free_total_velocity_m_s'][-1]
                assert h['active'][-1].all() and np.isfinite(x).all() and np.isfinite(v).all()
                saved_time=float(h['time_s'][-1]);complete=bool(h.attrs['complete'])
            break
        except (OSError,AssertionError):continue
    else:raise RuntimeError('Could not obtain a consistent copy of saved frames')
    pairs=contact_candidates(x,r,1e-11,'ckdtree')
    graph=coo_matrix((np.ones(len(pairs)),(pairs[:,0],pairs[:,1])),shape=(len(x),len(x))).tocsr()
    _,labels=connected_components(graph,directed=False)
    sizes=np.bincount(labels)
    timings=[];values=[]
    for method in [constrain_pair_velocities,component_projection]:
        durations=[]
        for _ in range(10):
            start=time.perf_counter();value=method(x,r,v,1e-11,3000,'ckdtree');durations.append(time.perf_counter()-start)
        timings.append(float(np.median(durations)));values.append(value)
    row=dict(source=str(source.relative_to(ROOT)),source_complete=complete,saved_time_s=saved_time,
        stable_source_stat=before,copy_sha256=sha(destination),snapshot=str(destination.relative_to(ROOT)),
        contact_pairs=len(pairs),component_sizes=sorted(sizes[sizes>1].tolist(),reverse=True),
        global_component_median_s=timings,maximum_velocity_difference_m_s=float(np.max(abs(values[0]-values[1]))),
        relative_L2=float(np.linalg.norm(values[0]-values[1])/max(np.linalg.norm(values[0]),1e-30)))
    rows.append(row);write(out/'ROWS.json',rows)
finish(out,dict(status='PASS_WITH_SCOPE',rows=rows,production_deployed=False,
    scope='Stable byte-copy snapshots of partial active results; source files only copied at OS level. HDF5 validation reads copies. Contact benchmarks only; not full trajectories, accepted secondary outcomes or restart checkpoints.'))
