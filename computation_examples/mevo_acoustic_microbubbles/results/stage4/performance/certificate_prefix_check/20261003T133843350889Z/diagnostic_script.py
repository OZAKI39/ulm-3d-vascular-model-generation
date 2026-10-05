import json,shutil
from pathlib import Path
import h5py,numpy as np
from scripts.stage4.common import ROOT,new_run,finish,sha
old=ROOT/'results/stage4/primary_secondary/20261003T115216806643Z'
new=ROOT/'results/stage4/primary_secondary/20261003T132343245938Z'
out=new_run('performance/certificate_prefix_check',dict(reference=str(old.relative_to(ROOT)),candidate=str(new.relative_to(ROOT)),scope='Exact saved-prefix comparison only; both full secondary cases remain incomplete'))
shutil.copyfile(__file__,out/'diagnostic_script.py')
rows=[]
for label in ['B2','B3']:
    a=next(old.glob('*Z_'+label+'_seed*'));b=next(new.glob('*Z_'+label+'_seed*'))
    source=b/'particles/trajectories.h5';snapshot=out/(label+'_candidate_snapshot.h5')
    for _ in range(5):
        before=source.stat();shutil.copyfile(source,snapshot);after=source.stat()
        if (before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=(after.st_size,after.st_mtime_ns,after.st_ctime_ns):continue
        try:
            with h5py.File(snapshot) as h:
                count=len(h['time_s']);assert all(len(h[k])==count for k in h if k not in ['particle_id','radii_m'])
            break
        except (OSError,AssertionError):continue
    else:raise RuntimeError('No consistent saved candidate copy')
    with h5py.File(a/'particles/trajectories.h5') as reference,h5py.File(snapshot) as candidate:
        frames=min(len(reference['time_s']),len(candidate['time_s']))
        equality={k:bool(np.array_equal(reference[k][:] if k in ['particle_id','radii_m'] else reference[k][:frames],candidate[k][:] if k in ['particle_id','radii_m'] else candidate[k][:frames],equal_nan=True)) for k in reference}
        rows.append(dict(case=label,reference=str(a.relative_to(ROOT)),candidate=str(b.relative_to(ROOT)),common_frames=frames,
            last_compared_time_s=float(reference['time_s'][frames-1]),datasets_bitwise_equal=equality,
            candidate_snapshot=str(snapshot.relative_to(ROOT)),candidate_snapshot_sha256=sha(snapshot),
            source_stat_unchanged=[before.st_size,before.st_mtime_ns,before.st_ctime_ns]))
finish(out,dict(status='PASS' if all(all(r['datasets_bitwise_equal'].values()) for r in rows) else 'FAIL',rows=rows,
    scope='Only the common saved prefix was compared bitwise. Active source HDF5 was byte-copied without opening it through HDF5; copies were validated. No full-trajectory or physical validity claim.'))
