from pathlib import Path
import json,h5py
R=Path(__file__).resolve().parents[1]
with h5py.File(R/'remote_raw/RAW_performance_beta0.5.h5') as f:
 g=next(g for g in f.values() if g.attrs['N_wall']==1000);peak=max(q.attrs['peak_rss_bytes'] for q in g.values() if isinstance(q,h5py.Group))
peak=int(peak)
est=1.1*4*peak+.5*2**30
record=dict(status='FROZEN_BEFORE_N2000_TEST',reason='Replace initial six-dense-array estimate with measured quadratic peak-memory scaling plus 10 percent and 0.5 GiB margin; physical and main resource caps unchanged',measured_N1000_peak_bytes=peak,N2000_estimated_peak_bytes=est,N2000_estimated_peak_GiB=est/2**30,unchanged_cap_GiB=20,estimated_wall_seconds='4x assembly and 8x factor: approximately 80 seconds; unchanged 180-second bound',initial_N2000_preflight_rejection_retained=True)
assert est<20*2**30
(R/'contracts/PERFORMANCE_N2000_MEASURED_MEMORY_PREFLIGHT.json').write_text(json.dumps(record,indent=2)+'\n');print(record)
