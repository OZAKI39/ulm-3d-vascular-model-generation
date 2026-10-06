from pathlib import Path
import csv,json,hashlib,numpy as np
b=Path('/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin');d=b/'population_v1/analysis/stabilization_local_conservation_review_20261006';archive=np.load(d/'STABILIZATION_TERMS.npz');z={k:archive[k] for k in archive.files};count=0
with (d/'STABILIZATION_TERMS.csv').open(newline='') as f:
 r=csv.reader(f);header=next(r);assert header[11:]==list(z['columns'])
 for i,row in enumerate(r):
  assert row[0]==str(z['field'][i]) and int(row[1])==int(z['pid'][i]) and row[2]==str(z['kind'][i]) and int(row[3])==i
  assert int(row[4])==int(z['cell'][i]) and int(row[5])==int(z['phase'][i]) and row[10]==str(z['valid'][i])
  expected=np.r_[z['x_m'][i],z['arc_m'][i],z['values'][i]];actual=np.array([float(x) for x in row[6:10]+row[11:]])
  assert np.array_equal(expected,actual,equal_nan=True),i
  count+=1
assert count==len(z['cell'])
out=dict(status='PASS',all_rows_compared=count,CSV_bytes=(d/'STABILIZATION_TERMS.csv').stat().st_size,NPZ_bytes=(d/'STABILIZATION_TERMS.npz').stat().st_size,CSV_sha256=hashlib.sha256((d/'STABILIZATION_TERMS.csv').read_bytes()).hexdigest(),NPZ_sha256=hashlib.sha256((d/'STABILIZATION_TERMS.npz').read_bytes()).hexdigest(),reason='Full CSV duplicates all NPZ columns and per-row identity; keep lossless NPZ, source exporter and summary. CSV remains local/server.')
Path('/tmp/brava_github_sync_20261006/CSV_NPZ_EQUIVALENCE.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
