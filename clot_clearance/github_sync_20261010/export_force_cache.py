from pathlib import Path
import json,hashlib
import numpy as np
from pd_clot.lammps_fragmentation.point_cloud_delivery import Data,OUT,RESULTS
from analysis.y_clot_whole_vessel_4k_v1.render import Forces
P=Path('/workspace/y_clot');dest=P/'github_sync_20261010';d=Data(OUT);f=Forces(d);frames=d.schedule();rows=[];hashes=[]
for r in frames:
 X,attached,_,_=d.state(r);row=f.state(r,X,attached);rows.append(row);hashes.append(hashlib.sha256(X.astype('<f8').tobytes()).hexdigest())
keys=['load_N','load_origin_m','drag_N','drag_origin_m','drag_active','valid_CFD','representative_atom_id']
arrays={k:np.asarray([r[k] for r in rows]) for k in keys};arrays.update(frame=np.arange(len(frames)),time_s=np.array([r['time_s'] for r in frames]),xyz_sha256=np.asarray(hashes))
np.savez_compressed(dest/'DISPLAY_FORCE_CACHE.npz',**arrays)
sha=lambda p:hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
meta=dict(scope='Read-only extraction of existing display forces at all 501 accepted saved frames; no new simulation, rounding, interpolation, or replacement of full per-node force audit.',frames=501,representative_atom_id=f.atom,source_force_audit_sha256=sha(f.path),source_frame_manifest_sha256=sha(RESULTS/'frame_manifests/y_clot_point_cloud_early_flow_coupling.json'),cache_sha256=sha(dest/'DISPLAY_FORCE_CACHE.npz'),max_net_force_history_error_N=f.max_sum_error,omitted_full_audit_bytes=f.path.stat().st_size)
(dest/'DISPLAY_FORCE_CACHE_PROVENANCE.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta))
