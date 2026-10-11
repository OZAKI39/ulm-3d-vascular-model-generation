"""Inventory the new run, select compact evidence, retain all raw data on server."""
from analysis.dual_track_v1.common import *
import datetime,tarfile

def main():
 roots=[A,B,FA,FB,COMMON,ROOT/'analysis/dual_track_v1'];excluded_names={'SERVER_ARTIFACT_MANIFEST.json','COMPACT_MANIFEST.json','LOCAL_RETRIEVAL_AUDIT.json','compact_results.tar.gz'}
 files=sorted({p for d in roots for p in d.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in excluded_names})
 rows=[];selected=[];omitted=[]
 npz_names={'equilibrium.npz','reference_topology.npz','geometry_arrays.npz','CLOT_SURFACE_TRACTIONS.npz','COMMON_POSITION_COMPARISON.npz'}
 text_ext={'.json','.md','.html','.txt','.csv','.conf','.xml','.inp','.in','.patch','.diff','.py','.cpp','.cxx'}
 for p in files:
  row=dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p));rows.append(row)
  take=(p.suffix in text_ext and row['bytes']<12*1024**2) or p.suffix in {'.png','.svg'} or (p.suffix=='.log' and row['bytes']<6*1024**2) or (p.name in npz_names and row['bytes']<20*1024**2)
  (selected if take else omitted).append(row)
 raw=dict(generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),server_root=str(ROOT),new_run_only=True,file_count=len(rows),total_bytes=sum(v['bytes'] for v in rows),files=rows,self_and_transfer_manifests_excluded_to_avoid_recursive_hashing=True)
 save(COMMON/'SERVER_ARTIFACT_MANIFEST.json',raw)
 p=COMMON/'SERVER_ARTIFACT_MANIFEST.json';selected.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)))
 compact=dict(selected_files=len(selected),selected_bytes=sum(v['bytes'] for v in selected),files=selected,omitted_files=len(omitted),omitted_bytes=sum(v['bytes'] for v in omitted),omitted=omitted,scope='New run only. Scientific arrays retain original precision. Native body.data, full mesh/VTU/restarts, shared libraries and intermediate native geometry remain on server. Source/inputs/config/logs/reports/figures/accepted saved static states/traction arrays retained locally.',archive='compact_results.tar.gz',manifest_self_not_hashed=True)
 save(COMMON/'COMPACT_MANIFEST.json',compact)
 with tarfile.open(COMMON/'compact_results.tar.gz','w:gz',compresslevel=3) as tf:
  for row in selected:tf.add(ROOT/row['path'],arcname=row['path'],recursive=False)
  p=COMMON/'COMPACT_MANIFEST.json';tf.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
 print(json.dumps(dict(raw_files=len(rows),raw_bytes=raw['total_bytes'],selected_files=len(selected),selected_bytes=compact['selected_bytes'],archive_bytes=(COMMON/'compact_results.tar.gz').stat().st_size)))
if __name__=='__main__':main()
