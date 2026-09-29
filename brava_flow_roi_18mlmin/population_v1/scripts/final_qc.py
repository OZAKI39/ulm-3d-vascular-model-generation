"""Final raw radius/ID/hash/exit-endpoint checks and encoded video inspection."""
from pathlib import Path
import argparse,json,csv,gzip,hashlib,subprocess
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--run-name',required=True);args=ap.parse_args()
root=Path(__file__).resolve().parents[1];run=root/'runs'/args.run_name;analysis=root/'analysis'/args.run_name
manifest=json.loads((run/'MANIFEST.json').read_text());identity=manifest['identity'];audit=json.loads((analysis/'INDEPENDENT_AUDIT.json').read_text());summary=json.loads((run/'SUMMARY.json').read_text())
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
config=json.loads((root/'data/environment.json').read_text());input_hashes={}
for key,expected in config['sha256'].items():
 path=Path(config[key] if key in ['mesh','flow'] else config['boundaries'][key]);actual=digest(path);assert actual==expected,key;input_hashes[key]=actual
for name,expected in identity['code_files_sha256'].items():assert digest(root/'source/particle_3d/src/particle_3d/population_v1'/name)==expected,name
with gzip.open(run/'source.jsonl.gz','rt') as f:source=[json.loads(l) for l in f]
accepted=[r for r in source if r['particle_id'] is not None];radius=np.full(1+max(r['particle_id'] for r in accepted),np.nan)
for row in accepted:radius[row['particle_id']]=row['radius_m']
max_radius_error=0.;max_diameter_error=0.
for chunk in manifest['chunks']:
 with np.load(run/chunk['path']) as d:
  max_radius_error=max(max_radius_error,float(np.max(abs(d['radius']-radius[d['particle_id']]))))
  max_diameter_error=max(max_diameter_error,float(np.max(abs(d['diameter']-2*d['radius']))))
assert max_radius_error==max_diameter_error==0.
rows=list(csv.DictReader((analysis/'trajectory_catalog.csv').open()));outside_terminal_ids=[]
for row in rows:
 pid=int(row['particle_id'])
 with np.load(analysis/'tracks'/f'mb_{pid:06d}.npz') as d:
  bad=np.flatnonzero(d['tetra_id']<0)
  if len(bad):
   assert len(bad)==1 and bad[0]==len(d['tetra_id'])-1 and row['status']=='EXITED',(pid,bad)
   outside_terminal_ids.append(pid)
ro=audit['wall_roundoff_m'];wall_ok=audit['independent_wall_max_error_m']<=ro and audit['independent_wall_min_gap_m']>=2e-9-ro
assert wall_ok
assert audit['birth_admission_lower_bound_max_violation_m']<=ro
assert audit['maximum_held_velocity_path_error_m']<=ro
assert audit['independent_outlet_pass'] and audit['source_stop_rule_reproduced']
videos=[]
for pose in ['0','15']:
 path=root/f'candidate_{pose}/animations/microbubble_physical_time.mp4';meta=json.loads((root/f'candidate_{pose}/data/render_manifest.json').read_text())
 info=json.loads(subprocess.check_output(['/usr/bin/ffprobe','-v','error','-count_frames','-show_entries','stream=codec_type,width,height,r_frame_rate,nb_read_frames,duration','-of','json',str(path)],text=True));v=next(s for s in info['streams'] if s['codec_type']=='video')
 assert int(v['nb_read_frames'])==576 and (v['width'],v['height'])==(1920,1080) and v['r_frame_rate']=='24/1',v
 assert not meta['preview_only'] and meta['all_particles_included'] and meta['text_language']=='English'
 assert np.allclose(meta['layout']['projected_size_ratios'],[1.3,1.3],atol=1e-10,rtol=0)
 assert meta['production_git_commit']==identity['git_commit'] and meta['flow_sha256']==identity['flow_sha256']
 videos.append(dict(candidate=pose,path=str(path.relative_to(root)),sha256=digest(path),bytes=path.stat().st_size,stream=v,
  GPU_EGL=any('NVIDIA' in s for s in meta['render_backend']),encoder=meta['encoder'],main_scale=meta['zoom']))
result=dict(data_integrity_pass=True,full_physical_trajectory_success=summary['scientific_integrity_pass'],
 scientific_scope='Data/implementation consistency does not reclassify numerical failures or censoring as successful transits',
 flow_and_boundary_hashes_unchanged=input_hashes,production_code_bytes_match_identity=True,
 maximum_radius_drift_from_original_source_m=max_radius_error,maximum_diameter_vs_2radius_error_m=max_diameter_error,
 interior_tetra_ids_valid=True,terminal_outlet_endpoint_outside_tetra_ids=outside_terminal_ids,
 exact_wall_gap_parity_and_handoff_with_documented_roundoff_pass=wall_ok,
 held_velocity_path_consistency_with_geometry_roundoff_pass=True,videos=videos)
(root/'validation/final_delivery_qc.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
