"""Verify exact saved states and render with the archived display-force cache."""
import argparse,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from pd_clot.lammps_fragmentation.point_cloud_delivery import Data,OUT,RESULTS
from analysis.y_clot_whole_vessel_4k_v1.render import Scene,decorate

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);parser.add_argument('--verify-only',action='store_true');parser.add_argument('--frame',type=int,default=52);args=parser.parse_args()
 args.output.mkdir(parents=True,exist_ok=False)
 here=Path(__file__).resolve().parent;path=here/'DISPLAY_FORCE_CACHE.npz';provenance=json.loads((here/'DISPLAY_FORCE_CACHE_PROVENANCE.json').read_text());assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==provenance['cache_sha256']
 data=Data(OUT);schedule=data.schedule();manifest=RESULTS/'frame_manifests/y_clot_point_cloud_early_flow_coupling.json';assert hashlib.file_digest(manifest.open('rb'),'sha256').hexdigest()==provenance['source_frame_manifest_sha256']
 recorded=json.loads(manifest.read_text())['frames'];cache=np.load(path,allow_pickle=False);assert len(schedule)==len(recorded)==len(cache['frame'])==501;assert 0<=args.frame<len(schedule)
 snapshot=None
 for r,ref in zip(schedule,recorded):
  X,a,b,c=data.state(r);k=r['frame'];xyz=hashlib.sha256(X.astype('<f8').tobytes()).hexdigest();assert xyz==ref['xyz_sha256']==cache['xyz_sha256'][k]
  assert hashlib.sha256(c.astype('<i8').tobytes()).hexdigest()==ref['component_id_sha256'];assert int(a.sum())==ref['attached'] and int(b.sum())==ref['surviving_bonds']
  assert cache['time_s'][k]==r['time_s']
  if r['time_s']>=.0011605:assert not cache['load_N'][k].any()
  if k==args.frame:snapshot=(r,X,a)
 report=dict(status='PASS',saved_frames_verified=501,native_assets_verified=len(data.sources.files),full_precision_current_trajectory=True,new_scientific_solves=0,force_cache_sha256=provenance['cache_sha256'])
 if not args.verify_only:
  r,X,a=snapshot;k=r['frame'];force={key:cache[key][k] for key in ['load_N','load_origin_m','drag_N','drag_origin_m','drag_active','valid_CFD','representative_atom_id']}
  scene=Scene(data);im=decorate(scene,scene.update(X,a,force),r,force);im.save(args.output/'preview_4k.png');report.update(rendered_frame=k,resolution=list(im.size),x_axis_angle_deg=scene.x_angle_deg);scene.p.close()
 (args.output/'REPLAY_VERIFICATION.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
if __name__=='__main__':main()
