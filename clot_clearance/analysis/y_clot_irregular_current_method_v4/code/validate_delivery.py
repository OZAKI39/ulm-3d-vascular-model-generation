from common import *
from PIL import Image
import subprocess,shutil
def main():
 result=json.loads((O/'NATIVE_V4_RESULT.json').read_text());assert result['status'] in ['COMPLETED','NUMERICAL_FAILURE']
 freeze=json.loads((O/'V4_IMPLEMENTATION_FREEZE.json').read_text())
 for name,digest in freeze['sources'].items():assert sha(A/name)==digest,name
 assert sha(LIB)==freeze['library_sha256'] and sha(O/'provenance/reference.npz')==freeze['reference_sha256']
 checks={};manifest={}
 for p in sorted((F/'figures').glob('*.png')):
  im=Image.open(p);im.verify();assert min(im.size)>=800;checks[p.name]=dict(size=im.size,format='PNG')
 assert len(checks)>=7
 for path,w,h,fps in [('videos/VIDEO_01_Y_CLOT_V4_FRAGMENTATION.mp4',1920,1080,'24/1'),('rotation/VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_4K.mp4',3840,2160,'30/1'),('rotation/VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_1080P.mp4',1920,1080,'30/1')]:
  p=F/path;q=json.loads(subprocess.check_output(['ffprobe','-v','error','-select_streams','v:0','-show_entries','stream=width,height,nb_frames,r_frame_rate,codec_name','-of','json',str(p)]))['streams'][0]
  assert q['width']==w and q['height']==h and q['r_frame_rate']==fps and int(q['nb_frames'])>=24
  if 'rotation/' in path:assert int(q['nb_frames'])==540
  subprocess.run(['ffmpeg','-v','error','-i',str(p),'-f','null','-'],check=True,capture_output=True);checks[path]=dict(q,full_decode='PASS')
 pv=json.loads((O/'DYNAMICS_VIDEO_PROVENANCE.json').read_text());seen={};last_time=-1.;ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count'])
 for frame in pv['frames']:
  for path,stat in zip(frame['sources'],frame['stats']):
   assert stat['world_Z_horizontal_error_px']<1e-7
   p=R/path
   with np.load(p) as z:
    t=float(z['time']);assert z['x'].dtype==np.float64 and z['d'].dtype==np.float64 and t<=result['accepted_time_s']+1e-14 and t>=last_time
    assert np.array_equal(z['x'][n:],ref['X'][n:]);assert abs(t-stat['time_s'])<1e-14
    last_time=t;seen[path]=dict(time_s=t,step=int(z['step']),sha256=sha(p))
 assert abs(last_time-result['accepted_time_s'])<1e-14
 rot=json.loads((O/'rotation/VISUALIZATION_VALIDATION.json').read_text());assert rot['status']=='PASS' and rot['layout']['max_Z_horizontal_error_px']==0 and rot['layout']['clipping_count']==0 and not rot['layout']['visible_text_collisions']
 load=json.loads((O/'LOAD_VISUALIZATION_AUDIT.json').read_text());assert load['force_replay_max_N']<1e-15 and load['force_times_area_error_N']<1e-15
 dump(O/'VIDEO_STATE_SOURCE_HASHES.json',seen)
 for p in [*F.rglob('*'),*A.glob('*'),D/'REVIEW_Y_CLOT_IRREGULAR_CURRENT_METHOD_V4.md',D/'REPRODUCE.md']:
  if p.is_file() and '__pycache__' not in str(p):manifest[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p))
 dump(O/'ARTIFACT_MANIFEST.json',manifest)
 validation=dict(status='PASS',physics_freeze_unchanged=True,artifact_checks=checks,actual_saved_states_checked=len(seen),full_video_decodes=3,world_Z_up=True,displacement_scale=1,load_replay='PASS',rotation_frames=540)
 dump(O/'DELIVERY_VALIDATION.json',validation);shutil.copy2(O/'DELIVERY_VALIDATION.json',D/'DELIVERY_VALIDATION.json')
 final=json.loads((O/'FINAL_STATUS.json').read_text());final['delivery_validation']='PASS';dump(O/'FINAL_STATUS.json',final);dump(D/'FINAL_STATUS.json',final)
 print('DELIVERY_VALIDATION_PASS',flush=True)
if __name__=='__main__':main()
