from sv13n_support import *
def test_five_requested_development_figures():
 from PIL import Image
 d=accepted('visuals');assert len(d['figures'])==5
 assert {Path(f['path']).name for f in d['figures']}=={'gpu_stack_pipeline.png','cuda_ghost_status.png','gpu_solver_progress.png','gpu_steady_convergence.png','gpu_final_field.png'}
 for f in d['figures']:
  p=ROOT/f['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==f['sha256']
  with Image.open(p) as im:assert im.width>=1200 and im.height>=700
