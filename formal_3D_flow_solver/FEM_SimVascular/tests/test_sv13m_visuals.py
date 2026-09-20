from sv13m_support import *
from sv_validation.sv13m import *
def test_required_artifacts_and_not_run_visibility():
 d=accepted('visuals')
 from PIL import Image
 for f in d['figures']:
  image=Image.open(ROOT/f['path']);image.verify()
  assert f['title'] and f['data_source']
 assert len(d['figures'])==(10 if read('svmp_gpu_smoke')['status']=='PASS' else 5)
