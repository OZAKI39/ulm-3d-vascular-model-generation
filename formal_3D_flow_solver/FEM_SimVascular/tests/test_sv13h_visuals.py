from sv13h_support import *
from PIL import Image
def test_ten_required_images_have_verified_hashes():
    d=load('visuals');assert len(d['figures'])==10
    names=['cuda_version_strategy','gpu_build_pipeline','cuda_runtime_smoke','petsc_gpu_backend',
        'svmultiphysics_gpu_smoke','gpu_cpu_solution_difference','gpu_residency','gpu_transfer_cost','gpu_memory','cpu_gpu_runtime']
    assert {Path(f['path']).stem for f in d['figures']}==set(names)
    for f in d['figures']:
        p=ROOT/f['path'];assert sha256(p)==f['sha256']
        with Image.open(p) as im:assert im.width>=1200 and im.height>=700;im.verify()
def test_unmeasured_panels_have_no_numeric_bars():
    d=load('visuals')
    unrun=[f for f in d['figures'] if f['label']=='NOT RUN']
    assert len(unrun)==7
    assert all(not f['has_measurements'] and f['numeric_bars']==0 for f in unrun)

