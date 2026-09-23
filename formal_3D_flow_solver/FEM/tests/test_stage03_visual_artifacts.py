from stage03_helpers import *
from PIL import Image
from fem3d.audit import sha256

def test_all_eleven_wsl_figures_exist_with_measured_sources_and_full_resolution():
    names=['real_geometry_and_bc','velocity_global','velocity_slices','pressure_global','pressure_sections','flux_balance','outlet_flow_split','residual_cell_flow_check','residual_cell_locations','solver_resource_usage','velocity_gradient_slice']
    m=read('reports/stage03/visualization_manifest.json')
    assert m['generated_on']=='WSL' and set(m['figures'])=={n+'.png' for n in names}
    for name in names:
        path=ROOT/'reports/stage03'/f'{name}.png'
        assert sha256(path)==m['figures'][path.name]['sha256']
        with Image.open(path) as im:assert im.width>=1000 and im.height>=600
