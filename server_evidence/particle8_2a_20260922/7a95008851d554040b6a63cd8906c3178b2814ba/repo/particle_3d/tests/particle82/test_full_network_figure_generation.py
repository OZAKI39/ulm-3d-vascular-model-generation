import json
from PIL import Image

def test_all_twelve_figures_bound_to_same_scene(saved_scene):
    s=saved_scene;files=sorted((s.root/'figures').glob('*.png'))
    assert [p.name[:2] for p in files]==[f'{i:02d}' for i in range(12)]
    for p in files:
        im=Image.open(p);im.load();assert im.width>=1600 and im.height>=1000
        receipt=json.loads((s.root/'figure_sources'/(p.stem+'.json')).read_text())
        assert receipt['source_scene_sha256']==s.sha256
        assert any('Frozen flow full-vessel' in line for line in receipt['labels'])
    expected='08_full_network_ulm_reconstruction.png' if s.catalog['all_outlets_observed'] else '08_observed_network_coverage.png'
    assert (s.root/'figures'/expected).exists()
