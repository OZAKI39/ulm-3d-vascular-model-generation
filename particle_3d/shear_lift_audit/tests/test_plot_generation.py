from pathlib import Path
from PIL import Image
import json
ROOT=Path(__file__).resolve().parents[1]

def test_required_figures_readable_and_nonempty():
    names=['00_stage_scope','01_local_shear_rate_map','02_slip_velocity_distribution','03_dimensionless_regime_map','04_candidate_lift_magnitude_distribution','05_lift_vs_drag','06_lift_to_drag_along_representative_trajectories','07_lift_vs_wall_lubrication','08_validity_by_wall_distance','09_branch_point_audit','10_radius_sensitivity','11_scientific_limitations']
    for n in names:
        p=ROOT/'figures'/(n+'.png')
        assert p.stat().st_size>12000
        with Image.open(p) as im:
            assert im.width>=1400 and im.height>=800
            assert len(im.convert('RGB').getcolors(im.width*im.height))>100
