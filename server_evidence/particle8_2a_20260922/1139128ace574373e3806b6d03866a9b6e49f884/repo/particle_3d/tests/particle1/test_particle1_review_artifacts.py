from pathlib import Path
import json
import csv
import numpy as np
import pytest
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]/'reports/particle1'
PAIRS=[('00_particle1_scope_and_provenance.png','00_particle1_scope_and_provenance.json'),
       ('01_single_mb_size_provenance.png','01_single_mb_size_provenance.json'),
       ('02_uniform_flow_trajectory.png','02_uniform_flow_trajectory.csv'),
       ('03_pure_rotation_validation.png','03_pure_rotation_validation.csv'),
       ('04_single_step_dataflow.png','04_single_step_example.json'),
       ('05_real_vessel_single_mb_trajectory.png','05_real_trajectory_dt2.csv'),
       ('06_real_vessel_single_mb_diagnostics.png','06_real_vessel_single_mb_diagnostics.csv'),
       ('07_validation_timestep_comparison.png','07_validation_timestep_comparison.csv')]


@pytest.mark.parametrize('figure,data',PAIRS)
def test_permanent_high_resolution_png_and_plot_source(figure,data):
    with Image.open(ROOT/'figures'/figure) as im:
        assert im.format=='PNG' and im.width>=1500 and im.height>=1000
        assert np.std(np.asarray(im))>10
    p=ROOT/'data'/data
    if p.suffix=='.json': assert json.loads(p.read_text())
    else:
        assert b'\r' not in p.read_bytes()
        with p.open(newline='') as f: assert len(list(csv.DictReader(f)))>=3


def test_all_steps_have_chinese_explanation_and_sample_metadata_binding():
    from particle_3d.audit import sha256
    for index in range(8):
        paths=list(ROOT.glob(f'{index:02d}_*.md'))
        assert len(paths)==1 and '检查说明' in paths[0].read_text()
    meta=json.loads((ROOT/'data/01_single_mb_sample.metadata.json').read_text())
    assert meta['population_sha256']==sha256(ROOT/'data/01_single_mb_sample.csv')
    assert meta['N']==1 and meta['formal_simulation_population'] is False
