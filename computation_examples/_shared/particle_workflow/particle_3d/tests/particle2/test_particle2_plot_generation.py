import runpy
import shutil
from pathlib import Path
from PIL import Image


def test_synthetic_figures_regenerate_from_permanent_sources(tmp_path,p2_repo,p2_data):
    functions=runpy.run_path(str(p2_repo/"particle_3d/scripts/generate_particle2_report.py"))
    data=tmp_path/"data";figures=tmp_path/"figures";data.mkdir();figures.mkdir()
    for name in ["04_static.csv","04_static_metrics.json","05_rotation.csv","05_rotation_metrics.json","selected_geometries.json"]:
        shutil.copyfile(p2_data/name,data/name)
    # runpy functions keep their own shared globals, separate from returned map.
    scope=functions["stage4"].__globals__
    scope.update(REPORT=tmp_path,DATA=data,FIGURES=figures)
    functions["stage4"]();functions["stage5"]()
    for path in figures.glob("*.png"):
        with Image.open(path) as im:assert im.width>1500 and im.height>900
    assert len(list(figures.glob("*.png")))==2
