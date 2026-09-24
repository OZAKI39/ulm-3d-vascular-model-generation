import importlib.util
from pathlib import Path
from PIL import Image
import json

def test_all_figures_and_raw_sources(report):
 files=sorted((report/'figures').glob('*.png')); assert len(files)==15
 for i,path in enumerate(files):
  assert path.name.startswith(f'{i:02d}_')
  with Image.open(path) as image: assert image.width>=1500 and image.height>=700 and image.getextrema()
  metadata=json.loads((report/'data'/(path.stem+'_plot.json')).read_text())
  for raw in metadata['source_files']: assert (report/'data'/raw).exists()

def test_plotter_can_regenerate_from_saved_data(repo,report,tmp_path):
 import shutil
 folder=tmp_path/'plots'; (folder/'data').mkdir(parents=True); (folder/'figures').mkdir()
 for name in ['05_orientation.csv','05_orientation_statistics.json']: shutil.copyfile(report/'data'/name,folder/'data'/name)
 path=repo/'particle_3d/scripts/plot_particle7.py'; spec=importlib.util.spec_from_file_location('p7_plotter_test',path); module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
 module.main(folder,[5]); assert len(list((folder/'figures').glob('*.png')))==1
