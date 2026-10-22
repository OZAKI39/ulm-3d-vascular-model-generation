import sys,importlib.metadata as md,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from vcommon import OUT,dump,identity,sha
sources={'torch':'https://github.com/pytorch/pytorch','torchvision':'https://github.com/pytorch/vision','numpy':'https://github.com/numpy/numpy','scipy':'https://github.com/scipy/scipy','opencv-python-headless':'https://github.com/opencv/opencv-python','OpenPIV':'https://github.com/OpenPIV/openpiv-python','pycpd':'https://github.com/siavashk/pycpd','norfair':'https://github.com/tryolabs/norfair','av':'https://github.com/PyAV-Org/PyAV','matplotlib':'https://github.com/matplotlib/matplotlib','jsonschema':'https://github.com/python-jsonschema/jsonschema','pyzmq':'https://github.com/zeromq/pyzmq','psutil':'https://github.com/giampaolo/psutil','h5py':'https://github.com/h5py/h5py'}
rows=[];dest=OUT/'evidence/licenses';dest.mkdir(parents=True,exist_ok=True)
for name in sources:
 d=md.distribution(name);files=[]
 for file in d.files or []:
  if any(s in str(file).lower() for s in ['.dist-info/licenses/','.dist-info/license','.dist-info/copying']):
   p=d.locate_file(file)
   if p.is_file():
    target=dest/(name+'__'+str(file).replace('/','__'));target.write_bytes(p.read_bytes());files.append({'file':str(target.relative_to(OUT)),'sha256':sha(target)})
 rows.append(dict(package=name,version=d.version,official_source=sources[name],license_expression=d.metadata.get('License-Expression') or d.metadata.get('License'),license_files=files,git_commit='NOT_AVAILABLE_FROM_INSTALLED_WHEEL',tag='release version recorded; source tag not independently verified',used_for_formal_results=True,weights='randomly initialized F10V models; no pretrained weight download'))
import torch,torchvision
for r in rows:
 if r['package']=='torch':r['git_commit']=torch.version.git_version
 if r['package']=='torchvision':r['git_commit']=getattr(torchvision.version,'git_version','NOT_AVAILABLE')
dump(OUT/'evidence/dependencies.json',dict(identity=identity(),packages=rows,torch_cuda=torch.version.cuda,python=sys.version,optional={'NeuFlow_v2':{'status':'NOT_RUN','official_source':'https://github.com/neufieldrobotics/NeuFlow_v2','reason':'Core Direct6/GeoMotion, uncertainty and real causal replay take priority within bounded campaign; no weights used'},'SynthPix':{'status':'NOT_RUN','official_source':'https://github.com/antonioterpin/synthpix','reason':'Official principal interface is image pairs; continuous triple/identity adapter and JAX environment not validated within this core campaign'},'DeepTrack2':{'status':'NOT_RUN','reason':'Optional; no need to duplicate real-patch generator'}}))
