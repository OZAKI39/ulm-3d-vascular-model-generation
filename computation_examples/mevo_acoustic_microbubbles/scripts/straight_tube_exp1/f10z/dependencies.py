"""Capture installed distribution identities rather than guess upstream commits."""
import importlib.metadata as meta,hashlib,json,sys
from pathlib import Path
USE={'PyMuPDF':'PDF page decode, text and raster review','numpy':'Numerical arrays and linear algebra','scipy':'Distances, sparse graph, original linear algebra','h5py':'Raw and exported HDF5','osqp':'Available reference backend; production active_set unchanged','warp-lang':'Existing GPU field and wall kernels','PyYAML':'Existing configuration loading','matplotlib':'Data-driven PNG/PDF figures','Pillow':'Saved-state video rasterization','av':'H264 encoding and complete PTS/decode validation','pytest':'Meaningful geometry/interface regression tests','scikit-learn':'Installed candidate; not used for fitting this delivery','psutil':'Installed locally; final remote monitoring uses stdlib'}
def collect():
 rows=[]
 for name,use in USE.items():
  try:
   d=meta.distribution(name);m=d.metadata;files=d.files or [];licenses=[]
   for f in files:
    if 'license' in str(f).lower() or 'copying' in str(f).lower():
     p=Path(d.locate_file(f))
     if p.is_file():licenses.append(dict(path=str(f),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
   record=next((f for f in files if str(f).endswith('.dist-info/RECORD')),None)
   rows.append(dict(name=name,version=d.version,license_expression=m.get('License-Expression'),license_metadata=m.get('License'),license_classifiers=m.get_all('Classifier',[]) and [s for s in m.get_all('Classifier',[]) if 'License' in s],license_files=licenses,project_urls=m.get_all('Project-URL',[]),git_commit='NOT_EMBEDDED_IN_INSTALLED_DISTRIBUTION',record_sha256=hashlib.sha256(Path(d.locate_file(record)).read_bytes()).hexdigest() if record else None,use=use))
  except meta.PackageNotFoundError:rows.append(dict(name=name,status='NOT_INSTALLED',use=use))
 return dict(python=sys.version,prefix=sys.prefix,packages=rows,compatibility='Actual campaign/tests/render logs are the compatibility evidence; no unsupported version claim',conditional_not_used=['PyTorch world model','pytorch_mppi','botorch','Gymnasium','xarray'])
if __name__=='__main__':print(json.dumps(collect(),indent=2))

