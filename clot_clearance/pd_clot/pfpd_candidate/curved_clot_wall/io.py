"""Saved-array interface for this candidate only."""
from pathlib import Path
import json,hashlib
import numpy as np
from ..kernel import Family
from ..higher_order_rk.mechanics_candidate import Reference,GradientOperator,VariationalBAQP

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,Path):return str(x)
    raise TypeError(type(x).__name__)

def dump(path,data):Path(path).write_text(json.dumps(data,default=plain,indent=2,allow_nan=False)+'\n')
def read(path):return json.loads(Path(path).read_text())
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def array_digest(a):return hashlib.sha256(np.asarray(a).tobytes()).hexdigest()

def model(directory,delta):
    d=Path(directory);mesh=dict(np.load(d/'mesh.npz'));X=mesh['X'];V=mesh['volume']
    pairs=np.load(d/'pairs.npy',mmap_mode='r');s=np.load(d/'source.npy',mmap_mode='r');t=np.load(d/'target.npy',mmap_mode='r');g=np.load(d/'g.npy',mmap_mode='r')
    f=Family.from_arrays(X,V,pairs,delta)
    c=VariationalBAQP(Reference.from_family(f),GradientOperator('PD2',s,t,g))
    return mesh,f,c
