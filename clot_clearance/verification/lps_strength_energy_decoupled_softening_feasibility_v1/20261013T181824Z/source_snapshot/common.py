from pathlib import Path
import json,hashlib,numpy as np
C=json.loads(Path(__file__).with_name("config.json").read_text())
R=Path(C["server"]); O=R/"verification"/C["task"]/C["run"]; F=R/"results"/C["task"]/C["run"]
OLD=R/"verification/lps_native_state_identity_and_fracture_bugfix_v1/20261011T203015Z"
DX=50e-6; K=1e5/3; G=940.; S0=.15; DELTA=150.5e-6
EXPECTED="00b246ebe4224e648c27b74afb64fc22b9cea4ba6ffb23f7c40d96c37ca22ae9"
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""): h.update(b)
 return h.hexdigest()
def dump(p,z):
 def cv(v):
  if isinstance(v,np.ndarray):return v.tolist()
  if isinstance(v,np.generic):return v.item()
  if isinstance(v,Path):return str(v)
  raise TypeError(type(v).__name__)
 Path(p).write_text(json.dumps(z,indent=2,ensure_ascii=False,default=cv,allow_nan=False))
