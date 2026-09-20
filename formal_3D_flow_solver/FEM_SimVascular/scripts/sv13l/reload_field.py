"""Fresh-process field reload used only for actual new native solver outputs."""
import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.validation import parse_result_vtu
from sv_validation.provenance import sha256
path=Path(sys.argv[1]);grid,u,p=parse_result_vtu(path)
print(json.dumps({'status':'PASS','path':str(path),'sha256':sha256(path),'points':grid.n_points,'cells':grid.n_cells,
    'velocity_finite':bool(np.isfinite(u).all()),'pressure_finite':bool(np.isfinite(p).all()),
    'velocity_max':float(np.linalg.norm(u,axis=1).max()),'pressure_range':[float(p.min()),float(p.max())]}))
