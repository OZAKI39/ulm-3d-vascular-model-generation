#!/usr/bin/env python3
"""Read old remote FEM environment fingerprints, with outputs in the SV0 project only."""
import hashlib
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,now,write_json
old=ROOT.parent/'formal_3D_flow_solver_FEM_lzy';assert old!=ROOT
before=json.loads((old/'outputs/stage03/cleanup/remote.json').read_text())
files={str(p.relative_to(old)):sha256(p) for p in (old/'remote/.env/conda-meta').glob('*.json')}
files['remote/.env/bin/python']=sha256(old/'remote/.env/bin/python')
h=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest()
result={'timestamp':now(),'status':'PASS' if h==before['fem_fingerprint_sha256'] else 'FAIL','file_count':len(files),
    'original_fingerprint_sha256':before['fem_fingerprint_sha256'],'current_fingerprint_sha256':h,
    'old_fem_root':str(old),'old_fem_environment_used':False,'python_or_conda_installations_in_old_environment':False,
    'scope':'Old FEM Python executable and all conda package records compared with its own Stage 3 cleanup fingerprint'}
write_json(ROOT/'outputs/sv0/environment/old_remote_fem_preservation.json',result);print(json.dumps(result,indent=2));assert result['status']=='PASS'
