#!/usr/bin/env python3
import json
import shutil
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import inventory,git_state,sha256,write_json,now
OLD=ROOT.parent/'FEM';O=ROOT/'reports/sv0';I=ROOT/'inputs/fem_reference'
assert not (O/'old_fem_baseline.json').exists()
state_before=git_state(OLD);baseline=inventory(OLD);baseline['git']=state_before
write_json(O/'old_fem_baseline.json',baseline)
# Stage 3 carries a JSON execution snapshot tied to the unchanged formal YAML.
snapshot_path=OLD/'inputs/stage03/reference_condition.json';snapshot=json.loads(snapshot_path.read_text());config=snapshot['config']
assert sha256(OLD/'configs/stage03_reference_vascular.yaml')==snapshot['yaml_sha256']
selected=Path(config['mesh']['source']);assert selected.as_posix()=='outputs/stage01_7/selected'
assert sha256(OLD/selected/'mesh/volume_mesh.npz')==config['mesh']['volume_mesh_sha256']
chosen={'volume_mesh.npz':selected/'mesh/volume_mesh.npz','exterior_surface.npz':selected/'surface/tagged_surface_si.npz',
        'port_contract.json':selected/'planar_port_contract_v2.json','baseline_quality.json':selected/'qc/volume_quality.json',
        'selection.json':selected/'metadata/selection.json','stage00_source_contract.json':Path('reports/stage00/source_contract.json'),
        'stage03_reference_vascular.yaml':Path('configs/stage03_reference_vascular.yaml'),
        'stage03_reference_condition.json':Path('inputs/stage03/reference_condition.json'),
        'stage03_pressure_support.json':Path('reports/stage03/singularity_diagnosis.json')}
records={}
for name,relative in chosen.items():
    source=OLD/relative;target=I/name;assert not target.exists();h=sha256(source)
    shutil.copy2(source,target);copied=sha256(target);assert copied==h
    records[name]={'source_path':str(source),'source_sha256':h,'copied_path':str(target.relative_to(ROOT)),
                   'copied_sha256':copied,'size':source.stat().st_size,'copy_time':now()}
assert git_state(OLD)==state_before
write_json(I/'MANIFEST.json',{'timestamp':now(),'status':'PASS','old_fem_root':str(OLD),'source_git':state_before,
                            'files':records,'total_copied_bytes':sum(r['size'] for r in records.values())})
oldremote=json.loads((OLD/'configs/remote.json').read_text())
remote={'host':oldremote['host'],'port':oldremote['port'],'identity_file':oldremote['identity_file'],
        'root':'/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy','old_fem_root':'/workspace/formal_3D_flow_solver_FEM_lzy'}
write_json(ROOT/'remote/connection.local.json',remote)
print('Old FEM baseline:',len(baseline['files']),'entries; HEAD',state_before['head'])
print('Copied',len(records),'formal input files:',sum(r['size'] for r in records.values()),'bytes, all SHA256 equal.')
