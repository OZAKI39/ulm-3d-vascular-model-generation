"""Check renamed paths, physics, code equivalence and installed metadata."""
from pathlib import Path
import ast
import base64
import csv
import hashlib
import importlib
import importlib.metadata
import json
import os
import subprocess
import xml.etree.ElementTree as ET

C=Path('/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular')
R=Path(__file__).resolve().parent
B=json.loads((R/'before.json').read_text())
replacements=json.loads((R/'path_changes.json').read_text())['code_replacements']

def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def moved(rel):
    for old,new in sorted(B['moves'].items(),key=lambda row:len(row[0]),reverse=True):
        if rel==old or rel.startswith(old+'/'):return C/(new+rel[len(old):])
    return C/rel

changes={}
unchanged_binary=0
for rel,h in B['hashes'].items():
    p=moved(rel);assert p.is_file(),p
    if sha(p)!=h:changes[rel]=str(p.relative_to(C))
    if p.suffix in ('.npz','.npy','.vtp','.vtu','.png','.mp4','.bin','.csv') or rel.startswith('external/'):
        assert sha(p)==h,rel
        unchanged_binary+=1
for rel,state in B['git'].items():
    p=moved(rel) if rel!='.' else C
    assert subprocess.check_output(['git','-C',str(p),'rev-parse','HEAD'],text=True).strip()==state['head']
    assert sha(p/'.git/index')==state['index_sha256']

def tree(text):
    t=ast.parse(text)
    if t.body and isinstance(t.body[0],ast.Expr) and isinstance(t.body[0].value,ast.Constant) and isinstance(t.body[0].value.value,str):
        t.body.pop(0)  # Only module documentation changed in three helpers.
    return ast.dump(t,include_attributes=False)

checked_code=0
for rel in B['hashes']:
    if not rel.endswith('.py') or rel.startswith('external/'):continue
    old=(R/'before_files'/rel).read_text();expected=old
    if rel.startswith(('mesh_generate/','solver_support/')):
        for a,b in replacements.items():expected=expected.replace(a,b)
        for a,b in [('sv1_1','solver_acceptance'),('sv1_2','transient_validation'),('sv1_3','steady_convergence'),
                    ('sv_flow_petsc.xml','flow_solver_petsc.xml'),("reports/sv1'","reports/mesh_and_flow'"),
                    ('logs/mesh_and_flow_2','logs/transient_validation'),('logs/mesh_and_flow_3','logs/steady_convergence'),
                    ('SV1_GEOMETRY_IMPORT','SURFACE_MODEL_IMPORT')]:expected=expected.replace(a,b)
    elif rel in ['rotate_visualization/'+n for n in ['render_visualization.py','render_surface_fields.py','prepare_streamlines.py','validate_results.py']]:
        expected=expected.replace('SV_MESH/','solver_mesh/').replace("'SV_MESH'","'solver_mesh'")
    elif rel=='rotate_visualization/prepare_surface_data.py':
        expected=expected.replace("str(flow/'src')","str(HERE.parent/'mesh_generate/src')").replace('from sv_validation.postprocess import','from vascular_validation.postprocess import')
    assert tree(moved(rel).read_text())==tree(expected),rel
    checked_code+=1

policy=json.loads((C/'configs/mesh_policy.json').read_text())
old_policy=json.loads((R/'before_files/configs/mesh_policy.json').read_text())
assert dict(policy,formal_mesh_name='SV_MESH')==old_policy
prov=json.loads((C/'reports/mesh_and_flow/mesh_policy_provenance.json').read_text())
assert sha(C/prov['path'])==prov['sha256']
for rel in ['configs/sv_flow.xml','outputs/sv1/vascular_flow/solver.xml']:
    original=(R/'before_files'/rel).read_text()
    assert moved(rel).read_text()==original.replace('../SV_MESH/','../solver_mesh/')
    parsed=ET.parse(moved(rel))
    for node in parsed.findall('.//Mesh_file_path')+parsed.findall('.//Face_file_path'):
        assert (C/'outputs/mesh_and_flow/vascular_flow'/node.text).resolve().is_file(),node.text

for name in ['vascular_validation.geometry','vascular_validation.mesh_diagnostics','vascular_validation.postprocess',
             'vascular_validation.execution','flow_solver_support.solver_checks','flow_solver_support.transient_checks',
             'flow_solver_support.steady_monitor','flow_solver_support.parallel_checks','flow_solver_support.checkpoint_checks',
             'flow_solver_support.flow_parser']:
    m=importlib.import_module(name);assert Path(m.__file__).is_relative_to(C)
assert importlib.metadata.distribution('vascular-workflow-tools').version=='0.1.0'
site=C/'.venv/lib/python3.13/site-packages'
dist=site/'vascular_workflow_tools-0.1.0.dist-info'
for rel,digest,size in csv.reader((dist/'RECORD').open()):
    if digest:
        data=(site/rel).read_bytes()
        assert digest=='sha256='+base64.urlsafe_b64encode(hashlib.sha256(data).digest()).decode().rstrip('=')
        assert int(size)==len(data)

old_names=[]
for base,dirs,files in os.walk(C):
    dirs[:]=[d for d in dirs if d not in ('.git','.venv','external','__pycache__')]
    old_names += [str(Path(base)/name) for name in dirs+files if name.lower().startswith('sv')]
assert not old_names,old_names
for old in B['moves']:assert not (C/old).exists(),old

result=dict(all_pass=True,files_checked=len(B['hashes']),unchanged_binary_and_vendor_files=unchanged_binary,
            code_AST_equivalence_checks=checked_code,physics_parameters_unchanged=True,
            git_heads_and_indexes_unchanged=True,old_owned_sv_paths=old_names,
            imports_and_editable_metadata_pass=True,changes=changes)
(R/'renaming_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='changes'},indent=2))
