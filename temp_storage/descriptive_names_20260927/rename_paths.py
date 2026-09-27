"""Rename owned project paths, preserving numerical and third-party interfaces."""
from pathlib import Path
import base64
import csv
import hashlib
import json
import os
import shutil
import subprocess

C = Path('/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular')
R = Path(__file__).resolve().parent
V = C / 'rotate_visualization'

def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def dump(p, d):
    p.write_text(json.dumps(d, indent=2, ensure_ascii=False) + '\n')

moves = {
    'mesh_generate/src/sv_validation': 'mesh_generate/src/vascular_validation',
    'solver_support/src/sv_solver_support': 'solver_support/src/flow_solver_support',
    'mesh_generate/scripts/sv_api_probe.py': 'mesh_generate/scripts/probe_meshing_api.py',
    'mesh_generate/scripts/sv_inspect_api.py': 'mesh_generate/scripts/inspect_meshing_api.py',
    'mesh_generate/scripts/sv_generate_mesh.py': 'mesh_generate/scripts/generate_tetra_mesh.py',
    'mesh_generate/scripts/sv_import_model.py': 'mesh_generate/scripts/import_surface_model.py',
    'mesh_generate/scripts/run_sv_python.py': 'mesh_generate/scripts/run_meshing_python.py',
    'configs/sv_flow.xml': 'configs/flow_solver.xml',
    'configs/sv_reference.yaml': 'configs/reference_physics.yaml',
    'outputs/sv1': 'outputs/mesh_and_flow',
    'reports/sv1': 'reports/mesh_and_flow',
    'rotate_visualization/input_data/SV_MESH': 'rotate_visualization/input_data/solver_mesh',
    'external/svMultiPhysics-reuse': 'external/flow_solver_source',
}
site = C / '.venv/lib/python3.13/site-packages'
old_dist = 'sv_feasibility_validation-0.1.0.dist-info'
new_dist = 'vascular_workflow_tools-0.1.0.dist-info'
old_pth = '__editable__.sv_feasibility_validation-0.1.0.pth'
new_pth = '__editable__.vascular_workflow_tools-0.1.0.pth'
moves[str((site/old_dist).relative_to(C))] = str((site/new_dist).relative_to(C))
moves[str((site/old_pth).relative_to(C))] = str((site/new_pth).relative_to(C))
if (C / 'logs/sv1').exists():
    moves['logs/sv1'] = 'logs/mesh_and_flow'
assert not (R / 'before.json').exists(), 'Already started'
for old,new in moves.items():
    assert (C/old).exists(), old
    assert not (C/new).exists(), new

hashes = {}
for base,dirs,files in os.walk(C):
    dirs[:] = [d for d in dirs if d not in ('.git','.venv','__pycache__','SimVascularDistribution')]
    for name in files:
        p=Path(base)/name
        if p.is_file():
            rel=p.relative_to(C)
            hashes[str(rel)] = sha(p)
            # Only small owned source/metadata; no duplicate datasets/vendor tree.
            if rel.parts[0]!='external' and p.suffix in ('.py','.sh','.md','.json','.xml','.yaml','.toml','.txt'):
                target=R/'before_files'/rel
                target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(p,target)
for p in [site/old_pth, *(site/old_dist).iterdir()]:
    target=R/'before_files'/p.relative_to(C)
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(p,target)
git_state={str(p.relative_to(C)):dict(
    head=subprocess.check_output(['git','-C',str(p),'rev-parse','HEAD'],text=True).strip(),
    index_sha256=sha(p/'.git/index')) for p in [C,C/'external/svMultiPhysics-reuse']}
all_moves=dict(moves)
all_moves['outputs/sv1/SV_MESH']='outputs/mesh_and_flow/solver_mesh'
dump(R/'before.json',dict(root=str(C),moves=all_moves,hashes=hashes,git=git_state))

for old,new in moves.items():
    (C/old).rename(C/new)
(C/'outputs/mesh_and_flow/SV_MESH').rename(C/'outputs/mesh_and_flow/solver_mesh')

replacements = {
    'sv_validation': 'vascular_validation',
    'sv_solver_support': 'flow_solver_support',
    'sv_api_probe.py': 'probe_meshing_api.py',
    'sv_inspect_api.py': 'inspect_meshing_api.py',
    'sv_generate_mesh.py': 'generate_tetra_mesh.py',
    'sv_import_model.py': 'import_surface_model.py',
    'run_sv_python.py': 'run_meshing_python.py',
    'sv_environment': 'meshing_environment',
    'SV1_API_PROBE_JSON': 'MESH_API_PROBE_JSON',
    'SV1_MESH_GENERATION': 'TETRA_MESH_GENERATION',
    'SV1_MESH_FACTOR': 'MESH_EDGE_FACTOR',
    'outputs/sv1/': 'outputs/mesh_and_flow/',
    'reports/sv1/': 'reports/mesh_and_flow/',
    'logs/sv1': 'logs/mesh_and_flow',
    'SV_MESH': 'solver_mesh',
}
edited=[]
for root in ['mesh_generate','solver_support']:
    for p in (C/root).rglob('*.py'):
        text=p.read_text();new=text
        for old,value in replacements.items(): new=new.replace(old,value)
        # Bare report paths and dormant helper defaults also use functional names.
        for old,value in [('sv1_1','solver_acceptance'),('sv1_2','transient_validation'),('sv1_3','steady_convergence'),
                          ('sv_flow_petsc.xml','flow_solver_petsc.xml'),("reports/sv1'","reports/mesh_and_flow'")]:
            new=new.replace(old,value)
        if new!=text:p.write_text(new);edited.append(str(p.relative_to(C)))

for p in [C/'configs/flow_solver.xml',C/'outputs/mesh_and_flow/vascular_flow/solver.xml']:
    p.write_text(p.read_text().replace('../SV_MESH/','../solver_mesh/'))

# Rename only the display/path label in the frozen mesh policy, retain every physical value.
p=C/'configs/mesh_policy.json';original_policy=json.loads(p.read_text())
policy=dict(original_policy,formal_mesh_name='solver_mesh');dump(p,policy)
prov=C/'reports/mesh_and_flow/mesh_policy_provenance.json'
d=json.loads(prov.read_text());d['sha256']=sha(p);dump(prov,d)
prov=C/'reports/mesh_and_flow/flow_config_provenance.json'
d=json.loads(prov.read_text());d['xml_sha256_before_path_rename']=d['xml_sha256']
d['xml_sha256']=sha(C/'configs/flow_solver.xml');d['path_rename_only']=True;dump(prov,d)

# Renderer edits are strictly literal input path substitutions, preserving the style.
renderers=['render_visualization.py','render_surface_fields.py']
for name in renderers+['prepare_streamlines.py','validate_results.py']:
    p=V/name;p.write_text(p.read_text().replace('SV_MESH/','solver_mesh/').replace("'SV_MESH'","'solver_mesh'"))
p=V/'prepare_surface_data.py'
text=p.read_text().replace("str(flow/'src')", "str(HERE.parent/'mesh_generate/src')")
text=text.replace('from sv_validation.postprocess import', 'from vascular_validation.postprocess import')
# The source case belongs to another frozen worktree and retains its real directory name.
p.write_text(text)

def rename_local_mesh_paths(obj):
    if isinstance(obj,dict):
        return {k.replace('SV_MESH/','solver_mesh/'):rename_local_mesh_paths(v) for k,v in obj.items()}
    if isinstance(obj,list):return [rename_local_mesh_paths(v) for v in obj]
    if isinstance(obj,str) and (obj.startswith('SV_MESH/') or obj.startswith('input_data/SV_MESH/')):
        return obj.replace('SV_MESH/','solver_mesh/')
    return obj

for rel in ['input_data/streamlines/SOURCE_LOCK.json','INPUT_MANIFEST.json',
            'results/MEDIA_VALIDATION.json','results/surface_fields/MEDIA_VALIDATION.json']:
    p=V/rel;dump(p,rename_local_mesh_paths(json.loads(p.read_text())))

render_provenance={}
for name in renderers:
    old=(R/'before_files/rotate_visualization'/name).read_text()
    current=(V/name).read_text()
    assert current==old.replace('SV_MESH/','solver_mesh/').replace("'SV_MESH'","'solver_mesh'")
    render_provenance[name]=dict(previous_sha256=hashes['rotate_visualization/'+name],
                                current_sha256=sha(V/name),only_mesh_directory_literal_changed=True)
for rel,script in [('results/MEDIA_VALIDATION.json','render_visualization.py'),
                   ('results/surface_fields/MEDIA_VALIDATION.json','render_surface_fields.py')]:
    p=V/rel;d=json.loads(p.read_text())
    d['render_script_sha256_before_path_rename']=d['script_sha256']
    d['script_sha256']=sha(V/script)
    if 'shared_renderer_sha256' in d:
        d['shared_renderer_sha256_before_path_rename']=d['shared_renderer_sha256']
        d['shared_renderer_sha256']=sha(V/'render_visualization.py')
    d['path_rename_only']=True;d['media_regenerated_during_path_rename']=False
    dump(p,d)
dump(V/'PATH_RENAME_PROVENANCE.json',dict(renderers=render_provenance,
    input_path_map={'input_data/SV_MESH':'input_data/solver_mesh'},
    media_regenerated=False, numerical_and_presentation_code_unchanged=True,
    historical_audit_snapshots_preserved=True, receipt=str(R/'before.json')))

# The input inventory describes the renamed local bundle, with unchanged dataset hashes.
p=V/'INPUT_MANIFEST.json';d=json.loads(p.read_text())
for rel,entry in d['files'].items():
    f=V/rel;assert f.is_file(),f
    entry['sha256']=sha(f);entry['bytes']=f.stat().st_size
dump(p,d)

# Install metadata is kept coherent without downloading/reinstalling dependencies.
p=C/'pyproject.toml';p.write_text(p.read_text().replace('sv-feasibility-validation','vascular-workflow-tools'))
dist=site/new_dist
p=dist/'METADATA';p.write_text(p.read_text().replace('Name: sv-feasibility-validation','Name: vascular-workflow-tools'))
(dist/'top_level.txt').write_text('flow_solver_support\nvascular_validation\n')
rows=list(csv.reader((dist/'RECORD').open()))
for row in rows:
    row[0]=row[0].replace(old_dist,new_dist).replace(old_pth,new_pth)
    if row[1]:
        data=(site/row[0]).read_bytes()
        row[1]='sha256='+base64.urlsafe_b64encode(hashlib.sha256(data).digest()).decode().rstrip('=')
        row[2]=str(len(data))
with (dist/'RECORD').open('w',newline='') as stream:csv.writer(stream).writerows(rows)
dump(R/'path_changes.json',dict(moves=all_moves,code_replacements=replacements,edited=edited))
print(json.dumps(dict(renamed_paths=len(all_moves),edited_code_files=len(edited),receipt=str(R)),indent=2))
