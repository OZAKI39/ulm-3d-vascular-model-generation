"""Copy frozen official/vascular inputs into SV1.3J and record permitted XML differences."""
import json, shutil, sys, tarfile
import xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,write_json
R=ROOT/'reports/sv1_3j';O=ROOT/'outputs/sv1_3j'
stage=O/'flow_inputs';assert not stage.exists();stage.mkdir()
official=stage/'official_fluid_gpu_smoke';official.mkdir()
records=[]
for f in json.loads((ROOT/'reports/sv1/official_smoke_inputs.json').read_text())['files']:
    original=ROOT/f['path'];assert sha256(original)==f['sha256']
    rel=original.relative_to(ROOT/'outputs/sv1/official_fluid_smoke')
    source=ROOT/'outputs/sv1_1/official_fluid_smoke/solver.xml' if str(rel)=='solver.xml' else original
    dest=official/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dest)
    records.append({'source':str(source.relative_to(ROOT)),'destination':str(dest.relative_to(stage)),'sha256':sha256(dest)})
mesh=stage/'SV_MESH';mesh.mkdir()
mesh_sources=[ROOT/'outputs/sv1/SV_MESH/mesh-complete.mesh.vtu',*sorted((ROOT/'outputs/sv1/SV_MESH/mesh-surfaces').glob('*.vtp'))]
for source in mesh_sources:
    dest=mesh/source.relative_to(ROOT/'outputs/sv1/SV_MESH');dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dest)
    records.append({'source':str(source.relative_to(ROOT)),'destination':str(dest.relative_to(stage)),'sha256':sha256(dest)})
source=ROOT/'configs/sv1_2/sv_flow.xml';tree=ET.parse(source)
changes=[]
for tag,value in [('Continue_previous_simulation','false'),('Number_of_time_steps','20')]:
    elem=tree.find('.//'+tag);changes.append({'tag':tag,'before':elem.text,'after':value});elem.text=value
dest=stage/'vascular_proof.xml';ET.indent(tree);tree.write(dest,encoding='utf-8',xml_declaration=True)
def physics(element):
    return (element.tag,tuple(sorted(element.attrib.items())),(element.text or '').strip(),tuple(physics(c) for c in element if c.tag not in ('Continue_previous_simulation','Number_of_time_steps')))
assert physics(ET.parse(source).getroot())==physics(tree.getroot())
records.append({'source':str(source.relative_to(ROOT)),'source_sha256':sha256(source),'destination':'vascular_proof.xml','sha256':sha256(dest),'authorized_XML_changes':changes})
archive=O/'flow_inputs.tar.gz'
with tarfile.open(archive,'w:gz') as t:
    for p in stage.iterdir():t.add(p,arcname=p.name)
options=json.loads((ROOT/'configs/sv1_3/petsc_options.json').read_text())
write_json(ROOT/'configs/sv1_3j/petsc_options.json',options)
write_json(R/'flow_input_manifest.json',{'archive':str(archive.relative_to(ROOT)),'archive_sha256':sha256(archive),'files':records,
    'official_XML_source':'SV1.1 PETSc-adapted official fluid/newtonian input; unchanged',
    'vascular_XML_differences':changes,'vascular_mesh_physics_BC_dt_nonlinear_method_unchanged':True})
print('Frozen inputs packaged; vascular XML differs only in t=0 restart flag and 20-step end time')
