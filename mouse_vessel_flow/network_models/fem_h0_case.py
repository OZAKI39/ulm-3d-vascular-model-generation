"""Create a distinct FEM input case, permitting only outlet pressure changes."""
from pathlib import Path
import json
import shutil
import xml.etree.ElementTree as ET
import numpy as np
from .audit import sha256, write_json

OUTLETS = ('OUTLET_01', 'OUTLET_02', 'OUTLET_03')


def serialize_fixed_pressure_xml(original, values):
    if len(values) != 3 or not np.isfinite(values).all():
        raise ValueError('Exactly three finite cap pressures are required')
    root = ET.fromstring(original)
    for name, value in zip(OUTLETS, values):
        bc = root.find(f".//Add_BC[@name='{name}']")
        if bc is None or bc.findtext('Type') != 'Neu':
            raise ValueError('Expected existing Neumann outlet '+name)
        bc.find('Value').text = format(float(value), '.17g')
    return ET.tostring(root, encoding='utf-8', xml_declaration=True)


def audit_xml_change(original, updated, values):
    old, new = ET.fromstring(original), ET.fromstring(updated)
    fields = []
    for name, value in zip(OUTLETS, values):
        path = f".//Add_BC[@name='{name}']/Value"
        a, b = old.find(path), new.find(path)
        if float(b.text) != float(value):
            raise ValueError('Serialized cap pressure mismatch')
        fields.append(dict(field=path, old=float(a.text), new=float(b.text),
                           numerically_changed=float(a.text) != float(b.text)))
        b.text = a.text
    if ET.tostring(old) != ET.tostring(new):
        raise ValueError('Unexpected non-outlet-pressure configuration change')
    return fields


def create_case(original, target, bc_path):
    original, target = Path(original), Path(target)
    if target.exists():
        raise FileExistsError('Never overwrite an existing case: '+str(target))
    bc = json.loads(Path(bc_path).read_text())
    if bc.get('status') != 'PASS':
        raise ValueError('Pressure transfer must pass before creating a FEM case')
    values = [next(r['pressure_cap_shifted_Pa'] for r in bc['ports'] if r['port']==name)
              for name in ('O1','O2','O3')]
    old_xml = (original/'run/solver.xml').read_bytes()
    new_xml = serialize_fixed_pressure_xml(old_xml, values)
    fields = audit_xml_change(old_xml, new_xml, values)
    target.mkdir(); (target/'run').mkdir(); (target/'reports').mkdir()
    shutil.copytree(original/'SV_MESH', target/'SV_MESH')
    for name in ('policy.json', 'run/PETSC_OPTIONS.txt'):
        shutil.copyfile(original/name, target/name)
    (target/'run/solver.xml').write_bytes(new_xml)
    names = sorted(str(p.relative_to(target)) for p in (target/'SV_MESH').rglob('*') if p.is_file())
    names += ['policy.json', 'run/PETSC_OPTIONS.txt']
    unchanged = {name: sha256(original/name) for name in names}
    assert all(sha256(target/name)==digest for name,digest in unchanged.items())
    hashes = {name:sha256(target/name) for name in names+['run/solver.xml']}
    write_json(target/'input_hashes.json', hashes)
    return dict(status='PASS', original_case=str(original), new_case=str(target),
                fields=fields, only_outlet_pressure_fields_allowed=True,
                unchanged_input_hashes=unchanged, original_solver_sha256=sha256(original/'run/solver.xml'),
                new_solver_sha256=sha256(target/'run/solver.xml'), source_BC_SHA=sha256(bc_path),
                formulation='unchanged P1/P1 + VMS',
                inlet_Q_preserved_m3s=1.551359160885543e-14,
                network_measured_target_Q_m3s=1.551359160440232e-14,
                inlet_target_note='Preserve production XML exactly; measured network target differs by 2.87e-10 relative')
