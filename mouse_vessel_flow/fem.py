"""Native FEM execution interface; Python never replaces svMultiPhysics/PETSc mathematics."""
from pathlib import Path
import json,shutil,subprocess
import xml.etree.ElementTree as ET
from .paths import case_config,project_path,new_output
from .provenance import sha256

def inspect_case(name='roi_only'):
    cfg=case_config(name);xml=project_path(cfg['solver_xml']);root=ET.parse(xml)
    equation=root.find('.//Add_equation[@type="fluid"]')
    return dict(case=name,solver_xml=cfg['solver_xml'],solver_sha256=sha256(xml),
        density_kg_m3=float(equation.findtext('Density')),
        dynamic_viscosity_pa_s=float(equation.findtext('Viscosity/Value')),
        dt_s=float(root.findtext('.//Time_step_size')),
        boundary_conditions={bc.get('name'):{c.tag:c.text for c in bc} for bc in equation.findall('Add_BC')},
        linear_solver=ET.tostring(equation.find('LS'),encoding='unicode'),
        native_environment=json.loads(project_path(cfg['environment']).read_text())['native_cfd'])

def prepare_case(name,output):
    """Copy exact frozen inputs into a NEW native case; no solver is launched."""
    cfg=case_config(name);base=project_path(cfg['data_root'])/'flow'
    target=new_output(output);(target/'run').mkdir();(target/'reports').mkdir()
    reference=project_path(cfg['fem_reference'])/'frozen_reference/SV_MESH'
    shutil.copytree(reference,target/'SV_MESH')
    for relative in ('policy.json','run/solver.xml','run/PETSC_OPTIONS.txt'):
        if not (base/relative).is_file():raise FileNotFoundError(base/relative)
        shutil.copyfile(base/relative,target/relative)
    # Every XML path is resolved against the newly prepared run directory.
    xml=ET.parse(target/'run/solver.xml')
    for tag in ('Mesh_file_path','Face_file_path'):
        for node in xml.iter(tag):
            path=(target/'run'/node.text).resolve()
            if not path.is_relative_to(target) or not path.is_file():raise ValueError('Unresolved XML resource: '+str(path))
    manifest={str(p.relative_to(target)):sha256(p) for p in sorted(target.rglob('*')) if p.is_file()}
    (target/'input_hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (target/'prepared_inputs.json').write_text(json.dumps(dict(parent=name,files=manifest,solver_executed=False),indent=2)+'\n')
    return target

def run_case(directory,build_json):
    """Run the original monitored native solver with explicit external build metadata.

    Expects executable, PETSc_build and runtime_library_path from an audited
    installation. The local executable/PETSc/runtime hash gates are enforced;
    historical server binaries retain their separate original identities.
    """
    from .solver_support.native_run import run_native
    directory=Path(directory).resolve()
    proof=json.loads((directory/'prepared_inputs.json').read_text())
    for relative,expected in proof['files'].items():
        if sha256(directory/relative)!=expected:raise ValueError('Prepared inputs changed: '+relative)
    build=json.loads(Path(build_json).resolve().read_text())
    return run_native(directory,build)
