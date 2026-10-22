"""Explicit native FEM preparation; Python does not implement another CFD solver."""
from pathlib import Path
import json, shutil, xml.etree.ElementTree as ET
from .paths import ROOT, read_json, write_json, fresh_directory, sha256
from .boundary import estimate

def native_status(root=ROOT):
    config=read_json(root/'configs/cfd/native.json')
    if not config.get('executable'):
        return dict(status='BLOCKED',reason=config['reason'])
    for name in ['executable','petsc_library']:
        if sha256(config[name])!=config[name+'_sha256']:
            raise ValueError('Native binary identity changed: '+name)
    if sha256(config['wrapper'])!=config['wrapper_sha256']:
        raise ValueError('Native runtime wrapper identity changed')
    return dict(status='PASS_BINARY_IDENTITY',executable=config['executable'],
                CFD_executed_this_refactor=False,scientific_accuracy_certified=False)

def prepare_case(output, root=ROOT, mpi_ranks=8):
    """New centreline-pressure case. Historical calibration/final XMLs are read-only."""
    if type(mpi_ranks) is not int or mpi_ranks<1:raise ValueError('Positive MPI rank count required')
    output=fresh_directory(output,root)
    (output/'run').mkdir();(output/'reports').mkdir()
    template=root/'cases/balanced_pressure_final'
    xml=ET.parse(template/'run/solver.xml')
    mesh=root/'inputs/mesh/SV_MESH'
    xml.find('.//Mesh_file_path').text=str(mesh/'mesh-complete.mesh.vtu')
    for face in xml.findall('.//Add_face'):
        face.find('Face_file_path').text=str(mesh/'mesh-surfaces'/(face.get('name')+'.vtp'))
    design,edges=estimate(root)
    pressure={f'OUTLET_{i:02d}':design['outlet_pressures_Pa'][f'O{i}'] for i in (1,2,3)}
    for bc in xml.findall('.//Add_BC'):
        if bc.get('name') in pressure:
            assert bc.findtext('Type')=='Neu'
            bc.find('Value').text=format(pressure[bc.get('name')],'.17g')
    ET.indent(xml,space='  ');xml.write(output/'run/solver.xml',encoding='utf-8',xml_declaration=True)
    shutil.copy2(template/'run/PETSC_OPTIONS.txt',output/'run/PETSC_OPTIONS.txt')
    policy=read_json(template/'policy.json')
    policy.pop('solver_build_manifest',None)
    policy.update(MPI_ranks=mpi_ranks,boundary_origin='CENTRELINE_POISEUILLE_RESISTANCE_ESTIMATE',
                  mode='pressure',mesh_arrays=str(mesh/'mesh_arrays.npz'),
                  native_config=str(root/'configs/cfd/native.json'))
    write_json(output/'policy.json',policy);write_json(output/'boundary_design.json',design)
    files=[output/'run/solver.xml',output/'run/PETSC_OPTIONS.txt',output/'policy.json']+list(mesh.rglob('*.vtp'))+[mesh/'mesh-complete.mesh.vtu',mesh/'mesh_arrays.npz']
    write_json(output/'input_hashes.json',{str(p):sha256(p) for p in files})
    write_json(output/'PLAN.json',dict(executed=False,scientific_status='UNSOLVED_NEW_CASE',
        pressure_Pa=pressure,actual_3D_outlet_split='NOT_RUN',MPI_ranks=mpi_ranks,
        authoritative_mesh=str(mesh),no_mesh_symlink=True,reference_template='Historical pressure case; only new outlet pressures/paths/explicit ranks changed'))
    return output

def solve(case, root=ROOT):
    from .native_runner import run
    native_status(root)
    return run(Path(case).resolve(),read_json(root/'configs/cfd/native.json'))
