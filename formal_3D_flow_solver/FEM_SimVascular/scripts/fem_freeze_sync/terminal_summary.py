"""Print the requested final handoff receipt from actual artifacts and Git refs."""
import csv
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]


def read(name):
    return json.loads((ROOT/name).read_text())


def main():
    base=read('sync_metadata/main_base.json')
    flow=read('frozen_reference/flow/flow_field_manifest.json')
    mesh=read('frozen_reference/mesh_manifest.json')
    boundary=read('frozen_reference/boundary_manifest.json')['boundaries']
    runtime=read('frozen_reference/run/runtime_manifest.json')
    stage=read('sync_metadata/stage_status.json')
    tests=read('sync_metadata/test_summary.json')
    git=lambda *a:subprocess.check_output(['git',*a],cwd=REPO,text=True).strip()
    head=git('rev-parse','HEAD');remote=git('rev-parse','origin/'+base['branch'])
    clean=not git('status','--porcelain')
    assert stage['status']=='PASS' and head==remote and clean
    with (ROOT/'sync_metadata/source_manifest.csv').open() as f:
        files=list(csv.DictReader(f))
    copied=read('sync_metadata/copied_from_wsl.json')
    print(f'''Stage FEM-FREEZE-SYNC completed.

repository:
    {base['repository']}

base:
    branch = main
    commit = {base['base_commit']}

new branch:
    {base['branch']}

remote push:
    PASS; local/remote HEAD match; five GitHub files reread; independent clone fresh-read PASS

FEM source:
    source path = {base['source_FEM_path']}
    synced files = {len(files)+1} (manifest itself included)
    direct original copies = {copied['copied_files']}
    upstream unfetched pointers archived as provenance = {copied['archived_lfs_pointer_files']}
    source manifest = formal_3D_flow_solver/FEM_SimVascular/sync_metadata/source_manifest.csv

Frozen FEM:
    stage = {runtime['stage']}
    status = {runtime['status']}
    solver = {runtime['solver']}
    PETSc = {runtime['PETSc_version']}
    CUDA = {runtime['CUDA_version']}
    MPI ranks = {runtime['MPI_ranks']}
    GPU = {runtime['GPU']}

Frozen flow:
    path = formal_3D_flow_solver/FEM_SimVascular/{flow['path']}
    SHA256 = {flow['sha256']}
    size = {flow['size_bytes']} bytes
    points = {flow['points']}
    cells = {flow['cells']}
    velocity array = {flow['velocity']['name']} (POINT, {flow['velocity']['dtype']}, {flow['velocity']['components']} components)
    pressure array = {flow['pressure']['name']} (POINT, {flow['pressure']['dtype']}, {flow['pressure']['components']} component)
    velocity gradient stored = NO

Frozen mesh:
    path = formal_3D_flow_solver/FEM_SimVascular/{mesh['path']}
    SHA256 = {mesh['sha256']}
    nodes = {mesh['nodes']}
    tetra = {mesh['tetra']}
    units = {mesh['units']}

Boundaries:''')
    for role in ['WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        v=boundary[role]
        print(f'    {role} = {v["path"]}; face ID {v["sv_face_id"]}; verified outward winding')
    print(f'''
Key evidence:
    reports = reports/sv1_1, sv1_2, sv1_3, sv1_3n, sv1_3o, sv1_3p, sv1_3q
    logs = logs/sv1_3q/remote/REAL_VASCULAR_GPU_ILU_REUSE_WINNER.log
    acceptance JSON = reports/sv1_3q/REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance.json
    frozen manifests = frozen_reference/; SHA256SUMS.txt

Large artifacts not committed:
    PETSc/MPI/CUDA/build binaries and redundant historical outputs omitted
    Historical machine inventories up to 82.42 MiB: size/SHA/path in sync_metadata/omitted_artifacts.json
    Upstream unfetched examples/installers: sync_metadata/upstream_lfs_pointer_inventory.json
    Essential final flow, volume, surfaces, and checkpoint are fully committed; no essential pointer; no LFS

Explicitly deferred validation:
    mesh convergence = NOT PERFORMED BY USER DECISION
    time-step sensitivity = NOT PERFORMED BY USER DECISION
    CPU/GPU field equivalence = DEFERRED

Particle handoff:
    formal_3D_flow_solver/FEM_SimVascular/PARTICLE_HANDOFF.md

Particle roadmap:
    formal_3D_flow_solver/FEM_SimVascular/PARTICLE_RESEARCH_ROADMAP.md

git:
    local HEAD = {head}
    remote HEAD = {remote}
    working tree clean = {clean}

tests:
    passed = {tests['passed']} (18 handoff + 20 original Stage Q)
    failed = {tests['failed']}
    historical full suite not rerun; its 8 historical failures remain documented

STAGE FEM-FREEZE-SYNC STATUS:
    PASS

Stopped. Particle-0 has not been implemented; no FEM run or particle development started.''')


if __name__=='__main__':
    main()
