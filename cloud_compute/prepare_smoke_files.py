"""Derive a bounded cloud entry from the exact archived continuous worker, preserving physics."""
from common import *
import ast,difflib
BASE=Path('/home/lzy/projects/cloud_compute')
ROOT=Path('/home/lzy/projects/cloud_upload/20260910T212439Z/source/mirheo_starter')

def main():
    physics=ROOT/'py_scripts/single_rbc_benchmark/physics.py';text=physics.read_text();tree=ast.parse(text)
    names=['read_off','write_off','order_vertices','geometry']
    segments=[ast.get_source_segment(text,n) for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(segments)==4
    (BASE/'cloud_geometry.py').write_text('"""Unmodified function bodies extracted from archived physics.py; no legacy workflow imports."""\nfrom pathlib import Path\nimport math\nimport numpy as np\n\n'+'\n\n'.join(segments)+'\n')
    source=ROOT/'py_scripts/single_rbc_repair/continuous_worker.py';old=source.read_text();new=old
    replacements={
        'from py_scripts.single_rbc_benchmark.physics import read_off, write_off, order_vertices, geometry':'from cloud_geometry import read_off, write_off, order_vertices, geometry',
        "os.environ.get('SINGLE_RBC_AUTHORIZED')!=s['plan_sha256']":"os.environ.get('CLOUD_JOB_AUTHORIZED')!=s['cloud_plan_sha256']",
        "'AUTHORIZED_TWO_RANK_LAUNCH_REQUIRED'":"'CLOUD_AUTHORIZED_TWO_RANK_LAUNCH_REQUIRED'",
        'from .coupling import bind_bouncers':'from py_scripts.single_rbc_repair.coupling import bind_bouncers',
        'from .continuous_protocol import evolve_stage':'from py_scripts.single_rbc_repair.continuous_protocol import evolve_stage',
        "for pv,name in ((fluid,'outer'),(inner,'inner'),(rbc,'rbc')):":"# Cloud smoke observes native membrane frames only; no full-fluid frame stream.\n        for pv,name in ((rbc,'rbc'),):",
        "done=evolve_stage(u,count,dt,sample);sync()":"done=evolve_stage(u,count,dt,sample);sync()\n        if compute and int(u.getState().current_step)!=count:raise RuntimeError('NATIVE_STEP_COUNT_MISMATCH')"}
    for a,b in replacements.items():assert new.count(a)==1,a;new=new.replace(a,b)
    (BASE/'cloud_repair_worker.py').write_text(new)
    patch=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='archived/continuous_worker.py',tofile='cloud/cloud_repair_worker.py'))
    (BASE/'cloud_worker_adaptation.patch').write_text(patch)
    write(BASE/'records/worker_derivation.json',dict(original_path=str(source),original_sha256=sha(source),cloud_worker_sha256=sha(BASE/'cloud_repair_worker.py'),
        geometry_source_path=str(physics),geometry_source_sha256=sha(physics),extracted_unchanged_functions=names,
        changes=['Cloud-specific launch identity replaces old experiment gate','Import only four unchanged geometry helpers','Absolute package import for unchanged coupling and continuous protocol',
            'Save native membrane frames plus final fluid statistics/probes; omit full-fluid snapshot stream for this smoke task','Assert actual native current_step at return'],
        physics_formulas_changed=False,old_authorization_reused=False))
    print(json.dumps(dict(status='SMOKE_ENTRIES_PREPARED',worker_sha256=sha(BASE/'cloud_repair_worker.py'))))
if __name__=='__main__':main()
