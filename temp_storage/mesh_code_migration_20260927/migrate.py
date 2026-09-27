"""One-time path migration; preserve data and record original bytes."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess

C = Path('/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular')
RECEIPT = Path(__file__).resolve().parent

def sha(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()

moves = {}
for name in ('prepare_surface.py', 'sv_import_model.py', 'sv_generate_mesh.py',
             'audit_mesh.py', 'freeze_mesh_policy.py', 'mesh_figures.py',
             'run_sv_python.py', 'sv_api_probe.py', 'sv_inspect_api.py'):
    moves['scripts/' + name] = C / 'mesh_generate/scripts' / name
solver_modules = {
    'sv11.py': 'solver_checks.py',
    'sv12.py': 'transient_checks.py',
    'sv13.py': 'steady_monitor.py',
    'sv13m.py': 'parallel_checks.py',
    'sv13n.py': 'checkpoint_checks.py',
}
for path in (C / 'src/sv_validation').glob('*.py'):
    dest = (C / 'solver_support/src/sv_solver_support' / solver_modules[path.name]
            if path.name in solver_modules
            else C / 'mesh_generate/src/sv_validation' / path.name)
    moves[str(path.relative_to(C))] = dest
moves.update({
    'scripts/flow_figures.py': C / 'solver_support/scripts/flow_figures.py',
    'scripts/sv13n/use_gpu13_env.sh': C / 'solver_support/use_gpu_env.sh',
    'scripts/sv13q/flow_parser.py': C / 'solver_support/src/sv_solver_support/flow_parser.py',
    'scripts/sv13q/sv13q_reuse.h': C / 'solver_support/include/preconditioner_reuse.h',
    'scripts/sv13q/patch_source.py': RECEIPT / 'reference_recipes/create_preconditioner_patch.py',
    'scripts/sv13q/setup_remote.py': RECEIPT / 'reference_recipes/build_remote_solver.py',
    'configs/sv1_3q/policy.json': C / 'solver_support/configs/preconditioner_reuse_policy.json',
    'external/sv13q/svMultiPhysics-reuse': C / 'external/svMultiPhysics-reuse',
})
assert not (RECEIPT / 'before.json').exists(), 'Migration already started'
assert not any((C / 'mesh_generate').iterdir()), 'Target mesh directory is not empty'
for old, new in moves.items():
    assert (C / old).exists(), old
    assert not new.exists(), new

# Data, artifacts and dependency source hashes allow a byte-level post-move check.
protected = {}
for root in ('inputs', 'outputs', 'reports', 'configs', 'rotate_visualization',
             'external/sv13q/svMultiPhysics-reuse'):
    for base, dirs, files in os.walk(C / root):
        dirs[:] = [d for d in dirs if d not in ('.git', '__pycache__', '.pytest_cache')]
        for name in files:
            p = Path(base) / name
            if p.is_file():
                protected[str(p.relative_to(C))] = sha(p)
original_code = {}
for old in moves:
    p = C / old
    if p.is_file() and p.suffix in ('.py', '.sh', '.h'):
        original_code[old] = sha(p)
        backup = RECEIPT / 'before_code' / old
        backup.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, backup)
metadata = [C / 'pyproject.toml', C / 'README.md',
            C / '.venv/lib/python3.13/site-packages/__editable__.sv_feasibility_validation-0.1.0.pth']
for p in metadata:
    backup = RECEIPT / 'before_metadata' / p.relative_to(C)
    backup.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(p, backup)
source = C / 'external/sv13q/svMultiPhysics-reuse'
git_state = {str(p.relative_to(C)): {
    'head': subprocess.check_output(['git', '-C', str(p), 'rev-parse', 'HEAD'], text=True).strip(),
    'index_sha256': sha(p / '.git/index'),
} for p in (C, source)}
record = dict(root=str(C), moves={str(C/k): str(v) for k,v in moves.items()},
              original_code=original_code, protected=protected, git=git_state)
(RECEIPT / 'before.json').write_text(json.dumps(record, indent=2) + '\n')

for old, new in moves.items():
    new.parent.mkdir(parents=True, exist_ok=True)
    (C / old).rename(new)

# All Python helpers continue resolving data paths against the original case root.
for p in (C / 'mesh_generate/scripts').glob('*.py'):
    text = p.read_text()
    text = text.replace('Path(__file__).resolve().parents[1]', 'Path(__file__).resolve().parents[2]')
    text = text.replace("str(ROOT / 'src')", "str(ROOT / 'mesh_generate/src')")
    text = text.replace('os.path.dirname(os.path.dirname(os.path.abspath(__file__)))',
                        'os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))')
    if p.name == 'run_sv_python.py':
        text = text.replace("script = ROOT / 'scripts' / args.script", "script = Path(__file__).resolve().parent / args.script")
    p.write_text(text)
for p in (C / 'mesh_generate/src/sv_validation').glob('*.py'):
    text = p.read_text().replace('Path(__file__).resolve().parents[2]', 'Path(__file__).resolve().parents[3]')
    p.write_text(text)
for p in (C / 'solver_support/src/sv_solver_support').glob('*.py'):
    text = p.read_text().replace('Path(__file__).resolve().parents[2]', 'Path(__file__).resolve().parents[3]')
    text = text.replace('from .provenance import', 'from sv_validation.provenance import')
    text = text.replace('from .validation import', 'from sv_validation.validation import')
    text = text.replace('from .sv12 import', 'from .transient_checks import')
    text = text.replace('from .sv13m import', 'from .parallel_checks import')
    p.write_text(text)
p = C / 'solver_support/scripts/flow_figures.py'
text = p.read_text().replace('Path(__file__).resolve().parents[1]', 'Path(__file__).resolve().parents[2]')
p.write_text(text.replace("str(ROOT/'src')", "str(ROOT/'mesh_generate/src')"))
(C / 'solver_support/src/sv_solver_support/__init__.py').write_text('"""Solver log, checkpoint and numerical acceptance helpers."""\n')

# Functional header filename only; keep C++ implementation and runtime protocol intact.
source = C / 'external/svMultiPhysics-reuse'
header = source / 'Code/Source/solver/sv13q_reuse.h'
header.rename(header.with_name('preconditioner_reuse.h'))
cpp = source / 'Code/Source/solver/petsc_impl.cpp'
text = cpp.read_text()
assert text.count('#include "sv13q_reuse.h"') == 1
cpp.write_text(text.replace('#include "sv13q_reuse.h"', '#include "preconditioner_reuse.h"'))

p = C / 'pyproject.toml'
p.write_text(p.read_text().replace('["src"]', '["mesh_generate/src", "solver_support/src"]'))
pth = C / '.venv/lib/python3.13/site-packages/__editable__.sv_feasibility_validation-0.1.0.pth'
pth.write_text(str(C / 'mesh_generate/src') + '\n' + str(C / 'solver_support/src') + '\n')

# Remove empty old stage containers, never data or unknown files.
for root in ('scripts', 'src', 'configs/sv1_3q', 'external/sv13q'):
    for base, dirs, files in os.walk(C / root, topdown=False):
        p = Path(base)
        if p.name == '__pycache__':
            shutil.rmtree(p)
        elif not any(p.iterdir()):
            p.rmdir()
print(json.dumps(dict(moved_entries=len(moves), protected_files=len(protected), receipt=str(RECEIPT)), indent=2))
