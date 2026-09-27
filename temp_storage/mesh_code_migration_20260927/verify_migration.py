"""Read-only production checks; write audit outputs in this receipt only."""
from pathlib import Path
import ast
import hashlib
import importlib.util
import json
import os
import subprocess
import sys

R = Path(__file__).resolve().parent
C = Path('/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular')
before = json.loads((R / 'before.json').read_text())

def sha(p):
    return hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()

def translated(relative):
    old = str(C / relative)
    if old in before['moves']:
        return Path(before['moves'][old])
    prefix = 'external/sv13q/svMultiPhysics-reuse/'
    if relative.startswith(prefix):
        relative = 'external/svMultiPhysics-reuse/' + relative[len(prefix):]
        relative = relative.replace('/sv13q_reuse.h', '/preconditioner_reuse.h')
    return C / relative

changes = []
for rel, expected in before['protected'].items():
    p = translated(rel)
    assert p.is_file(), p
    if sha(p) != expected:
        changes.append(str(p.relative_to(C)))
assert changes == ['external/svMultiPhysics-reuse/Code/Source/solver/petsc_impl.cpp'], changes
cpp = C / changes[0]
old_cpp = cpp.read_text().replace('#include "preconditioner_reuse.h"', '#include "sv13q_reuse.h"')
assert hashlib.sha256(old_cpp.encode()).hexdigest() == before['protected'][
    'external/sv13q/svMultiPhysics-reuse/Code/Source/solver/petsc_impl.cpp']

for old, state in before['git'].items():
    p = translated(old) if old != '.' else C
    assert subprocess.check_output(['git','-C',str(p),'rev-parse','HEAD'],text=True).strip() == state['head']
    assert sha(p / '.git/index') == state['index_sha256']

for old, new in before['moves'].items():
    assert not Path(old).exists(), old
    assert Path(new).exists(), new

# All moved code bytes must match after undoing just documented path/import edits.
code_checked = 0
for old, expected in before['original_code'].items():
    p = Path(before['moves'][str(C / old)])
    source = p.read_text()
    if str(p).startswith(str(C / 'mesh_generate/scripts')):
        source = source.replace('Path(__file__).resolve().parents[2]', 'Path(__file__).resolve().parents[1]')
        source = source.replace("str(ROOT / 'mesh_generate/src')", "str(ROOT / 'src')")
        source = source.replace('os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))',
                                'os.path.dirname(os.path.dirname(os.path.abspath(__file__)))')
        source = source.replace('script = Path(__file__).resolve().parent / args.script', "script = ROOT / 'scripts' / args.script")
    elif str(p).startswith(str(C / 'mesh_generate/src')):
        source = source.replace('Path(__file__).resolve().parents[3]', 'Path(__file__).resolve().parents[2]')
    elif str(p).startswith(str(C / 'solver_support/src')):
        source = source.replace('Path(__file__).resolve().parents[3]', 'Path(__file__).resolve().parents[2]')
        source = source.replace('from sv_validation.provenance import', 'from .provenance import')
        source = source.replace('from sv_validation.validation import', 'from .validation import')
        source = source.replace('from .transient_checks import', 'from .sv12 import')
        source = source.replace('from .parallel_checks import', 'from .sv13m import')
    elif p.name == 'flow_figures.py':
        source = source.replace('Path(__file__).resolve().parents[2]', 'Path(__file__).resolve().parents[1]')
        source = source.replace("str(ROOT/'mesh_generate/src')", "str(ROOT/'src')")
    assert hashlib.sha256(source.encode()).hexdigest() == expected, old
    if p.suffix == '.py':
        ast.parse(p.read_text(), filename=str(p))
    code_checked += 1

old_names = []
for base, dirs, files in os.walk(C):
    dirs[:] = [d for d in dirs if d not in ('.git', '.venv', '__pycache__', 'SimVascularDistribution')]
    for name in dirs + files:
        if name.startswith(('sv13', 'sv1_3')):
            old_names.append(str(Path(base) / name))
assert not old_names, old_names

from sv_validation import execution, visuals, mesh_diagnostics
from sv_solver_support import solver_checks, transient_checks, steady_monitor, parallel_checks, checkpoint_checks, flow_parser
for module in (execution, visuals, solver_checks, transient_checks, steady_monitor):
    assert module.ROOT == C, (module.__name__, module.ROOT)
assert Path(mesh_diagnostics.__file__).is_relative_to(C / 'mesh_generate')
assert Path(checkpoint_checks.__file__).is_relative_to(C / 'solver_support')

result = dict(all_pass=True, protected_files_checked=len(before['protected']),
              data_and_visualization_bytes_unchanged=True, code_files_checked=code_checked,
              code_changes='Only location, ROOT resolution, imports, launcher script path and C++ header include',
              old_stage_paths_remaining=old_names, old_paths_removed=len(before['moves']),
              git_heads_and_indexes_unchanged=True, module_imports_pass=True,
              mesh_code_root=str(C / 'mesh_generate'), solver_support_root=str(C / 'solver_support'))
(R / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
