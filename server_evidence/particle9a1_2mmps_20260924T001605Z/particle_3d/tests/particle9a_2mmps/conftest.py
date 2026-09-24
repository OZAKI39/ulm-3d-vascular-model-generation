from pathlib import Path
import sys,ast,types
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))

@pytest.fixture(scope='session')
def legacy():
    path=Path(__file__).resolve().parents[2]/'reports/particle9a_2mmps/reference/legacy_particle_mobility.py'
    tree=ast.parse(path.read_text())
    # Execute the actual read-only reference functions without optional JIT.
    tree.body=[node for node in tree.body if not isinstance(node,ast.Try)]
    for node in ast.walk(tree):
        if isinstance(node,ast.FunctionDef):node.decorator_list=[]
    module=types.ModuleType('legacy_planar_reference')
    exec(compile(tree,str(path),'exec'),module.__dict__)
    return module
