import ast
import pytest
from particle_3d.lammps_bridge import validate_command
from particle_3d.lammps_state import PROPERTY_COMMAND

@pytest.mark.parametrize('command',['fix x all nve','fix x all nve/sphere','fix x all nve/asphere','fix x all langevin 1 1 1 1','fix x all brownian 1 1','fix x all viscous 1','fix x all addforce 1 0 0','run 1','run 100','run 0 every 1 "run 1"','pair_style lj/cut 1','velocity all set 1 0 0','mass * 2','include dynamics.in','run 0; run 1','run ${n}','atom_style sphere'])
def test_forbidden_command_rejected(command):
    with pytest.raises(ValueError):validate_command(command)

def test_production_command_literals_and_api_audit(repo):
    path=repo/'particle_3d/src/particle_3d/lammps_bridge.py'
    tree=ast.parse(path.read_text())
    for node in ast.walk(tree):
        if isinstance(node,ast.Constant) and isinstance(node.value,str):
            s=node.value
            if s=='pair_style zero ':validate_command(s+'1e-6') # literal prefix of the audited cutoff f-string
            elif s.startswith(('fix ','run ','pair_style ','mass ','atom_style ')):validate_command(s)
    validate_command(PROPERTY_COMMAND)
    assert path.read_text().count('._lmp.command(')==1
