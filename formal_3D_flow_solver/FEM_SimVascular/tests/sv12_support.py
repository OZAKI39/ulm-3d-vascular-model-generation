import pytest
from sv_validation.sv12 import REPORT,load

def artifact(name):
    if not (REPORT/(name+'.json')).exists():pytest.skip('SV1.2 gate not reached: '+name)
    return load(name)

def accepted():return artifact('accepted_solution')
