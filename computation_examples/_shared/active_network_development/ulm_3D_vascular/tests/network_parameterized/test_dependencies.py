"""Static import gates plus runtime process-launch prohibition during actual fitting."""
from pathlib import Path
import ast
import builtins
import subprocess
from network_1d0d.parameterized_hydraulics import ParameterizedHydraulicSpec,solve_parameterized_operating_point
from network_1d0d.parameter_fit import Targets,Observation,fit_parameters

ROOT=Path(__file__).resolve().parents[2]
FILES=[ROOT/'network_1d0d'/name for name in ['parameterized_hydraulics.py','parameter_fit.py','parameter_config.py','parameter_freeze.py']]

def imports():
    found=[]
    for f in FILES:
        for node in ast.walk(ast.parse(f.read_text())):
            if isinstance(node,ast.Import):found.extend(n.name for n in node.names)
            elif isinstance(node,ast.ImportFrom):found.append(node.module or '')
    return found

def test_no_particle_code_dependency():
    assert not any('particle' in name.lower() or 'microbubble' in name.lower() or 'rbc' in name.lower() for name in imports())

def test_no_cfd_solver_dependency():
    forbidden=('sv','svMultiPhysics','formal_3D_flow_solver','fem_h0_case','subprocess')
    assert not any(name==bad or name.startswith(bad+'.') for name in imports() for bad in forbidden)

def test_actual_fit_has_no_process_launch_or_external_field_read(geometry_cache,monkeypatch):
    r=solve_parameterized_operating_point(geometry_cache,ParameterizedHydraulicSpec(s_O1=.92,s_O2=1.07))
    target=Targets('SYNTHETIC_TARGET',tuple(Observation('flow_fraction',p,float(f),.005) for p,f in zip(('O1','O2','O3'),r.signed_outlet_fractions)),'runtime_gate_pure_0D')
    calls=[]
    def forbidden(*a,**k):calls.append(a);raise AssertionError('No external simulation/process allowed')
    monkeypatch.setattr(subprocess,'Popen',forbidden)
    import os
    monkeypatch.setattr(os,'system',forbidden)
    oldopen=builtins.open
    def guarded_open(path,*a,**kw):
        if isinstance(path,(str,Path)) and any(s in str(path) for s in ('frozen_flow','particle','microbubble','.vtu','.vtp','H0_solution_si')):
            raise AssertionError('Prohibited external field/trajectory read')
        return oldopen(path,*a,**kw)
    monkeypatch.setattr(builtins,'open',guarded_open)
    result=fit_parameters(geometry_cache,ParameterizedHydraulicSpec(),target)
    assert result['status']=='FIT_CONVERGED' and calls==[]
    assert result['CFD_calls']==result['particle_simulation_calls']==0
