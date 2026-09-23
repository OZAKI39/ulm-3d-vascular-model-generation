"""Permanent checks include negative cases so a false PASS is detectable."""
from pathlib import Path
import sys, copy, xml.etree.ElementTree as ET
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/flow_2mmps'))
from validate import check_bc, check_mass, Q, CASE
from sv_validation.validation import triangle_flux
from sv_validation.sv13 import stop_gate

def test_only_authorized_boundary_change():
    assert check_bc(CASE/'run/solver.xml',ROOT/'frozen_reference/run/solver.xml')

@pytest.mark.parametrize('role,field,value',[
    ('INLET','Value',str(Q)),('INLET','Value',str(-Q/2)),
    ('INLET','Zero_out_perimeter','false'),('OUTLET_01','Type','Dir'),
    ('OUTLET_02','Value','1'),('WALL','Value','0.01')])
def test_unapproved_boundary_change_rejected(tmp_path,role,field,value):
    tree=ET.parse(CASE/'run/solver.xml');tree.find(f".//Add_BC[@name='{role}']/{field}").text=value
    path=tmp_path/'bad.xml';tree.write(path)
    with pytest.raises(AssertionError):check_bc(path,ROOT/'frozen_reference/run/solver.xml')

def measured(incoming=Q,out_total=Q):
    return dict(signed_outward_boundary_flows_m3_s={'INLET':-incoming},
        outlet_flows_m3_s={f'OUTLET_0{i}':out_total*f for i,f in enumerate((.1,.7,.2),1)})

def test_mass_gate_accepts_conservation():assert check_mass(measured())

@pytest.mark.parametrize('qin,qout',[(Q*.99,Q*.99),(Q,Q*.99),(-Q,Q),(float('nan'),Q)])
def test_mass_gate_rejects_wrong_target_or_conservation(qin,qout):
    with pytest.raises(AssertionError):check_mass(measured(qin,qout))

def test_exact_oriented_triangle_flux():
    # xy triangle area 3, outward +z, linearly varying normal speed 1,2,3.
    points=np.array([[0.,0.,0.],[2.,0.,0.],[0.,3.,0.]])
    velocity=np.array([[0.,0.,1.],[0.,0.,2.],[0.,0.,3.]])
    assert triangle_flux(points,np.array([[0,1,2]]),velocity)==pytest.approx(6.)
    assert triangle_flux(points,np.array([[0,2,1]]),velocity)==pytest.approx(-6.)

def test_five_consecutive_intervals_required():
    policy=dict(steady_last_intervals=5,velocity_change_limit=1e-5,flow_change_limit=1e-6,mass_limit=1e-6)
    state=dict(epsilon_Q=0,epsilon_mass=0,velocity_finite=True,pressure_finite=True,wall_noslip_pass=True)
    good=[dict(E_u=1e-7,E_Q=1e-8) for _ in range(5)]
    assert not stop_gate(good[:4],state,0,0,policy)
    assert stop_gate(good,state,0,0,policy)
    bad=copy.deepcopy(good);bad[2]['E_u']=1e-3
    assert not stop_gate(bad,state,0,0,policy)
    assert not stop_gate(good,state,1,0,policy)
