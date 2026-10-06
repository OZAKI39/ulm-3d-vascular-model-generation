"""Focused checks of the new diagnostic adapter, never a production qualifier."""
from pathlib import Path
import dataclasses
import json
import numpy as np
import pytest
from scipy import sparse
from diagnostic_core import held_alignment, increment_loads, segment_cells, strict_probe

HERE=Path(__file__).resolve().parent


def test_held_velocity_alignment_and_birth_row():
    # Unequal accepted intervals; the birth row is deliberately not a motion RHS.
    t=np.array([.65,.6503,.6508,.65081])
    held=np.array([[.2,.1,0],[.3,.1,0],[.1,-.1,.2]])
    x=np.vstack([np.zeros(3),np.cumsum(np.diff(t)[:,None]*held,axis=0)])
    saved=np.vstack([[9,9,9],held])
    dt,recovered,error=held_alignment(t,x,saved)
    np.testing.assert_array_equal(recovered,held)
    assert np.max(abs(error))<1e-18
    assert np.max(abs(np.diff(x,axis=0)-dt[:,None]*saved[:-1]))>1e-3
    with pytest.raises(ValueError):held_alignment([0,0,1],x[:3],saved[:3])


@pytest.mark.parametrize('pid',[1189,1329,642,822])
def test_original_window_alignment(pid):
    with np.load(HERE/f'original_{pid}.npz') as z:
        _,_,error=held_alignment(z['physical_time_s'],z['x_m'],z['velocity_m_s'])
        assert np.max(abs(error))<16*np.finfo(float).eps*np.max(abs(z['x_m']))


def test_increment_reconstructs_vector_with_constraint_sign_and_residual():
    from particle_3d.resistance_solver import solve_resistance
    from particle_3d.resistance_assembly import ResistanceSystem
    diag=np.array([2,3,4,2e-12,3e-12,4e-12])
    R0=np.diag(diag);R1=R0.copy();R1[0,0]*=3
    U0=np.array([.1,.2,.3,2,3,4]);U1=np.array([.3,.1,.4,3,4,5])
    q0=np.array([.1,0,0,0,0,0]);q1=-q0
    residual=np.array([1e-12,0,0,0,0,1e-23])
    b0=R0@U0-q0;b1=R1@U1-q1-residual
    loads,target,_,_=increment_loads(R0,R1,b0,b1,U0,U1,q0,q1)
    system=ResistanceSystem((0,),sparse.csr_matrix(R1),b1,diag,U1,[],[])
    reconstructed=sum(solve_resistance(dataclasses.replace(system,rhs=load)).velocity for load in loads)
    np.testing.assert_allclose(reconstructed,U1-U0,rtol=2e-14,atol=2e-15)
    np.testing.assert_allclose(loads.sum(0),target,rtol=2e-14,atol=1e-24)
    # Wrong constraint sign must fail by a material amount, even when all
    # matrices are well conditioned after physical self scaling.
    wrong=loads.sum(0)-2*(q1-q0)
    assert abs(wrong[0]-target[0])>.1


def test_three_cells_on_one_segment_and_strict_interface_membership():
    from particle_3d.field import FrozenFEMField
    points=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1],[1,1,-1]],dtype=float)
    tetra=np.array([[0,1,2,3],[1,2,3,4],[1,2,4,5]])
    field=FrozenFEMField(points,tetra,points**2,np.zeros(6))
    a=np.array([.2,.2,.2]);b=np.array([.8,.6,-.1])
    spans,crossings=segment_cells(field,a,b)
    assert [v[2] for v in spans]==[0,1,2]
    assert len(crossings)==2
    # Two analytically known planes: x+y+z=1, then x+y-z=1.
    np.testing.assert_allclose([v['fraction'] for v in crossings],[4/7,8/13],rtol=1e-14)
    for c in crossings:
        for fraction in [1e-3,1e-4,1e-5]:
            left=a+(c['fraction']-fraction)*(b-a)
            right=a+(c['fraction']+fraction)*(b-a)
            strict_probe(field,left,c['left']);strict_probe(field,right,c['right'])
            with pytest.raises(ValueError):strict_probe(field,left,c['right'])
    with pytest.raises(ValueError):strict_probe(field,np.array([-1.,-1.,-1.]))
    # Inside the production roundoff allowance still is not a valid strict
    # two-sided interior probe; no epsilon-based extrapolation is permitted.
    with pytest.raises(ValueError):strict_probe(field,np.array([-1e-15,.2,.2]))


@pytest.mark.parametrize('pid',[1189,1329,642,822])
def test_saved_increment_closure_and_complete_real_segment_enumeration(pid):
    with np.load(HERE/f'measured_{pid}.npz') as z:
        delta=np.diff(z['U'],axis=0)
    with np.load(HERE/f'increments_{pid}.npz') as z:
        terms=z['generalized_velocity_contributions'].sum(1)
    np.testing.assert_allclose(terms[:,:3],delta[:,:3],rtol=1e-11,atol=1e-15)
    np.testing.assert_allclose(terms[:,3:],delta[:,3:],rtol=1e-10,atol=1e-10)
    with np.load(HERE/f'original_{pid}.npz') as z:
        intervals=len(z['physical_time_s'])-1
        if pid==822:peak=int(np.flatnonzero(z['physical_time_s']==.6505)[0])
    result=json.loads((HERE/f'result_{pid}.json').read_text())
    assert len(result['interfaces']['all_segment_partitions'])==intervals
    if pid==822:
        assert len(result['interfaces']['all_segment_partitions'][peak]['cells'])>2
