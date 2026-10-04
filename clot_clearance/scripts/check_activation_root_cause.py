"""Numerical checks directly relevant to the new component BC experiment."""
import subprocess
import sys
import numpy as np
from scripts.activation_root_cause import OUT,ROOT,read,save,arrays,gates_and_summary
from pd_clot.pfpd_candidate.activation_boundary_control import BoundaryProblem,FULL,AXIAL


def main():
    cases=read(OUT/'COMMON_STRAIN_EQUILIBRIA.json')
    assert len(cases)==8 and all(r['accepted'] for r in cases)
    zeros=[]
    for n in (5,7):
        for mode in (FULL,AXIAL):
            p=BoundaryProblem(n,0,mode);r=gates_and_summary(p,p.affine(),p.lower)
            assert r['accepted'];zeros.append(r)
    consistency=[]
    for mode in (FULL,AXIAL):
        p=BoundaryProblem(5,.003,mode)
        a=arrays(OUT/'equilibria'/f'C5/{mode}/0.00300000/STATE.npz');x=a['x'];s=a['s']
        st=p.state(x,s);work=-st.force[p.top,2].sum()*.00025
        h=1e-6;plus=BoundaryProblem(5,.003+h,mode);minus=BoundaryProblem(5,.003-h,mode)
        fd=(plus.state(plus.prescribe(x),s).energy-minus.state(minus.prescribe(x),s).energy)/(2*h)
        error=float(abs(fd/work-1));assert error<5e-6
        consistency.append(dict(mode=mode,state='accepted C5 equilibrium at epsilon .003',
            delta_epsilon=h,finite_difference_dU_depsilon_J=float(fd),reaction_times_L_J=float(work),relative_error=error,passed=True,
            convention='partial derivative at fixed free DOFs; valid reaction pairing at stationary state; no new equilibrium solve'))
    result=subprocess.run([sys.executable,'-B',str(ROOT/'tests/test_activation_boundary_control.py')],capture_output=True,text=True)
    save(OUT/'TARGETED_TESTS.json',dict(command='python -B tests/test_activation_boundary_control.py',
        exit_code=result.returncode,passed=5 if result.returncode==0 else None,output=result.stdout+result.stderr,
        scope='new BC adapter only; no full suite and no old scientific campaign'))
    assert result.returncode==0
    analyses=[read(p) for p in (OUT/'analysis').glob('C*.json')]
    max_F=max(max(r['affine']['Fbar_max_error'],r['affine']['Ftilde_max_error']) for r in analyses)
    max_CV=max(r['affine']['Psi_plus_all']['CV'] for r in analyses)
    max_kinematic=max(r['Ftilde_bond_kinematic_identity_max_error'] for r in analyses)
    assert max_F<2e-12 and max_CV<1e-8 and max_kinematic<1e-15
    record=dict(passed=True,accepted_common_strain_states=8,new_equilibrium_solves=6,reused_equilibria=2,
        max_force_residual=max(r['metrics']['force_residual'] for r in cases),minimum_J=min(r['minimum_J'] for r in cases),
        max_reaction_balance_relative=max(r['reaction_balance_relative'] for r in cases),
        boundary_exact=all(r['gates']['BC_exact'] for r in cases),phase_exact_zero=True,
        zero_load_cases=zeros,energy_reaction_consistency=consistency,affine_F_max_error=max_F,
        affine_driver_max_CV=max_CV,Ftilde_xi_equals_eta_max_error_m=max_kinematic,targeted_tests_passed=5,
        C7_stability='NOT_RUN',C9='NOT_RUN',candidate_activation='NOT_RUN_PHASE_A_INCONCLUSIVE')
    save(OUT/'NUMERICAL_GATES.json',record)
    print({k:v for k,v in record.items() if k not in ('zero_load_cases','energy_reaction_consistency')})


if __name__=='__main__':main()
