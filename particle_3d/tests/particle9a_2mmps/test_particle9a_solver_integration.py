import numpy as np
from particle_3d.particle9a_motion import Particle9AStepper
from particle_3d.particle81_simulation import SavedTrajectoryStepper
from particle_3d.particle65_cases import synthetic,records_for,MU
from particle_3d.lammps_neighbors import ValidationNeighborPolicy

def test_particle9a_solver_integration():
    shapes,_,wall,_=synthetic('wall')
    def normal(i,s,t):return np.array([0.,0.,-1e-4]),np.array([0.,0.,.3])
    args=(records_for(shapes,normal),ValidationNeighborPolicy(20e-6,.5e-6,'P9_TEST'),normal,MU)
    old=SavedTrajectoryStepper(*args,wall=wall)
    new=Particle9AStepper(*args,wall=wall,gradient_provider=lambda x:np.zeros((3,3)))
    for t in [.001,.004,.016,.032]:
        old.step_to(t);new.step_to(t)
        np.testing.assert_allclose(old.read()[0].position,new.read()[0].position,rtol=0,atol=1e-20)
        np.testing.assert_allclose(old.read()[0].velocity,new.read()[0].velocity,rtol=4e-15,atol=1e-20)
    assert new.planar_statistics['assembly_calls']>0
    assert new.minimum_gap>=2e-9-wall.roundoff_m
    assert new.rejected_count==old.rejected_count


def test_coupled_shear_reaches_accepted_velocity_and_quaternion(legacy):
    from particle_3d.planar_wall_hydrodynamics import legacy_planar_coefficients_si
    shapes,_,wall,_=synthetic('wall')
    gamma=500.;a=1e-6;xi=.01
    def shear(i,s,t):return np.array([gamma*(a+xi*a),0.,0.]),np.array([0.,gamma/2,0.])
    gradient=np.zeros((3,3));gradient[0,2]=gamma
    query=ValidationNeighborPolicy(20e-6,.5e-6,'P9_AFFINE_SOLVER_TEST')
    stepper=Particle9AStepper(records_for(shapes,shear),query,shear,MU,wall=wall,gradient_provider=lambda x:gradient)
    stepper.step_to(1e-6);p=stepper.read()[0]
    old=legacy.background_hydrodynamic_velocity_scalar(505.,0.,0.,gamma,0.,0.,0.,1.,1.,xi,MU)
    coeff=legacy_planar_coefficients_si(xi);fs,ty,_=legacy.wall_shear_load_scalar(MU,1.,xi,-gamma)
    base=1/(6*np.pi*MU*a);e=coeff['projection_entry_error']
    assert abs(p.velocity[0]-old[0]*1e-6)<=base*e*abs(ty*1e-18/a)+1e-15
    assert abs(p.omega[1]-old[2])<=base*e*abs(fs*1e-12)/a+1e-9
    assert p.omega[1]<0 and p.velocity[0]>0
    np.testing.assert_allclose(p.position,shapes[17].center_m+1e-6*p.velocity,rtol=0,atol=1e-20)
    assert not np.array_equal(p.q,stepper.samples[0][10:14])
