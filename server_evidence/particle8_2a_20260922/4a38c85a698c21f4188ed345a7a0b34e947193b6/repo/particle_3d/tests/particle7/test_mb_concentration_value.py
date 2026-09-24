from particle_3d.injection_population import ConstantMBConcentrationV0,C_MB,LinearProfile,FluxClock
import pytest
def test_constant(contract):
 p=ConstantMBConcentrationV0()
 assert p.number_concentration_m3(0)==p.number_concentration_m3(100)==C_MB==8.5e12
 assert contract['MB_NOMINAL_CONCENTRATION_ML']*1e6==C_MB
 def integrate(t): return 6*t+5*t*t+4*t**3/3
 q=LinearProfile([0,1],[2,4]); c=LinearProfile([0,1],[3,5]); clock=FluxClock(q,c)
 assert clock.cumulative(.7)==pytest.approx(integrate(.7),rel=1e-14)
 assert clock.time_at(integrate(.7))==pytest.approx(.7,abs=2e-15)

def test_plateau_uses_first_threshold_crossing():
 q=LinearProfile([0,1,2,3],[2,0,0,2]); clock=FluxClock(q)
 assert clock.time_at(1)==1
 assert clock.cumulative(1)==clock.cumulative(2)==1
 assert clock.time_at(1.25)==pytest.approx(2.5)
