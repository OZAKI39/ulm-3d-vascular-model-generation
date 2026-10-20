"""Only new fluid definition, load mapping and native SI port checks."""
from common import *
from native import native
from acoustic_streaming_shear import traction
import time
def main():
 start=time.perf_counter();ns,C=native(O/'checks/nonintegrating_port','inverse');s=ns.dynamics();results=[]
 for t in [.001,.0024,.0034]:
  ns.l.command(f'reset_timestep 0 time {t:.17g}');ns.l.command('run 0');s=ns.state();expected=ns.tool.force(t);uc,fc,_=C(ns.ref).evaluate(s['x'],s['d']);err=float(np.max(abs(s['force']-s['material_force']-fc-expected)))
  assert err<1e-15 and np.array_equal(s['x'],ns.ref['X'])
  results.append(dict(time_s=t,force_port_error_N=err,peak_shear_Pa=ns.tool.peak(t),no_integration=True))
 ns.close()
 rng=np.random.default_rng(1700);norm=rng.normal(size=(40,3));norm/=np.linalg.norm(norm,axis=1)[:,None];grad=rng.normal(size=(40,3,3));pressure=np.linspace(-8,8,40)
 shear,normal=traction(grad,pressure,norm);projection=float(np.max(abs(np.sum(shear*norm,axis=1))))
 pure,_=traction(np.zeros_like(grad),pressure,norm);assert projection<1e-14 and np.max(abs(pure))<1e-14
 dump(O/'KEY_CHECKS.json',dict(status='PASS',native_port=results,tangent_projection_error_Pa=projection,pure_pressure_shear_max_Pa=float(np.max(abs(pure))),fluid_check=json.loads((O/'FLUID_SHEAR_CHECK.json').read_text()),mapping_check=json.loads((O/'SMALL_SURFACE_MAPPING_AUDIT.json').read_text()),library_sha256=sha(LIB),native_integration_steps=0,wall_s=time.perf_counter()-start))
 print('KEY_CHECKS_PASS')
if __name__=='__main__':main()
