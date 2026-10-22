"""FNM cache convergence and genuine APECSS solves at every admitted frequency."""
import time,json,concurrent.futures,dataclasses
import h5py,numpy as np
from mevo_mb.acoustics.water_disc_proxy import WaterDiscs,CachedWaterBasis
from mevo_mb.bubbles.apecss import RadialCase,ApecssBubbleResponseModel,APECSS_COMMIT
from mevo_mb.bubbles.parameters import SonoVueParameterDistribution
from scripts.straight_tube_exp1.f10s.common import context as old_context
from scripts.straight_tube_exp1.f10ac_qd.common import *
from .interface import FREQUENCIES

def fields(f):
 begin=time.perf_counter();folder=AD/str(round(f));folder.mkdir(parents=True,exist_ok=True)
 elements=read(OLD/'CONFIG_SNAPSHOT.json')['geometry']['pzt'];model=WaterDiscs(elements,f)
 rng=np.random.default_rng(301810);q=rng.uniform(-40e-6,40e-6,(48,3));p,g=model.evaluate(q,96);p48,g48=model.evaluate(q,48)
 pref,gref=model.rayleigh_reference(q[:8],96,256)
 checks=dict(quadrature_pressure=float(abs(p-p48).max()),quadrature_gradient=float(abs(g-g48).max()/abs(g).max()),reference_pressure=float(abs(pref-p[:,:8]).max()),reference_gradient=float(abs(gref-g[:,:8]).max()/abs(g[:,:8]).max()))
 eps=1e-9;gd=np.stack([(model.evaluate(q+np.eye(3)[i]*eps,96)[0]-model.evaluate(q-np.eye(3)[i]*eps,96)[0])/(2*eps) for i in range(3)],-1)
 checks['analytic_gradient_finite_difference']=float(abs(g-gd).max()/abs(g).max())
 for name,step in [('coarse',3e-6),('fine',1.5e-6)]:
  axis=np.linspace(-60e-6,60e-6,round(120e-6/step)+1);points=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),-1).reshape(-1,3);pp=np.empty((4,len(points)),complex)
  for i in range(0,len(points),4096):pp[:,i:i+4096]=model.evaluate(points[i:i+4096],48)[0]
  meta=dict(model='HOMOGENEOUS_WATER_CONTROL_PROXY',frequency_Hz=f,sound_speed_m_s=1500,elements=elements,quadrature_order=48,harmonic_convention='Re[P exp(+i omega t)]',source_amplitude_meaning='rho*c*v_n Pa-equivalent; not calibrated face pressure',identity=identity())
  path=folder/(name+'.h5')
  with h5py.File(path,'x') as h:
   for s in 'xyz':h['axes/'+s+'_m']=axis
   h.create_dataset('pressure_basis',data=pp.reshape((4,)+(len(axis),)*3),compression='gzip');h.attrs['metadata_json']=json.dumps(meta)
  cp,cg=CachedWaterBasis(path).query_basis(q);checks[name+'_pressure']=float(abs(cp-p).max());checks[name+'_gradient']=float(abs(cg-g).max()/abs(g).max())
  if f==2500000.:
   old=CachedWaterBasis(OLD/'cache'/f'{name}.h5');checks[name+'_historical_pressure_exact']=bool(np.array_equal(old.pressures,pp.reshape(old.pressures.shape)))
 lines=np.linspace(-50e-6,50e-6,501);z=np.zeros((501,3));z[:,1]=lines;lp,lg=model.evaluate(z,96)
 # Pressure field XY and XZ; fields only, not a fabricated dynamical sample.
 axis=np.linspace(-60e-6,60e-6,121);grid=np.zeros((121,121,3));grid[:,:,0],grid[:,:,1]=np.meshgrid(axis,axis,indexing='ij');fp,fg=model.evaluate(grid.reshape(-1,3),96)
 np.savez_compressed(folder/'FIELD_DIAGNOSTICS.npz',points_m=q,pressure=p,gradient=g,line_y_m=lines,line_pressure=lp,line_gradient=lg,xy_axis_m=axis,xy_pressure=fp.reshape(4,121,121),center_pressure=model.evaluate([[0,0,0]],96)[0],center_gradient=model.evaluate([[0,0,0]],96)[1])
 passed=all(v<(.002 if 'gradient' in k else 1e-4) for k,v in checks.items() if 'exact' not in k)
 out=dict(frequency_Hz=f,checks=checks,PASS=passed,wall_s=time.perf_counter()-begin,identity=identity());write(folder/'FIELD_GATE.json',out);return out

def case(f,p,shell=1.,fine=False):
 par=SonoVueParameterDistribution(ROOT/'configs/stage2/sonovue_parameters.json').parameter_set(R)
 return RadialCase(R,f,float(p),par.shell_elasticity_N_m.value*shell,par.shell_viscosity_kg_s.value,.036,cycles=96 if fine else 48,samples_per_cycle=1024 if fine else 512,tolerance=5e-11 if fine else 1e-10,provenance='Unchanged historical radius-dependent SonoVue parameter manifold; elasticity multiplier '+str(shell)+'; sigma_init .036 hypothesis, not calibration')
def solve(args):
 f,p,shell,tag,fine=args;folder=AD/str(round(f))/f'shell_{shell:.1f}'/tag
 model=ApecssBubbleResponseModel(ROOT/'third_party/build/ultrasound_apecss')
 md,harm,cyc=model.solve(case(f,p,shell,fine),folder)
 return dict(pressure_Pa=float(p),case=md['case'],harmonics=harm,metrics=md['metrics'],wall_s=md['elapsed_s'],path=str(folder.relative_to(AD)))
def responses(f,shell):
 start=time.perf_counter();knots=np.array([.001,.01,.1,1,10,50,100,250,500,750,1000,1500,2000,3000,4500,6000,9000,12000.]);args=[(f,p,shell,f'node_{i:03}',False) for i,p in enumerate(knots)]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:nodes=list(pool.map(solve,args))
 values=np.array([r['harmonics'] for r in nodes]);tests=[]
 hold=(knots[:-1]+knots[1:])/2
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:refs=list(pool.map(solve,[(f,p,shell,f'hold_{i:03}',False) for i,p in enumerate(hold)]+[(f,1500,shell,'refined_1500',True)]))
 for i,r in enumerate(refs[:-1]):
  pred=(values[i]+values[i+1])/2;ref=r['harmonics'];j=np.arange(1,9);omega=2*np.pi*f
  err=float(np.linalg.norm(pred-ref)/max(np.linalg.norm(ref),1e-40));e=float(abs(np.sum(abs(omega*j*pred)**2)-np.sum(abs(omega*j*ref)**2))/max(np.sum(abs(omega*j*ref)**2),1e-80))
  tests.append(dict(pressure_Pa=r['pressure_Pa'],harmonic_relative_error=err,derivative_energy_relative_error=e,PASS=err<.01 and e<.02))
 refined=refs[-1];idx=list(knots).index(1500);convergence=float(np.linalg.norm(refined['harmonics']-values[idx])/np.linalg.norm(values[idx]))
 valid=all(n['metrics']['valid'] and n['metrics']['harmonic_energy_resolved'] for n in nodes+refs) and all(t['PASS'] for t in tests) and convergence<.01
 meta=dict(frequency_Hz=f,radius_m=R,shell_elasticity_multiplier=shell,validated=valid,apecss_commit=APECSS_COMMIT,executable_sha256=sha(ROOT/'third_party/build/ultrasound_apecss'),parameter_sha256=sha(ROOT/'configs/stage2/sonovue_parameters.json'),case=nodes[0]['case'],pressure_range_Pa=[knots[0],knots[-1]],frequency_interpolation=False,zero='exact quiescent',identity=identity())
 with h5py.File(AD/str(round(f))/f'response_shell_{shell:.1f}.h5','x') as h:
  h['pressure_Pa']=knots;h['harmonics_m3']=values;h.attrs['metadata_json']=json.dumps(meta)
 regression=None
 if f==2500000. and shell==1.:
  old=old_context()['response'];errs=[]
  for p,h in zip(knots,values):
   try:ref=old.query(np.array([complex(p)]))[0];errs.append(dict(pressure_Pa=p,error=float(np.linalg.norm(h-ref)/np.linalg.norm(ref))))
   except ValueError:pass
  regression=dict(points=errs,PASS=max(v['error'] for v in errs)<.02)
 out=dict(metadata=meta,nodes=nodes,holdouts=refs,checks=tests,refined_relative=convergence,PASS=valid,regression=regression,wall_s=time.perf_counter()-start)
 write(AD/str(round(f))/f'RESPONSE_{shell:.1f}.json',out);return out

def main():
 verify();start=time.perf_counter();results={}
 for f in [2500000.,2400000.,2600000.]:
  fld=fields(f);r=responses(f,1.);passed=fld['PASS'] and r['PASS'] and (not r['regression'] or r['regression']['PASS'])
  write(AD/str(round(f))/'PHYSICS_GATE.json',dict(PASS=passed,frequency_Hz=f,field=fld,response=r['metadata'],scope='field and radial-response admission; complete trajectories require separate numerical gates'))
  results[str(round(f))]=dict(field=fld,response_pass=r['PASS'],regression=r['regression'],PASS=passed)
  print('FREQUENCY',f,'PHYSICS',passed,flush=True)
 for f in FREQUENCIES:
  for shell in [.8,1.2]:responses(f,shell)
 write(AD/'BUILD_SUMMARY.json',dict(frequencies=results,wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
