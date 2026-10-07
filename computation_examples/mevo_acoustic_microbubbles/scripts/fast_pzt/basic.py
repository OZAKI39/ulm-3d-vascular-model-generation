"""Necessary small GPU/formula checks, then actual A/B aggregation ON/OFF."""
import argparse
from dataclasses import replace
import numpy as np
import warp as wp
from mevo_mb.fast_sim.models import Physics,ResponseSection,Basis,random_cloud,pair_forces_numpy
from mevo_mb.fast_sim.engine import FastSimulator,complex_vectors
from mevo_mb.fast_sim import kernels as k
from mevo_mb.bubbles.forces import primary_bjerknes
from .common import BASE,inputs,write,source_identity
from .experiment import run_case


def main(label):
    out=BASE/label; out.mkdir(parents=True,exist_ok=False); checks=[]
    p=Physics(); response=ResponseSection(inputs()[0],p); basis=Basis(inputs()[1],p.frequency_Hz)
    x=random_cloud(4,p.radius_m,12e-6,17,center=[0,80e-6,40e-6])
    sim=FastSimulator(x,p,response,basis)
    event={'time_s':0,'amplitude_Pa':[500,500,500,500],'phase_rad':[0,.2,.4,.6]}
    sim.control(event); sim.fields(); wp.synchronize()
    field=sim.last_field; pressure=field.interpolate_complex_pressure(x); gradient=field.interpolate_complex_gradient(x)
    volume=response.query(pressure)
    np.testing.assert_allclose(sim.pressure.numpy().view(np.complex128).ravel(),pressure,rtol=3e-12,atol=1e-10)
    np.testing.assert_allclose(sim.harmonics.numpy().view(np.complex128).reshape(4,8),volume,rtol=3e-12,atol=1e-30)
    np.testing.assert_allclose(sim.force.numpy(),primary_bjerknes(volume[:,0],gradient),rtol=3e-12,atol=1e-25)
    checks.append({'check':'real_basis_interpolation_phase_response_primary_numpy','status':'PASS','compile_s':sim.compile_s})
    primary=sim.force.numpy(); sim.pairs(); wp.synchronize()
    np.testing.assert_allclose(sim.force.numpy()-primary,pair_forces_numpy(x,volume,p),rtol=3e-12,atol=1e-25)
    checks.append({'check':'multiharmonic_full_pair_numpy','status':'PASS'})
    sim.reset(); off=dict(event,enabled=[False]*4)
    run=sim.run(20e-6,5e-6,[off],sample_interval_s=5e-6)
    assert run['summary']['status']=='COMPLETE' and np.array_equal(run['positions_m'][-1,0],x)
    checks.append({'check':'all_PZT_OFF_no_spurious_motion','status':'PASS'})
    xx=np.array([[-6e-6,0,0],[6e-6,0,0]])
    pair=FastSimulator(xx,p,response)
    for sign in [1,-1]:
        v=np.zeros((2,8),complex); v[:,0]=[1e-19,sign*1e-19]
        pair.harmonics.assign(complex_vectors(v)); pair.force.zero_(); pair.pairs(); wp.synchronize()
        f=pair.force.numpy(); np.testing.assert_allclose(f,pair_forces_numpy(xx,v,p),rtol=2e-12,atol=1e-26)
        assert np.sign(f[0,0])==sign; np.testing.assert_allclose(f.sum(axis=0),0,atol=1e-26)
    checks.append({'check':'attraction_repulsion_sign_and_reciprocity','status':'PASS',
                   'scope':'Manufactured in-phase/counterphase fundamental responses; identical bubbles in a truly uniform field attract'})
    pair.reset(); pair.control({'time_s':0,'pressure_Pa':1000}); pair.fields(); pair.pairs()
    pair.x.assign(np.array([[-1.005*p.radius_m,0,0],[1.005*p.radius_m,0,0]]))
    free=np.array([[1e-3,0,0],[-1e-3,0,0]]); pair.free.assign(free)
    pair.launch(k.local_mobility,2,[pair.x,pair.free,pair.free,pair.vel1,2,p.radius_m,.1,.001]); wp.synchronize()
    c=.5*((.5/.01-.45*np.log(.01)-3/56*.01*np.log(.01))-(.5/.1-.45*np.log(.1)-3/56*.1*np.log(.1)))
    expected=.5*(free+(free+c*free[::-1])/(1+c))
    np.testing.assert_allclose(pair.vel1.numpy(),expected,rtol=1e-11,atol=1e-17)
    assert np.linalg.norm(pair.vel1.numpy())<np.linalg.norm(free)
    checks.append({'check':'local_normal_mobility_damps_approach','status':'PASS'})
    overlap=np.array([[-.99*p.radius_m,0,0],[.99*p.radius_m,0,0]])
    pair.tmp1.assign(overlap)
    pair.launch(k.contact,2,[pair.tmp1,pair.tmp2,2,p.radius_m,.5,pair.diag,pair.error]); wp.synchronize()
    separated=pair.tmp2.numpy(); np.testing.assert_allclose(separated.mean(axis=0),overlap.mean(axis=0),atol=1e-22)
    assert np.sum(separated**2)>np.sum(overlap**2)
    checks.append({'check':'contact_is_separating_and_preserves_centroid','status':'PASS',
        'analytic_scope':'For equal Stokes mobility, each simultaneous pair correction is antisymmetric and sum_i x_i dot delta_i >= 0; Rg squared cannot decrease'})
    pair.reset(); run=pair.run(20e-6,3e-6,[{'time_s':0,'pressure_Pa':0},{'time_s':7e-6,'pressure_Pa':1000}],sample_interval_s=3e-6)
    assert run['summary']['status']=='COMPLETE' and 7e-6 in run['time_s']
    np.testing.assert_array_equal(run['positions_m'][run['time_s']<=7e-6,0],np.broadcast_to(xx,(np.count_nonzero(run['time_s']<=7e-6),2,3)))
    assert not np.array_equal(run['positions_m'][-1,0],xx)
    checks.append({'check':'noncommensurate_PZT_switch_alignment','status':'PASS','macro_steps':run['summary']['macro_step_count']})
    for query in [0.05,40001]:
        try: response.query(np.array([query])); raise AssertionError('Extrapolation accepted')
        except ValueError: pass
    try: Basis(inputs()[1],2.4e6); raise AssertionError('Wrong basis frequency accepted')
    except ValueError: pass
    checks.append({'check':'response_domain_and_basis_frequency_rejection','status':'PASS'})
    # Continuity at the exact transition, same signed complex correlations.
    v=response.query(np.array([1000,1000])); transition=4*p.radius_m
    forces=[pair_forces_numpy(np.array([[0.,0,0],[transition*(1+eps),0,0]]),v,p)[0,0] for eps in [-1e-7,0,1e-7]]
    assert max(abs(np.array(forces)/forces[1]-1))<3e-7
    checks.append({'check':'smooth_near_field_transition','status':'PASS','relative_span':float(np.ptp(forces)/forces[1])})
    write(out/'GPU_CHECKS.json',dict(status='PASS',checks=checks,source=source_identity()))
    scene=dict(count=32,seed=20261007,cloud_radius_m=20e-6,duration_s=.05,dt_s=20e-6,
               sample_interval_s=.0005,events=[{'time_s':0,'pressure_Pa':1000}])
    rows=[]
    for secondary in [False,True]:
        for n in [2,32]:
            cfg=dict(scene,count=n,physics={'secondary':secondary})
            if n==2: cfg['positions_m']=xx.tolist()
            rows.append(run_case(out/f'N{n}_secondary_{secondary}',cfg))
    write(out/'RESULT.json',dict(complete=True,checks=checks,rows=rows))


if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('--label',required=True); main(a.parse_args().label)
