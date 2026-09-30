"""Independent field-import and COM-drag verification fixtures, not clot runs."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import Delaunay
import pyvista as pv
from pd_clot.geometry import make_cloud
from pd_clot.regularization.flow import FlowSample
from pd_clot.regularization.transport import fragment_drag_half_step
from pd_clot.regularization.resolved_flow import ResolvedStreamingFieldProvider

ROOT=Path(__file__).resolve().parents[1]


class Uniform:
    def sample(self,x,t):return FlowSample(np.tile([.001,0,0],(len(x),1)),np.zeros(len(x)),np.zeros((len(x),3,3)))


def drag(out):
    c=json.loads((ROOT/'configs/streaming_regularized_demo_coarse.json').read_text());g=make_cloud(c['clot'])
    ids=np.flatnonzero(~g.fixed)[:12];cat=np.zeros(len(g.X),int);cat[ids]=1;m=g.volume*c['clot']['density_kg_m3'];M=m[ids].sum();V=g.volume[ids].sum()
    d=(6*V/np.pi)**(1/3);mu=c['pipe']['viscosity_Pa_s'];rho=c['pipe']['density_kg_m3'];T=.003
    def derivative(t,v):
        slip=.001-v[0];Re=rho*d*abs(slip)/mu
        return [3*np.pi*mu*d/M*(1+.15*Re**.687)*slip]
    ref=solve_ivp(derivative,[0,T],[0.],rtol=1e-12,atol=1e-15).y[0,-1]
    rows=[]
    for dt in [1e-5,5e-6]:
        x=g.X.copy();v=np.zeros_like(x);impulse=np.zeros(3);trajectory=[];last=None
        for k in range(round(T/dt)):
            v,diag=fragment_drag_half_step(g,x,v,m,[dict(particle_ids=ids.tolist(),fragment_id=0,attached_to_base=False)],cat,Uniform(),k*dt,dt,c)
            x[ids]+=dt*v[ids];last=diag[0];impulse+=last['integrated_drag_impulse_kg_m_s'];trajectory.append(v[ids[0],0])
        error=abs(v[ids[0],0]-ref)
        assert np.allclose(impulse,np.sum(m[:,None]*v,axis=0),atol=1e-24,rtol=1e-11)
        rows.append(dict(dt_s=dt,final_com_velocity_m_s=float(v[ids[0],0]),reference_com_velocity_m_s=float(ref),absolute_velocity_error_m_s=error,final_Re=last['Reynolds_number'],regime=last['drag_regime']))
        np.savez_compressed(out/f'drag_dt_{dt:g}.npz',time_s=np.arange(1,len(trajectory)+1)*dt,com_velocity_x_m_s=trajectory)
    assert rows[1]['absolute_velocity_error_m_s']<.6*rows[0]['absolute_velocity_error_m_s']
    result=dict(status='PASS',scope='Preclassified component COM-drag verification only; no artificial PD fracture or fragment creation',cases=rows,torque='Not modeled')
    (out/'DRAG_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')


def imported(out):
    A=np.array([[.2,.3,0],[0,-.1,.4],[.1,0,-.1]]);b=np.array([.001,.002,.003])
    xyz=np.indices((3,3,3)).reshape(3,-1).T*.0005
    cells=Delaunay(xyz).simplices
    # Delaunay on a symmetric lattice may contain zero-volume tetrahedra.
    det=np.linalg.det(xyz[cells[:,1:]]-xyz[cells[:,:1]])
    cells=cells[np.abs(det)>1e-20]
    mesh=pv.UnstructuredGrid(np.column_stack((np.full(len(cells),4),cells)).ravel(),np.full(len(cells),10,np.uint8),xyz*1000)
    mesh['velocity']=(xyz@A.T+b)*1000;mesh['pressure']=(400+xyz@np.array([100.,200.,300.]))/1000
    file=out/'affine_volume_mm.vtu';mesh.save(file)
    metadata=dict(length_unit='mm',velocity_unit='mm/s',pressure_unit='kPa',time_unit='ms',normal_convention='outward_solid',snapshot_time=0,
                  arrays=dict(velocity='velocity',pressure='pressure'),validation=dict(maximum_relative_divergence=.001))
    (out/'affine_volume_units.json').write_text(json.dumps(metadata,indent=2)+'\n')
    rng=np.random.default_rng(324);query=rng.uniform(.00001,.00099,(100,3))
    provider=ResolvedStreamingFieldProvider(file,metadata,.003,query,out/'IMPORT_DIAGNOSTIC.json');value=provider.sample(query,0)
    errors=dict(velocity=float(np.abs(value.velocity-(query@A.T+b)).max()),pressure=float(np.abs(value.pressure-(400+query@np.array([100.,200.,300.]))).max()),gradient=float(np.abs(value.gradient-A).max()))
    assert errors['velocity']<1e-12 and errors['pressure']<1e-9 and errors['gradient']<1e-10
    (out/'IMPORT_VERIFICATION.json').write_text(json.dumps(dict(status='PASS',scope='Independent affine fixture, not real microbubble CFD',maximum_absolute_SI_errors=errors),indent=2)+'\n')
    # End-to-end solver dispatch uses a separate bounded, weak, uniform fixture.
    from pd_clot.regularization.runner import run
    from pd_clot.regularization.audit import audit_run
    xyz=np.indices((2,2,2)).reshape(3,-1).T*.004-.002
    cells=Delaunay(xyz).simplices
    cube=pv.UnstructuredGrid(np.column_stack((np.full(len(cells),4),cells)).ravel(),np.full(len(cells),10,np.uint8),xyz)
    cube['velocity']=np.tile([.001,0,0],(len(xyz),1));cube['pressure']=np.ones(len(xyz))
    file=(out/'uniform_clot_domain.vtu').resolve();cube.save(file)
    meta=dict(length_unit='m',velocity_unit='m/s',pressure_unit='Pa',time_unit='s',normal_convention='outward_solid')
    meta_path=(out/'uniform_clot_units.json').resolve();meta_path.write_text(json.dumps(meta,indent=2)+'\n')
    c=json.loads((ROOT/'configs/streaming_regularized_demo_coarse.json').read_text())
    c['simulation']['number_of_macro_steps']=1;c['streaming'].update(provider='resolved_field',field_path=str(file),metadata_path=str(meta_path))
    c['verification_label']='Synthetic uniform imported-field adapter verification; not measured microstreaming'
    run(c,out/'runner_case',ROOT/'verification/regularization/calibration_coarse/CALIBRATION.json',quiet=True)
    audit_run(out/'runner_case')


def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['drag','import'],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    (drag if a.mode=='drag' else imported)(a.output)


if __name__=='__main__':main()
