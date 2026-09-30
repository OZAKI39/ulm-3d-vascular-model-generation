import argparse,csv,hashlib,json,time,traceback
from pathlib import Path
import numpy as np
from .geometry import make_cloud
from .mechanics import evaluate_cloud,prepare_shape,timestep_limit
from .damage import instantaneous,accumulate,particle_damage,fragments
from .streaming import make_provider,map_forces
from .output import write_vtk,write_json,write_collection


def simulate(c,out=None,root=None,quiet=False):
    if c['material']['model']!='neo_hookean':raise ValueError('Unsupported material; Ogden is a future extension')
    g=make_cloud(c['clot']);s=c['simulation'];dc=c['damage'];dt=s['dt_s']
    limit=timestep_limit(g,c)
    if dt<=0 or dt>limit:raise ValueError(f'Explicit dt={dt:g} s exceeds conservative estimate {limit:g} s')
    if s['number_of_macro_steps']<1 or s['output_interval']<1:raise ValueError('Positive macro steps/output interval required')
    if s['representative_cycles']<1 or s['representative_frequency_Hz']<=0:raise ValueError('Positive resolved representative cycle required')
    period=1/s['representative_frequency_Hz'];nstep=int(round(period/dt))
    if abs(nstep*dt-period)>period*1e-10:raise ValueError('dt must divide representative period exactly')
    if dc['mode'] not in ['instantaneous_stretch','cyclic_accumulation']:raise ValueError('Unknown damage mode')
    b=np.ones(len(g.pairs));x=g.X.copy();v=np.zeros_like(x)
    mass=g.volume*c['clot']['density_kg_m3'];provider=make_provider(c,root)
    shape=prepare_shape(g,b,c['safety']);t=0.;N=0;start=time.perf_counter()
    history=[];states=[];collections=[];max_accel=0.;global_minJ=1.
    if out:
        out=Path(out)
        out.mkdir(parents=True,exist_ok=False)
        write_json(out/'CONFIG.json',c)
        source_root=Path(__file__).resolve().parents[1]
        files=[*source_root.glob('pd_clot/*.py'),*source_root.glob('vendor/*.py')]
        write_json(out/'IDENTITY.json',dict(config_sha256=hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest(),
            code_sha256={str(p.relative_to(source_root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            model='3D_NOSB_PD_ONE_WAY_SYNTHETIC',dt_estimate_s=limit,particle_count=len(x),bond_count=len(b)))

    def force_state(time_s):
        mech=evaluate_cloud(g,x,b,c['material'],c['safety'],shape)
        tr=provider.traction(g.face_position,g.face_normal,time_s)
        force=mech[0]+map_forces(g,tr)
        acceleration=force/mass[:,None]
        acceleration[g.fixed]=0
        if not np.isfinite(acceleration).all():raise FloatingPointError('Nonfinite force/acceleration')
        return mech,tr,acceleration

    def record(step,q=None):
        mech,tr,_=force_state(t);D,Dmax=particle_damage(g,b)
        labels,degree,frag=fragments(g,b,dc['connectivity_threshold'])
        finiteJ=mech[2][shape[3]]
        row=dict(macro_step=step,represented_cycles=N,carrier_exposure_s=N/c['streaming']['ultrasound_frequency_Hz'],
            representative_time_s=t,mean_particle_damage=float(D.mean()),maximum_particle_damage=float(D.max()),
            damaged_particle_fraction=float(np.mean(D>dc['damaged_particle_threshold'])),
            broken_bond_fraction=float(np.mean(b==0)),maximum_displacement_um=float(np.linalg.norm(x-g.X,axis=1).max()*1e6),
            minimum_J=float(finiteJ.min()) if len(finiteJ) else None,
            strain_energy_J=float(np.dot(g.volume,mech[3])),stabilization_energy_J=float(np.dot(g.volume,mech[4])),
            kinetic_energy_J=float(.5*np.sum(mass[:,None]*v*v)),maximum_cycle_amplitude=float(q.max()) if q is not None else 0.,
            fixed_base_error_m=float(np.linalg.norm(x[g.fixed]-g.X[g.fixed],axis=1).max()),
            **frag)
        history.append(row)
        state=dict(x=x.copy(),v=v.copy(),integrity=b.copy(),damage=D.copy(),fragment_id=labels.copy(),step=step,cycles=N)
        states.append(state)
        if out and (step%s['output_interval']==0 or step==s['number_of_macro_steps']):
            fields=dict(particle_damage=D,max_bond_damage=Dmax,number_of_active_bonds=degree,
                        strain_energy_density=mech[3],stabilization_energy_density=mech[4],fragment_id=labels,
                        deformation_gradient_valid=shape[3].astype(np.uint8),J=np.nan_to_num(mech[2],nan=0.))
            write_vtk(out/'vtk',step,g,x,v,b,fields,tr);collections.append((step,N))
        return row

    try:
        record(0)
        damping=np.exp(-s['damping_s_inv']*dt/2)
        for macro in range(1,s['number_of_macro_steps']+1):
            minimum=np.full(len(b),np.inf);maximum=np.full(len(b),-np.inf)
            # Optional first-macro warmup excludes startup transients from the damage driver.
            warmup=int(s.get('warmup_cycles',1)) if macro==1 else 0
            cycles=s['representative_cycles']+warmup
            mech,tr,acc=force_state(t)
            for k in range(nstep*cycles):
                v*=damping;v+=.5*dt*acc;x+=dt*v
                x[g.fixed]=g.X[g.fixed];v[g.fixed]=0
                t+=dt
                stretch=np.linalg.norm(x[g.pairs[:,1]]-x[g.pairs[:,0]],axis=1)/g.length
                if k>=warmup*nstep:
                    minimum=np.minimum(minimum,stretch);maximum=np.maximum(maximum,stretch)
                if dc['enabled'] and dc['mode']=='instantaneous_stretch':
                    updated=instantaneous(b,stretch,dc['s1'],dc['s2'])
                    if np.any(updated!=b):b=updated;shape=prepare_shape(g,b,c['safety'])
                mech,tr,acc=force_state(t)
                peak=float(np.linalg.norm(acc,axis=1).max());max_accel=max(max_accel,peak)
                if peak>c['safety']['maximum_acceleration_m_s2']:raise FloatingPointError(f'Acceleration bound exceeded: {peak}')
                if np.linalg.norm(x-g.X,axis=1).max()>c['safety']['maximum_displacement_m']:
                    raise FloatingPointError('Displacement guard exceeded; follower loading/contact outside prototype scope')
                if np.any(shape[3]):global_minJ=min(global_minJ,float(mech[2][shape[3]].min()))
                v+=.5*dt*acc;v*=damping;v[g.fixed]=0
            q=.5*(maximum-minimum)
            if dc['enabled'] and dc['mode']=='cyclic_accumulation':
                b=accumulate(b,q,dc);shape=prepare_shape(g,b,c['safety'])
            N+=dc['DeltaN']
            row=record(macro,q)
            if not quiet:print(json.dumps({k:row[k] for k in ['macro_step','represented_cycles','mean_particle_damage','broken_bond_fraction','number_of_fragments','maximum_displacement_um']}),flush=True)
        summary=dict(status='COMPLETED_PROTOTYPE',elapsed_wall_s=time.perf_counter()-start,
            particle_count=len(x),bond_count=len(b),maximum_acceleration_m_s2=max_accel,minimum_J_all_substeps=global_minJ,
            dt_s=dt,dt_estimate_s=limit,physically_validated_thrombolysis=False,
            time_interpretation='Low-frequency mechanical proxy and independently counted carrier cycles; no resolved MHz motion',
            final=history[-1])
        if out:
            write_json(out/'SUMMARY.json',summary);write_json(out/'history.json',history)
            with (out/'summary.csv').open('w',newline='') as f:
                keys=[k for k,vv in history[0].items() if not isinstance(vv,list)]
                writer=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');writer.writeheader();writer.writerows(history)
            np.savez_compressed(out/'states.npz',X=g.X,pairs=g.pairs,fixed=g.fixed,volume=g.volume,
                x=np.array([z['x'] for z in states]),v=np.array([z['v'] for z in states]),
                integrity=np.array([z['integrity'] for z in states]),damage=np.array([z['damage'] for z in states]),
                fragment_id=np.array([z['fragment_id'] for z in states]),cycles=np.array([z['cycles'] for z in states]))
            write_collection(out,collections)
        return dict(cloud=g,x=x,v=v,integrity=b,history=history,states=states,summary=summary)
    except Exception as error:
        if out:
            write_json(out/'FAILED.json',dict(error=repr(error),representative_time_s=t,represented_cycles=N,
                maximum_acceleration_m_s2=max_accel,traceback=traceback.format_exc()))
            np.savez_compressed(out/'failed_state.npz',x=x,v=v,integrity=b,X=g.X,pairs=g.pairs)
        raise


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();simulate(json.loads(args.config.read_text()),args.output,args.config.resolve().parent)


if __name__=='__main__':main()
