"""Ordered hard gates for the intact variational candidate; no equilibrium/phase."""
import csv,gc,hashlib,json,resource,threading,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import psutil
from pd_clot.pfpd_candidate.baqp import evaluate as historical_evaluate
from pd_clot.pfpd_candidate.energy import Material,constitutive
from pd_clot.pfpd_candidate.kernel import Family,radial
from pd_clot.pfpd_candidate.nonlocal_boundary.geometry import block_family
from pd_clot.pfpd_candidate.nonlocal_boundary.nonaffine_manufactured import QuadraticAxial,QuadraticShear
from pd_clot.pfpd_candidate.nonlocal_boundary.kinematic_energy_diagnostics import continuum_inner_energy
from pd_clot.pfpd_candidate.higher_order_rk.rk_weights import directed_topology,build_weights
from pd_clot.pfpd_candidate.higher_order_rk.rk_kinematics import gradient
from pd_clot.pfpd_candidate.higher_order_rk.mechanics_candidate import Reference,GradientOperator,VariationalBAQP,c1_operator
from pd_clot.pfpd_candidate.higher_order_rk.variational_energy import bond_state
from pd_clot.pfpd_candidate.higher_order_rk.adjoint_accumulation import delta_bond_gradient
from pd_clot.pfpd_candidate.higher_order_rk.validation_states import rotation,smooth_direction,target
from scripts.verify_pfpd_baqp import deformation as historical_deformation

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
OLD=ROOT/'verification/higher_order_baqp_rk_oracle_v1'
EARLIER=ROOT/'verification/baqp_near_boundary_kinematic_energy_audit_v1'
CFG=json.loads((OUT/'CANDIDATE_CONFIGURATION.json').read_text())
MAT=Material('svk',1000.,2000./3)
H=CFG['h_m'];L=CFG['primary_L_m'];DELTA=CFG['delta_m']

def plain(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v))
def save(n,v):(OUT/n).write_text(json.dumps(v,indent=2,default=plain,allow_nan=False)+'\n')
def csvsave(n,rows):
    with (OUT/n).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for v in iter(lambda:f.read(1024*1024),b''):h.update(v)
    return h.hexdigest()
def array_sha(a):return hashlib.sha256(np.asarray(a).tobytes()).hexdigest()
def norm(a):return np.linalg.norm(a,axis=tuple(range(1,a.ndim)))
def check(ok,message):
    if not ok:
        save('HARD_STOP.json',dict(reason=message,thresholds_unchanged=True,production_unchanged=True))
        raise RuntimeError(message)
def force_scale(force):return max(CFG['force_characteristic_stress_Pa']*H**2,float(norm(force).max()))
def energy_scale(U,V):return max(abs(U),MAT.mu*.004**2*float(V.sum()))
def metrics(d,V,scale):
    e=norm(d)
    return dict(max_absolute=float(e.max()),max_normalized=float(e.max()/scale),
        RMS_normalized=float(np.sqrt(np.average(e*e,weights=V))/scale))
def balance(name,model,n,x,force,length):
    scale=force_scale(force);net=force.sum(0);torque=np.cross(x-x.mean(0),force).sum(0)
    return dict(state=name,model=model,N=n,F_ref_N=scale,L_ref_m=length,
        force_sum_x_N=net[0],force_sum_y_N=net[1],force_sum_z_N=net[2],
        current_moment_x_Nm=torque[0],current_moment_y_Nm=torque[1],current_moment_z_Nm=torque[2],
        force_residual=float(np.linalg.norm(net)/scale),moment_residual=float(np.linalg.norm(torque)/(scale*length)))
def check_balance(rows):
    csvsave('GLOBAL_BALANCE_RESULTS.csv',rows)
    check(all(r['force_residual']<=CFG['global_force_gate'] and r['moment_residual']<=CFG['global_moment_gate'] for r in rows),'Global force/current moment balance failed')
def slices(count):
    for start in range(0,count,CFG['bond_block']):yield slice(start,start+CFG['bond_block'])

def run():
    timings=[];balances=[]
    # Verify pre-existing evidence before constructing any candidate force.
    protected=json.loads((OUT/'PROTECTED_INPUT_HASHES.json').read_text())
    check(all(sha(ROOT/p)==h for p,h in protected.items()),'Protected input changed')
    check(sha(ROOT/CFG['primary_PD2_cache'])==CFG['primary_PD2_cache_sha256'],'PD2 cache hash mismatch')
    with np.load(OLD/'REFERENCE_GEOMETRY.npz') as a:
        f=Family.from_arrays(a['X'],a['volume'],a['pairs'],DELTA,a['omega'])
        source,target_ids=a['source'],a['target'];saved_m=a['bond_weight']
    with np.load(ROOT/CFG['primary_PD2_cache']) as a:weights=a['g']
    ref=Reference.from_family(f)
    check(np.array_equal(ref.m,saved_m),'Current intact energy prefactor differs from accepted saved weights')
    op=GradientOperator('PD2',source,target_ids,weights)
    pd2=VariationalBAQP(ref,op,MAT,CFG['bond_block'],CFG['directed_block'])
    # C1 generic regression first, on a manageable physical cube with unchanged h/delta.
    _,small=block_family(CFG['small_N'],CFG['small_L_m'],DELTA)
    ss,tt,ptr=directed_topology(small.pairs,len(small.X));sref=Reference.from_family(small)
    c1small=VariationalBAQP(sref,c1_operator(small,ss,tt),MAT)
    c1states=dict(affine=historical_deformation(small.X,False),
        quadratic_shear=QuadraticShear().position(small.X))
    v=smooth_direction(small.X,CFG['small_L_m'],CFG['smooth_random_seed'])
    v*=CFG['small_random_amplitude_over_h']*H/norm(v).max()
    c1states['seeded_smooth']=small.X+v
    es={};fs={};arrays={};saved_positions={}
    for name,x in c1states.items():
        t=time.perf_counter();old=historical_evaluate(small,x,MAT);new=c1small.evaluate(x)
        difference=abs(new.energy-old.energy)/abs(old.energy)
        es[name]=dict(historical_J=old.energy,variational_J=new.energy,relative_difference=difference,
            bitwise_equal=old.energy==new.energy,PASS=difference<=CFG['C1_energy_gate'])
        scale=force_scale(old.force);m=metrics(new.force-old.force,small.volume,scale)
        fs[name]=dict(**m,F_ref_N=scale,bitwise_equal=np.array_equal(new.force,old.force),
            max_component_difference_N=float(np.abs(new.force-old.force).max()),
            PASS=m['max_normalized']<=CFG['C1_force_max_gate'] and m['RMS_normalized']<=CFG['C1_force_RMS_gate'])
        arrays.update({name+'_historical_force':old.force,name+'_variational_force':new.force,name+'_x':x})
        balances.append(balance(name,'C1',7,x,new.force,CFG['small_L_m']))
        timings.append(dict(stage='C1_regression',state=name,wall_seconds=time.perf_counter()-t))
    np.savez_compressed(OUT/'C1_FORCE_REGRESSION.npz',**arrays,X=small.X,volume=small.volume)
    save('C1_ENERGY_REGRESSION.json',es);save('C1_FORCE_REGRESSION_SUMMARY.json',fs)
    check(all(v['PASS'] for v in es.values()),'C1 energy convention regression failed; no PD2 force allowed')
    check(all(v['PASS'] for v in fs.values()),'C1 force regression failed; no PD2 force allowed')
    check_balance(balances)
    print('C1 ENERGY/FORCE PASS',max(v['relative_difference'] for v in es.values()),max(v['max_normalized'] for v in fs.values()),flush=True)
    del old,new,arrays
    # Small verification cache uses exactly the accepted local-system routine, once.
    dw=radial(np.linalg.norm(small.X[tt]-small.X[ss],axis=1)/DELTA)
    small_cache=build_weights(small.X,small.volume,ss,tt,ptr,dw,DELTA,9,CFG['workers'])
    check(np.all(small_cache['rank']==9),'Small verification support rank deficient')
    np.savez_compressed(OUT/'N7_PD2_REFERENCE_CACHE.npz',**small_cache)
    np.savez_compressed(OUT/'N7_REFERENCE_GEOMETRY.npz',X=small.X,volume=small.volume,pairs=small.pairs,
        omega=small.omega,source=ss,target=tt,ptr=ptr)
    smallop=GradientOperator('PD2',ss,tt,small_cache['g']);p7=VariationalBAQP(sref,smallop,MAT)
    small_cache_info=dict(N=7,nodes=len(small.X),bonds=len(small.pairs),min_rank=int(small_cache['rank'].min()),
        max_condition=float(small_cache['condition'].max()),cache_sha256=sha(OUT/'N7_PD2_REFERENCE_CACHE.npz'),
        geometry_sha256=sha(OUT/'N7_REFERENCE_GEOMETRY.npz'),used_in='tiny delta-Ftilde / selected-DOF implementation oracles only')
    del small_cache
    x=c1states['quadratic_shear'];v=smooth_direction(small.X,CFG['small_L_m'],CFG['direction_seeds'][0]);eps=H*CFG['delta_Ftilde_step_rel']
    dF=gradient(v,ss,tt,smallop.g);F=p7.nodes(x);Fp=p7.nodes(x+eps*v);Fm=p7.nodes(x-eps*v)
    sl=slice(None);analytic=delta_bond_gradient(sref,v,dF,sl)
    plus=bond_state(sref,x+eps*v,Fp,sl,MAT)[0];minus=bond_state(sref,x-eps*v,Fm,sl,MAT)[0]
    fd=(plus-minus)/(2*eps);scale=float(norm(analytic).max())
    identity=dict(step_m=eps,max_absolute_residual_per_m=float(norm(fd-analytic).max()),
        derivative_scale_per_m=scale,max_normalized_residual=float(norm(fd-analytic).max()/scale),gate=CFG['delta_Ftilde_gate'])
    identity['PASS']=identity['max_normalized_residual']<=CFG['delta_Ftilde_gate']
    save('DELTA_FTILDE_IDENTITY.json',identity)
    check(identity['PASS'],'Analytic delta-Ftilde identity failed')
    tensors=np.array([np.eye(3),[[1,.03,0],[0,1,0],[0,0,1]],np.diag([1.02,.99,1.01]),rotation()])
    _,P=constitutive(tensors,MAT);fdP=np.zeros_like(P);e=CFG['constitutive_F_step']
    for a in range(3):
        for b in range(3):
            plus=tensors.copy();minus=tensors.copy();plus[:,a,b]+=e;minus[:,a,b]-=e
            fdP[:,a,b]=(constitutive(plus,MAT)[0]-constitutive(minus,MAT)[0])/(2*e)
    scale=max(CFG['force_characteristic_stress_Pa'],float(norm(P).max()))
    constit=dict(tensors=tensors,P=P,central_FD_P=fdP,step=e,max_absolute_error_Pa=float(norm(fdP-P).max()),
        max_normalized_error=float(norm(fdP-P).max()/scale),gate=CFG['constitutive_gate'])
    constit['PASS']=constit['max_normalized_error']<=CFG['constitutive_gate']
    save('CONSTITUTIVE_STRESS_CHECK.json',constit);check(constit['PASS'],'Existing constitutive stress FD mismatch')
    # Full N19 cached operator regression against accepted full-population arrays.
    kinematics={};states={}
    for name in ['axial','shear','cross']:
        x,exactF=target(name,f.X,L);actual=pd2.nodes(x)
        with np.load(OLD/(name.upper()+'_NODAL_STATE.npz')) as a:
            m=metrics(actual-a['PD2'],f.volume,.004)
            target_same=np.array_equal(x,a['x']);exact_same=np.array_equal(exactF,a['F_exact'])
        truth=metrics(actual-exactF,f.volume,.004)
        kinematics[name]=dict(previous_PD2_difference=m,analytic_error=truth,
            target_bitwise_equal=target_same,analytic_F_bitwise_equal=exact_same)
        check(target_same and exact_same and m['max_normalized']<=CFG['kinematic_max_gate'] and
            truth['max_normalized']<=CFG['kinematic_max_gate'] and truth['RMS_normalized']<=CFG['kinematic_RMS_gate'],'Cached PD2 oracle regression failed')
    save('PD2_KINEMATIC_REGRESSION.json',kinematics)
    save('PD2_CACHE_IDENTITY.json',dict(status='IDENTICAL',primary_cache_sha256=sha(ROOT/CFG['primary_PD2_cache']),
        g_array_sha256=array_sha(op.g),geometry_sha256=sha(OLD/'REFERENCE_GEOMETRY.npz'),
        old_configuration_sha256=sha(OLD/'PD2_CONFIGURATION.json'),bond_prefactors_bitwise_equal=True,
        reference_arrays_writeable=False,small_verification_cache=small_cache_info))
    print('KINEMATIC / DELTA-FTILDE / CONSTITUTIVE PASS',identity['max_normalized_residual'],flush=True)
    # Now and only now evaluate primary PD2 forces on all N19 physical nodes.
    names=['finite_affine_shear','axial','shear','cross','affine_axial','affine_shear']
    def mechanical(name):
        x,Fexact=target(name,f.X,L);t=time.perf_counter();state=pd2.evaluate(x)
        return name,x,state,time.perf_counter()-t
    with ThreadPoolExecutor(CFG['workers']) as pool:
        for name,x,state,elapsed in pool.map(mechanical,names):
            states[name]=(x,state);balances.append(balance(name,'PD2',19,x,state.force,L))
            timings.append(dict(stage='PD2_energy_force',state=name,wall_seconds=elapsed))
            np.savez_compressed(OUT/(name.upper()+'_STATE.npz'),x=x,Fbar=state.Fbar,force=state.force,energy_J=state.energy)
    check_balance(balances)
    # Rigid-state gates plus the historical superposed-rotation covariance check.
    Q=rotation(CFG['rotation_angle_degrees'],CFG['rotation_axis']);c=np.array(CFG['rigid_translation_m'])
    objectivity={}
    for name,x,Fexact in [('translation',f.X+c,np.eye(3)),('rotation',f.X@Q.T+c,Q)]:
        t=time.perf_counter();state=pd2.evaluate(x);node=float(norm(state.Fbar-Fexact).max());bond=0.;strain=0.
        for sl in slices(len(f.pairs)):
            Ft=bond_state(ref,x,state.Fbar,sl,MAT)[0]
            bond=max(bond,float(norm(Ft-Fexact).max()))
            strain=max(strain,float(norm(.5*(Ft.swapaxes(-1,-2)@Ft-np.eye(3))).max()))
        maxforce=float(norm(state.force).max());ffloor=max(1e-12*CFG['force_characteristic_stress_Pa']*H**2,1e-15)
        row=dict(energy_J=state.energy,max_force_N=maxforce,node_F_max_error=node,bond_F_max_error=bond,
            Green_strain_max=strain,force_gate_N=ffloor,force_over_characteristic=maxforce/(4*H**2),
            energy_over_characteristic=state.energy/(MAT.mu*.004**2*sum(f.volume)))
        row['PASS']=node<=CFG['rigid_F_max_gate'] and bond<=CFG['rigid_F_max_gate'] and strain<=CFG['rigid_strain_max_gate'] and state.energy<=CFG['rigid_energy_absolute_gate_J'] and maxforce<=ffloor
        objectivity[name]=row;save('OBJECTIVITY_RESULTS.json',objectivity)
        check(row['PASS'],'PD2 rigid '+name+' failed')
        timings.append(dict(stage='PD2_objectivity',state=name,wall_seconds=time.perf_counter()-t))
    x,base=states['finite_affine_shear'];rotated=pd2.evaluate(x@Q.T+c)
    er=abs(rotated.energy-base.energy)/base.energy;fr=metrics(rotated.force-base.force@Q.T,f.volume,force_scale(base.force))
    objectivity['superposed_rotation']=dict(energy_relative_difference=er,force_covariance=fr,PASS=er<=1e-10 and fr['max_normalized']<=1e-10)
    save('OBJECTIVITY_RESULTS.json',objectivity);check(objectivity['superposed_rotation']['PASS'],'Deformed-state objectivity/covariance failed')
    balances.append(balance('superposed_finite_shear','PD2',19,x@Q.T+c,rotated.force,L));check_balance(balances)
    print('PD2 OBJECTIVITY AND GLOBAL BALANCE PASS',flush=True)
    # Same kernel, volumes and exact m_b: affine and quadratic energy checks.
    c1=VariationalBAQP(ref,c1_operator(f,source,target_ids),MAT)
    affines={}
    for name in ['affine_axial','affine_shear']:
        x,pstate=states[name];Fc=c1.nodes(x);Uc=c1.energy(x);diffF=0.;diffW=0.;Wref=0.
        for sl in slices(len(f.pairs)):
            Ft1,W1,_=bond_state(ref,x,Fc,sl,MAT);Ft2,W2,_=bond_state(ref,x,pstate.Fbar,sl,MAT)
            diffF=max(diffF,float(norm(Ft1-Ft2).max()));diffW=max(diffW,float(np.abs(W1-W2).max()));Wref=max(Wref,float(np.abs(W1).max()))
        er=abs(pstate.energy-Uc)/abs(Uc)
        row=dict(U_C1_J=Uc,U_PD2_J=pstate.energy,relative_energy_difference=er,
            bond_F_max_difference=diffF,bond_F_difference_over_0p004=diffF/.004,
            bond_density_max_difference_Pa=diffW,bond_density_relative_difference=diffW/Wref)
        row['PASS']=er<=CFG['affine_energy_gate'] and diffF/.004<=CFG['affine_bond_F_gate'] and diffW/Wref<=CFG['affine_bond_density_gate']
        affines[name]=row
    save('AFFINE_ENERGY_COMPARISON.json',affines);check(all(r['PASS'] for r in affines.values()),'Affine C1-PD2 energy or bondwise regression failed')
    oldenergies=json.loads((EARLIER/'ENERGY_SOURCE_DECOMPOSITION.json').read_text());quad=[]
    for name in ['axial','shear','cross']:
        x,pstate=states[name];Uc=c1.energy(x);U1=0.
        for sl in slices(len(f.pairs)):
            i,j=f.pairs[sl].T;mid=.5*(f.X[i]+f.X[j]);Fexact=target(name,mid,L)[1]
            U1+=float(np.sum(ref.m[sl]*constitutive(Fexact,MAT)[0]))
        if name!='cross':
            field=QuadraticAxial() if name=='axial' else QuadraticShear()
            U2=continuum_inner_energy(field,MAT,L,0,32)[0];U264=continuum_inner_energy(field,MAT,L,0,64)[0]
            saved=oldenergies[name]['intact']
            historical={k:abs(v-saved[k])/energy_scale(saved[k],f.volume) for k,v in [('U0_J',Uc),('U1_J',U1),('U2_J',U2)]}
            check(max(historical.values())<=CFG['quadratic_U1_gate'],'Historical energy decomposition regression failed')
        else:
            # XZ field is Y-independent; tensor Gauss in X,Z times exact Y length.
            def gauss(order):
                z,w=np.polynomial.legendre.leggauss(order);z*=L/2;w*=L/2
                a,b=np.indices((order,order)).reshape(2,-1);X=np.column_stack([z[a],np.zeros(len(a)),z[b]])
                return float(np.sum(w[a]*w[b]*L*constitutive(target('cross',X,L)[1],MAT)[0]))
            U2,U264=gauss(32),gauss(64);historical=None
        scale=energy_scale(U1,f.volume);er=abs(pstate.energy-U1)/scale
        row=dict(field=name,U_C1_J=Uc,U_PD2_J=pstate.energy,U_exactkin_currentweights_J=U1,U_continuum_J=float(U2),
            kinematic_correction_J=pstate.energy-Uc,remaining_weight_quadrature_difference_J=pstate.energy-U2,
            U1_minus_U2_J=U1-U2,energy_scale_J=scale,PD2_vs_U1_normalized_error=er,
            Gauss32_vs64_relative=float(abs(U2-U264)/abs(U2)),PASS=er<=CFG['quadratic_U1_gate'])
        quad.append(row)
    csvsave('QUADRATIC_ENERGY_DECOMPOSITION.csv',quad)
    check(all(r['PASS'] and r['Gauss32_vs64_relative']<=1e-12 for r in quad),'PD2 energy vs exact-current-weight reference failed')
    print('AFFINE AND QUADRATIC ENERGY PASS',flush=True)
    # Exactly six directions; fixed physical epsilon = eps_rel*h. Save all 18 samples.
    def fd_case(arg):
        name,number,seed=arg;x,st=states[name];v=smooth_direction(f.X,L,seed)
        Dforce=float(-np.sum(st.force*v));Eref=energy_scale(st.energy,f.volume);out=[];t=time.perf_counter()
        for k,eps_rel in enumerate(CFG['FD_eps_rel']):
            eps=eps_rel*H;up=pd2.energy(x+eps*v);um=pd2.energy(x-eps*v);Dfd=(up-um)/(2*eps)
            denominator=max(abs(Dfd),abs(Dforce),Eref/H*np.linalg.norm(v))
            err=abs(Dfd-Dforce)/denominator;raw=abs(Dfd-Dforce)/max(abs(Dfd),abs(Dforce),1e-300)
            out.append(dict(state=name,direction=number,seed=seed,eps_rel=eps_rel,eps_m=eps,U_plus_J=up,U_minus_J=um,
                D_force_N=Dforce,D_fd_N=Dfd,E_ref_J=Eref,direction_norm=float(np.linalg.norm(v)),denominator_N=float(denominator),
                normalized_error=err,unfloored_relative_error=raw,middle_step=k==1,
                middle_gate_PASS=(err<=CFG['FD_middle_gate']) if k==1 else '',
                gross_gate_PASS=err<=CFG['FD_all_steps_gross_gate']))
        np.savez_compressed(OUT/f'FD_{name}_{number}_DIRECTION.npz',v=v)
        return out,dict(stage='FD_six_energy_calls',state=name+'_direction_'+str(number),wall_seconds=time.perf_counter()-t)
    tasks=[(name,k,seed) for name in CFG['FD_states'] for k,seed in enumerate(CFG['direction_seeds'])]
    fdrows=[]
    with ThreadPoolExecutor(CFG['workers']) as pool:
        for r,t in pool.map(fd_case,tasks):fdrows.extend(r);timings.append(t)
    csvsave('ENERGY_FORCE_DIRECTIONAL_DERIVATIVES.csv',fdrows)
    fdsummary=[]
    for name,k,seed in tasks:
        rr=[r for r in fdrows if r['state']==name and r['direction']==k];errors=[r['normalized_error'] for r in rr]
        status='PASS' if errors[1]<=CFG['FD_middle_gate'] and max(errors)<=CFG['FD_all_steps_gross_gate'] else 'FAIL'
        fdsummary.append(dict(state=name,direction=k,middle_normalized_error=errors[1],all_step_errors=errors,
            sequence='BELOW_PREDECLARED_NOISE_PLATEAU' if max(errors)<=CFG['FD_noise_plateau_gate'] else 'NO_GROSS_NONCONVERGENCE' if status=='PASS' else 'FAILED_SEQUENCE',status=status))
    save('ENERGY_FORCE_FD_SUMMARY.json',fdsummary)
    check(all(r['status']=='PASS' for r in fdsummary),'Energy-force FD failed; no changes to steps/weights/tolerances allowed')
    print('SIX DIRECTIONAL FD CHECKS PASS',max(r['normalized_error'] for r in fdrows if r['middle_step']),flush=True)
    # <=12 tiny coordinate derivatives, same fixed middle step; never build a Jacobian.
    x=target('cross',small.X,CFG['small_L_m'])[0];state=p7.evaluate(x);scale=force_scale(state.force)
    mid=CFG['small_L_m']/2
    targets=np.array([[0,0,0],[0,0,-mid],[-mid,0,-mid],[-mid,-mid,-mid]])
    ids=[int(np.linalg.norm(small.X-q,axis=1).argmin()) for q in targets];selected=[];eps=H*CFG['selected_DOF_step_rel']
    for node in ids:
        for component in range(3):
            xp=x.copy();xm=x.copy();xp[node,component]+=eps;xm[node,component]-=eps
            up,um=p7.energy(xp),p7.energy(xm);slow=-(up-um)/(2*eps);fast=state.force[node,component]
            selected.append(dict(node=node,component=component,eps_m=eps,U_plus_J=up,U_minus_J=um,
                slow_force_N=slow,adjoint_force_N=fast,F_ref_N=scale,normalized_error=abs(slow-fast)/scale))
    csvsave('SELECTED_DOF_FORCE_ORACLE.csv',selected)
    np.savez_compressed(OUT/'N7_SELECTED_STATE.npz',x=x,Fbar=state.Fbar,force=state.force,energy_J=state.energy)
    balances.append(balance('cross_selected_DOF','PD2',7,x,state.force,CFG['small_L_m']));check_balance(balances)
    check(len(selected)==12 and max(r['normalized_error'] for r in selected)<=CFG['selected_DOF_gate'],'Selected-DOF force derivative oracle failed')
    # Inputs and prior scientific products remain immutable.
    changed=[p for p,h in protected.items() if sha(ROOT/p)!=h]
    save('HISTORICAL_INTEGRITY.json',dict(protected=len(protected),changed=changed,PASS=not changed))
    check(not changed,'Historical inputs changed during candidate audit')
    save('CALL_TIMINGS.json',timings)
    save('FINAL_STATUS.json',dict(PD2_KINEMATICS='REUSED_VERIFIED',DELTA_FTILDE_IDENTITY='PASS',
        C1_VARIATIONAL_ENERGY_REGRESSION='PASS',C1_VARIATIONAL_FORCE_REGRESSION='PASS',PD2_RIGID_TRANSLATION='PASS',
        PD2_RIGID_ROTATION='PASS',PD2_GLOBAL_FORCE_BALANCE='PASS',PD2_GLOBAL_MOMENT_BALANCE='PASS',
        PD2_ENERGY_FORCE_FD='PASS',PD2_SELECTED_DOF_FORCE_ORACLE='PASS',PD2_AFFINE_ENERGY_REGRESSION='PASS',
        PD2_QUADRATIC_ENERGY_KINEMATIC_REFERENCE='PASS',PD2_OPTIONAL_MANUFACTURED_RECOVERY='NOT_RUN',
        HIGHER_ORDER_VARIATIONAL_MECHANICS_CANDIDATE='VERIFIED',PRODUCTION_BAQP='UNCHANGED',
        BOUNDARY_DATA_CONVERSION='NOT_REINTEGRATED_THIS_ROUND',PHASE_FORCE='NOT_IMPLEMENTED',C5_C7='NOT_RUN',
        FRACTURE='NOT_RUN',CLOT_GEOMETRY='NOT_RUN'))
    print('HIGHER_ORDER_VARIATIONAL_MECHANICS_CANDIDATE = VERIFIED',flush=True)

def main():
    stop=threading.Event();mon=dict(peak_process_tree_RSS_bytes=0,minimum_available_RAM_fraction=1.)
    def monitor():
        while not stop.is_set():
            proc=psutil.Process();rss=proc.memory_info().rss
            for child in proc.children(recursive=True):
                try:rss+=child.memory_info().rss
                except psutil.Error:pass
            mon['peak_process_tree_RSS_bytes']=max(mon['peak_process_tree_RSS_bytes'],rss)
            m=psutil.virtual_memory();mon['minimum_available_RAM_fraction']=min(mon['minimum_available_RAM_fraction'],m.available/m.total)
            stop.wait(.1)
    thread=threading.Thread(target=monitor,daemon=True);thread.start();start=time.perf_counter()
    try:run()
    finally:
        stop.set();thread.join()
        mon.update(wall_seconds=time.perf_counter()-start,main_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            CPU_workers=CFG['workers'],BLAS_threads=1,GPU='GPU_NOT_USED_BY_DESIGN',
            worker_model='3 independent state/direction threads; shared immutable reference arrays; no nested pools',
            full_Hessian_formed=False,equilibrium_run=False)
        save('PERFORMANCE_SUMMARY.json',mon)
    check(mon['minimum_available_RAM_fraction']>=CFG['RAM_headroom_min'],'RAM headroom below predeclared minimum')

if __name__=='__main__':main()
