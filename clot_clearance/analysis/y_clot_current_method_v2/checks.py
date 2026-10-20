"""Small implementation checks; no complete protocol or parameter search."""
from common import *
from native import native,predict
from pd_clot.lammps_fragmentation.native_y_wall_contact import CONTACT
import time

def main():
    started=time.perf_counter();case=O/'checks/native';ns,C=native(case,'biphasic')
    ref=ns.ref;model=energy_model(ref);x=ref['X'].copy();d=ns.field(3)
    assert not d.any() and not ns.field(2).any()
    s0=ns.static();independent=model.evaluate(x,d)
    assert np.max(abs(independent['force']-s0['material_force']))<1e-15
    assert np.max(abs(s0['material_force']))<1e-15 and abs(independent['U'])<1e-24
    tool=ns.tool;v=np.random.default_rng(20261019).normal(0,.003,x.shape);v[ns.n:]=0
    z=tool.evaluate(x,v,d,HOLD);active=z['active']
    assert active.any() and not np.any(z['force'][ns.n:]) and z['wall_clearance'].min()>20e-6
    ferr=float(np.max(abs(z['force'].sum(0)-z['face_force'].sum(0))))
    merr=float(np.max(abs(np.cross(x,z['force']).sum(0)-np.cross(z['centroid'],z['face_force']).sum(0))))
    perr=abs(float(np.sum(z['force']*v)-np.sum(z['face_force']*z['face_velocity'])))
    assert ferr<1e-16 and merr<1e-18 and perr<1e-16
    passive=tool.evaluate(x,v,d,END);assert passive['drag_power']<=0 and np.max(abs(passive['drive_force']))==0
    back=tool.boundary_class==2
    assert not z['active'][:tool.L][back].any()
    moved=x.copy();moved[:ns.n,2]-=10e-6
    still_attached=tool.geometry(moved,d)
    assert not still_attached['active'][:tool.L][back].any()
    released=d.copy();released[ref['bond_types']==2]=1
    free_geometry=tool.geometry(moved,released)
    released_back_count=int(free_geometry['active'][:tool.L][back].sum())
    assert released_back_count>0,'No reference attachment face became wet after genuine release and opening'
    tool.rollback_trial()
    # Distance-potential derivative probe uses the same original Y wall.
    q=tool.wall_geometry(x);probe=x.copy();probe[0]=q['nearest_wall_point_m'][0]+q['energy_normal'][0]*10e-6
    zero=np.zeros_like(x)
    def wall_energy(xx):
        qq=tool.wall_geometry(xx)
        return .5*CONTACT['stiffness_N_m']*np.dot(qq['overlap'],qq['overlap'])
    qq=tool.wall_geometry(probe);analytic=CONTACT['stiffness_N_m']*qq['overlap'][0]*qq['energy_normal'][0]
    numeric=[]
    for j in range(3):
        plus=probe.copy();minus=probe.copy();plus[0,j]+=1e-9;minus[0,j]-=1e-9
        numeric.append(-(wall_energy(plus)-wall_energy(minus))/(2e-9))
    wall_error=float(np.max(abs(np.array(numeric)-analytic)));assert wall_error<1e-9
    s=ns.dynamics();tool.commit(tool.evaluate(s['x'],s['velocity'],s['d'],s['time']))
    assert np.max(abs(s['external_force']))<1e-18
    states=[]
    for i in range(3):
        p=s;h=MIN_H;xp=p['x'].copy();xp[:ns.n]+=h*(p['velocity'][:ns.n]+h*p['force'][:ns.n]/(2*RHO*V))
        kp,dp=predict(model,ref,xp,p['kappa'],ns.field(13),ns.field(14));s=ns.step(h)
        assert np.max(abs(s['x']-xp))<2e-18 and np.max(abs(s['d']-dp))<1e-12
        assert np.array_equal(s['x'][ns.n:],ref['X'][ns.n:])
        ind=model.evaluate(s['x'],s['d']);assert np.max(abs(ind['force']-s['material_force']))<1e-12
        tool.commit(tool.evaluate(s['x'],s['velocity'],s['d'],s['time']));states.append(s)
        if i==1:
            ns.l.command(f'write_restart {case/"test.restart"}');tool.save_checkpoint(case/'surface_restart.npz',s)
    on=ns.field(13);off=ns.field(14);ns.close()
    rs,_=native(O/'checks/restart','biphasic',case/'test.restart');p=rs.dynamics();rs.tool.restore_checkpoint(case/'surface_restart.npz',p)
    assert np.array_equal(rs.field(13),on) and np.array_equal(rs.field(14),off)
    for key in ['x','velocity','d','kappa','atom_ids']:assert np.array_equal(p[key],states[1][key]),key
    continued=rs.step(MIN_H);errors={k:float(np.max(abs(continued[k]-states[2][k]))) for k in ['x','velocity','force','d','kappa']}
    assert errors['x']<2e-18 and errors['velocity']<1e-12 and errors['force']<1e-12 and errors['d']<1e-12
    rs.close()
    dump(O/'KEY_CHECKS.json',dict(status='PASS',initial_material_force_zero=True,initial_wall_contact=False,
        original_geometry_exact=True,strong_thresholds_native_exact=True,initial_L1_faces=int(active[:tool.L].sum()),
        initial_L2_faces=int(active[tool.L:].sum()),initial_area_m2=float((z['area']*z['weight']).sum()),
        pairing_force_error_N=ferr,pairing_moment_error_Nm=merr,pairing_power_error_W=perr,
        wall_potential_gradient_error_N=wall_error,restart_errors=errors,native_microsteps=4,formal_protocols=0,
        attached_back_faces_loaded=0,geometry_probe_released_open_back_faces=released_back_count,
        base_dt_s=BASE,min_dt_s=MIN_H,wall_s=time.perf_counter()-started))
    print('Y_PREFLIGHT_PASS',errors,flush=True)

if __name__=='__main__':main()

