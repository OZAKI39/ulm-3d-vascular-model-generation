"""Reuse the validated native VV/material; add only the follower force callback."""
from common import *
from current_surface import CurrentSurface
import traceback,ast,ctypes,types
from mpi4py import MPI
from lammps import lammps

def make_case(case):
    case=Path(case);(case/'provenance').mkdir(parents=True,exist_ok=True);(case/'keys').mkdir(exist_ok=True)
    for name in ['reference.npz','reference.native','body.data']:
        p=case/'provenance'/name
        if not p.exists():p.symlink_to(os.path.relpath(O/'provenance'/name,p.parent))
    return case

def factory(case):
    source=R/'analysis/v3_softening_pilot_v1/native.py';tree=ast.parse(source.read_text())
    classes=[x for x in tree.body if isinstance(x,ast.ClassDef) and x.name in ['Native','Contact']]
    ns=dict(globals(),O=Path(case));exec(compile(ast.Module(body=classes,type_ignores=[]),str(source),'exec'),ns)
    return ns['Native'],ns['Contact']

def native(case,protocol='uniform',restart=None):
    case=make_case(case);Base,Contact=factory(case)
    class Native(Base):
        def __init__(self):
            super().__init__('native',restart=restart,dt=BASE)
            if restart is None:
                d=self.ref['initial_cut'].astype(float);k=SC[self.ref['bond_types']]*d
                assert self.lib.lpsvar_action(self.l.lmp,1,0,self.ptr(k),self.ptr(d))==0
            self.tool=CurrentSurface(self.ref,protocol);self.ext=types.SimpleNamespace(error=None,last={});self.force_stage_dt=0.

        def step(self,h):
            if h!=self.dt:self.l.command(f'timestep {h:.17g}');self.dt=h
            self.force_stage_dt=h
            try:self.l.command('run 1 pre no post no')
            finally:self.force_stage_dt=0.
            return self.state()

        def dynamics(self):
            self.l.commands_string(f'group clot type 1\nfix integrate clot nve\nfix coordinator all LPS/VARIATIONAL/DYNAMIC {RHO:.17g}\nfix v3forces clot external pf/callback 1 1')
            def callback(caller,step,nlocal,tags,x,f):
                f[:]=0
                try:
                    ids=tags.astype(int)-1
                    assert nlocal==self.N and np.array_equal(np.sort(ids),np.arange(self.N))
                    current=np.empty((self.N,3));velocity=np.empty_like(current)
                    current[ids]=x;velocity[ids]=self.l.numpy.extract_atom('v')[:self.N]
                    t=float(self.l.get_thermo('time'));d=self.field(3);drift=velocity.copy();iterations=0
                    if self.force_stage_dt:
                        # Implicit velocity-dependent second half-kick, entirely
                        # inside the external callback. Native NVE is unchanged.
                        internal=np.empty_like(current);internal[ids]=self.l.numpy.extract_atom('f')[:self.N]
                        factor=self.force_stage_dt/(2*RHO*V)
                        for iterations in range(1,13):
                            z=self.tool.evaluate(current,velocity,d,t)
                            updated=drift.copy();updated[:self.n]+=factor*(internal[:self.n]+z['force'][:self.n])
                            err=np.max(abs(updated-velocity));velocity=updated
                            if err<2e-13*max(1.,float(np.max(abs(velocity)))):break
                        else:raise RuntimeError('Follower slip half-kick did not converge')
                    z=self.tool.evaluate(current,velocity,d,t);f[:]=z['force'][ids]
                    self.ext.last=dict(external_force=z['force'],external_drive_force=z['drive_force'],
                        external_drag_force=z['drag_force'],external_time=np.array(t),
                        force_stage_velocity=velocity,verlet_drift_velocity=drift,slip_iterations=np.array(iterations),surface_centroid=z['centroid'],
                        surface_normal=z['normal'],surface_area=z['area'],surface_active=z['active'],
                        surface_velocity=z['face_velocity'],surface_outer_velocity=z['outer_velocity'],
                        surface_traction=z['traction'],surface_force=z['face_force'])
                except BaseException:
                    self.ext.error=traceback.format_exc();f[:]=np.nan
            self.callback=callback;self.l.set_fix_external_callback('v3forces',callback)
            self.dynamic=True;self.l.command('run 0');return self.state()
    return Native(),Contact

# Import only the unchanged predictor's function definition, not an old runner.
_source=PA/'native.py';_tree=ast.parse(_source.read_text())
exec(compile(ast.Module(body=[n for n in _tree.body if isinstance(n,ast.FunctionDef) and n.name=='predict'],type_ignores=[]),str(_source),'exec'))
