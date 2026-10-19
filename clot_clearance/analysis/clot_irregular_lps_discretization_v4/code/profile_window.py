"""Eight accepted V3 microsteps from saved intact/damaged states; no baseline rerun."""
from pathlib import Path
import sys,os,time,json,ctypes,types,functools,contextlib,resource,argparse
import numpy as np
R=Path('/workspace/y_clot/project');TASK='clot_irregular_lps_discretization_v4';RUN='20261019T221000Z'
O=R/'verification'/TASK/RUN;OLD=R/'verification/clot_fragment_surface_coupling_v3/20261018T170400Z'
class Timers:
    def __init__(self):self.rows={};self.stack=[]
    @contextlib.contextmanager
    def scope(self,name):
        entry=[time.perf_counter(),time.process_time(),0.,0.];self.stack.append(entry)
        try:yield
        finally:
            w=time.perf_counter()-entry[0];c=time.process_time()-entry[1];self.stack.pop()
            row=self.rows.setdefault(name,[0.,0.,0]);row[0]+=w-entry[2];row[1]+=c-entry[3];row[2]+=1
            if self.stack:self.stack[-1][2]+=w;self.stack[-1][3]+=c
    def wrap(self,fn,name):
        @functools.wraps(fn)
        def call(*a,**k):
            with self.scope(name):return fn(*a,**k)
        return call
T=Timers()
def main(window,optimized=False):
    p0=O/('p1' if optimized else 'p0');p0.mkdir(exist_ok=True)
    frozen=O/'provenance/v3_source';sys.path.insert(0,str(frozen))
    import common as c
    c.O=O/'p0';c.A=frozen;c.LIB=O/'p0/build/liblammps.so'
    import current_surface as surf
    if optimized:
        sys.path.insert(0,str(R/'analysis'/TASK/'code'))
        from optimized_surface import install
        install(surf)
    for cls in [surf.CurrentSurface,surf.OldSurface]:
        for name,group in [('geometry','surface_geometry_accessibility'),('evaluate','surface_load_mapping'),('topology','surface_topology_and_components'),('commit','accepted_history_commit'),('save_checkpoint','disk_checkpoint_output'),('restore_checkpoint','snapshot_restore'),('rollback_trial','trial_rollback')]:
            if hasattr(cls,name):setattr(cls,name,T.wrap(getattr(cls,name),group))
    # Keep the exact V3 kernel. Import a saved FP64 state once, before dynamics.
    src=(frozen/'native.py').read_text()
    src=src.replace("d=self.ref['initial_cut'].astype(float);k=.30*d", "d=self.ref['initial_cut'].astype(float) if SAVED is None else SAVED['d'].copy();k=.30*d if SAVED is None else SAVED['kappa'].copy()")
    marker="            self.tool=CurrentSurface(self.ref,'uniform');"
    assert marker in src
    inject="""            if SAVED is not None:
                ids=self.l.numpy.extract_atom('id')[:self.N].astype(int)-1
                self.l.numpy.extract_atom('x')[:self.N]=SAVED['x'][ids]
                self.l.numpy.extract_atom('v')[:self.N]=SAVED['velocity'][ids]
                self.l.command(f'reset_timestep {int(SAVED["step"])} time {float(SAVED["time"]):.17g}')
"""
    src=src.replace(marker,inject+marker)
    path=OLD/'cases/biphasic/keys'/('initial.npz' if window=='intact' else 'frame_014415.npz')
    saved=c.load_npz(path)
    mod=types.ModuleType('native');mod.__dict__['SAVED']=saved;exec(compile(src,str(frozen/'native.py'),'exec'),mod.__dict__);sys.modules['native']=mod
    original=mod.native;holder={}
    def make(*args,**kwargs):
        with T.scope('initialization'):
            ns,C=original(*args,**kwargs)
            ns.lib.lpsvar_profile.argtypes=[ctypes.POINTER(ctypes.c_double),ctypes.c_int];ns.lib.lpsvar_profile.restype=ctypes.c_int
            for name,group in [('step','native_integrator_and_material'),('state','accepted_state_snapshot'),('static','trial_static_force'),('dynamics','native_setup')]:setattr(ns,name,T.wrap(getattr(ns,name),group))
            holder['native']=ns
            C.evaluate=T.wrap(C.evaluate,'independent_contact_audit')
            return ns,C
    mod.native=make;mod.predict=T.wrap(mod.predict,'damage_candidate_prediction')
    import fragment_recording
    for name in ['advance','prepare','component','set_state']:
        setattr(fragment_recording.FragmentFluidRecorder,name,T.wrap(getattr(fragment_recording.FragmentFluidRecorder,name),'fragment_fluid_diagnostics'))
    # Original run loop, shortened by accepted step count only. Output stays isolated.
    code=(frozen/'run_cases.py').read_text()
    if optimized:
        from optimized_surface import optimize_lineage_source
        code=optimize_lineage_source(code)
    code=code.replace("while prev['time']<END-1e-14:","while prev['time']<END-1e-14 and len(rows)<8:")
    code=code.replace("case=O/'cases'/name","case=PROFILE_CASE")
    code=code.replace(' component_row(s,True)'," component_row=T.wrap(component_row,'fragment_lineage_and_diagnostics')\n component_row(s,True)")
    code=code.replace("'COMPLETED' if not stop else 'NUMERICALLY_UNRESOLVED'","'PROFILE_WINDOW_COMPLETED' if not stop else 'NUMERICALLY_UNRESOLVED'")
    # Record native timers before closing the library instance.
    code=code.replace("ns.close()","capture_native(ns);ns.close()")
    ns={'__name__':'profiled_run','PROFILE_CASE':p0/'windows'/window/'biphasic','T':T}
    def capture_native(native):
        values=np.zeros(30);assert native.lib.lpsvar_profile(native.ptr(values),0)==30;holder['native_timers']=values.reshape(10,3)
    ns['capture_native']=capture_native
    exec(compile(code,str(frozen/'run_cases.py'),'exec'),ns)
    ns['save']=T.wrap(ns['save'],'disk_checkpoint_output')
    ns['graph'].topology=T.wrap(ns['graph'].topology,'material_connectivity')
    original_model=ns['energy_model']
    def model(ref):
        with T.scope('initialization'):m=original_model(ref)
        m.evaluate=T.wrap(m.evaluate,'independent_complete_energy_force');return m
    ns['energy_model']=model
    cpu=time.process_time();wall=time.perf_counter()
    with T.scope('run_loop_other'):summary=ns['run']('biphasic')
    wall=time.perf_counter()-wall;cpu=time.process_time()-cpu
    data=dict(window=window,optimized=optimized,source_state=str(path.relative_to(R)),source_sha256=c.sha(path),physical_state_import='Exact saved x,v,kappa,d and original IDs; new window cumulative counters start at zero.',summary=summary,wall_s=wall,cpu_s=cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
              python_exclusive={k:dict(wall_s=v[0],cpu_s=v[1],calls=v[2],mean_wall_s=v[0]/v[2],wall_fraction=v[0]/wall) for k,v in T.rows.items()},
              native_exclusive={k:dict(wall_s=float(v[0]),cpu_s=float(v[1]),calls=int(v[2]),mean_wall_s=float(v[0]/max(v[2],1))) for k,v in zip(['gather','evaluate','compute','damage_candidate','history_commit','sync_topology','initialize','state_output','serialize','deserialize'],holder['native_timers'])},
              accounting='Python scopes are exclusive. Native scopes subdivide opaque native calls and must NOT be added to the Python total. Prediction rejects happen before native commit; no state rollback is needed for these rejects.',timing_instrumentation='Two clocks per scope; native RAII timers. Same physical force and history expressions.')
    c.dump(p0/f'PERFORMANCE_{window.upper()}.json',data)
    print('PROFILE_COMPLETE',window,wall,summary['status'],flush=True)
    assert summary['status']=='PROFILE_WINDOW_COMPLETED',summary['stop_reasons']
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('window',choices=['intact','damaged']);p.add_argument('--optimized',action='store_true');a=p.parse_args();main(a.window,a.optimized)
