"""Current-step native topology, event-triggered COM drag, auditable SI forces.

No material modification, native hooks, prescribed velocities or future labels.
"""
import json, time, os, sys, shutil, subprocess, traceback, hashlib, re
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .config import ROOT, dump
from .y_clot_morphology_v3 import OUT as PRIOR, SOURCE, LIBDIR, DEPS, load_geometry
from .point_cloud_delivery import OUT, RESULTS
from .native_y_wall_contact import WallCallback, CONTACT
from .detached_fragment_metrics import graph_components
from .y_bifurcation_sources import digest
from .failure_map import csv_out

CASE=OUT/'cases/C1_GAUSSIAN_A200'
OLD=PRIOR/'cases/C1_GAUSSIAN_A200'
DT=.5e-6
MASS=1000*(50e-6)**3
TAU=100e-6
DTYPE=np.dtype([(k,'<i4') for k in ['step','atom_id','component_id','node_count','detach_step','valid_CFD','near_wall_uncertain','clot_present_flow_uncertain']]+[(k,'<f8',(3,)) for k in ['x_m','v_force_half_m_s','u_COM_m_s','drag_N','numerical_N','wall_N','load_N','force_before_flow_N']]+[('tau_s','<f8'),('tetra_id','<i4')])

def drag_law(count,velocity,fluid):
    """Identical instantaneous COM force in native and long transport stages."""
    count=np.asarray(count);v=np.asarray(velocity);u=np.asarray(fluid)
    diameter=(6*count*(50e-6)**3/np.pi)**(1/3)
    reynolds=1000*diameter*np.linalg.norm(v-u,axis=-1)/.001
    factor=np.where(reynolds<1000,1+.15*reynolds**.687,.44*reynolds/24)
    tau=np.where(count<=2,TAU,1000*diameter**2/(18*.001*factor))
    return MASS*count[...,None]*(u-v)/tau[...,None],tau,reynolds

class EarlyFlow:
    def __init__(self,lmp,comm,path,wall):
        self.lmp=lmp;self.comm=comm;self.path=path;self.wall=wall;self.g=load_geometry();self.n=self.g['clot_count'];self.N=len(self.g['X'])
        self.b=self.g['bonds'];self.types=self.g['bond_types'];self.keys=self.b[:,0]*self.N+self.b[:,1];self.alive=np.ones(len(self.b),bool)
        self.nbonds=-1;self.detached=np.zeros(self.n,bool);self.detach_step=np.full(self.n,-1,int);self.drag_step=np.full(self.n,-1,int)
        self.events=[];self.topology_events=[];self.history=[];self.calls=0;self.seconds=0;self.max_detached=0;self.max_component=0
        self.r=json.loads((path/'request.json').read_text());self.load=json.loads((path/'LOAD_PROVIDER.json').read_text())
        from .y_fragmentation_cases import tangent
        self.reference_tangent=tangent(self.g);self.ref_tree=cKDTree(self.g['X'][:self.n])
        self.file=(path/'detached_force_audit.bin').open('wb') if comm.rank==0 else None
        self.fullfile=(path/'saved_frame_force_audit.bin').open('wb') if comm.rank==0 else None
        self.audit_count=0;self.snapshot_count=0
        if comm.rank==0:dump(path/'FORCE_AUDIT_DTYPE.json',dict(descr=DTYPE.descr,bytes_per_row=DTYPE.itemsize,velocity_time_level='Native post_force uses Verlet half-step velocity; trajectory dump after final_integrate has full-step velocity.',query_position='Each current component COM; u repeated for member force distribution.',force_order='native SI conversion, Gaussian load, native WALL, flow + attached-only numerical damping, anchor setforce'))
    def callback(self,caller,step,nlocal,tags,x,f):
        try:self.compute(int(step),nlocal,tags,x,f)
        except BaseException:
            (self.path/f'EARLY_FLOW_ERROR_rank{self.comm.rank}.txt').write_text(traceback.format_exc());self.comm.Abort(41)
    def update_topology(self,step):
        count=int(self.lmp.extract_global('nbonds'))
        if count==self.nbonds:return False
        nb,raw=self.lmp.gather_bonds();a=np.ctypeslib.as_array(raw).reshape(-1,3)
        assert nb==count and (a[:,0]>0).all()
        ends=np.sort(a[:,1:3].astype(np.int64)-1,axis=1);key=ends[:,0]*self.N+ends[:,1];idx=np.searchsorted(self.keys,key)
        assert len(np.unique(key))==count and np.array_equal(self.keys[idx],key) and np.array_equal(self.types[idx],a[:,0])
        alive=np.zeros(len(self.b),bool);alive[idx]=True;assert not np.any(alive&~self.alive)
        broken=np.flatnonzero(self.alive&~alive)
        if self.nbonds<0:assert np.array_equal(alive,self.alive)
        else:self.topology_events.extend((step,int(i)) for i in broken)
        self.alive=alive;self.nbonds=count
        lab,sizes,attached=graph_components(self.n,self.b,self.types,alive);self.lab=lab;self.sizes=sizes;self.attached=attached
        self.cid=np.full(len(sizes),self.n+1,int);np.minimum.at(self.cid,lab,np.arange(1,self.n+1))
        det=~attached[lab];assert not np.any(self.detached&~det),'Removed bonds cannot reattach'
        new=det&~self.detached;self.detach_step[new]=step;self.detached=det
        self.members=[np.flatnonzero(lab==j) for j in np.flatnonzero(~attached)]
        self.free_components=np.flatnonzero(~attached)
        self.max_detached=max(self.max_detached,int(det.sum()));self.max_component=max(self.max_component,int(sizes[~attached].max(initial=0)))
        return True
    def compute(self,step,nlocal,tags,x,f):
        started=time.perf_counter();changed=self.update_topology(step);f[:]=0
        ids=tags.astype(int)-1;local=np.flatnonzero(ids<self.n);ii=ids[local]
        vel=self.lmp.numpy.extract_atom('v')[:nlocal]
        # Every member's current x, half-step v, actual wall F and pre-flow F.
        owned=np.zeros((self.n,13));owned[ii,:3]=x[local];owned[ii,3:6]=vel[local]
        owned[ii,6:9]=self.wall.last_force[local];owned[ii,9:12]=self.lmp.numpy.extract_atom('f')[:nlocal][local];owned[ii,12]=self.wall.clearance[ii]
        global_state=np.empty_like(owned);self.comm.Allreduce(owned,global_state)
        drag=np.zeros((self.n,3));numeric=-500*MASS*global_state[:,3:6];numeric[self.detached]=0
        if self.comm.rank==0:
            X=global_state[:,:3];V=global_state[:,3:6];u=np.zeros_like(X);valid=np.zeros(self.n,bool);tet=np.full(self.n,-1,int);tau=np.zeros(self.n);near=np.zeros(self.n,bool);uncertain=np.zeros(self.n,bool)
            if self.members:
                pos=np.array([X[z].mean(0) for z in self.members]);v=np.array([V[z].mean(0) for z in self.members]);counts=np.array([len(z) for z in self.members])
                sample=self.wall.query.flow.sample(pos);okay=sample['valid'];F=np.zeros_like(pos);tt=np.zeros(len(pos));reynolds=np.zeros(len(pos))
                F[okay],tt[okay],reynolds[okay]=drag_law(counts[okay],v[okay],sample['velocity_m_s'][okay])
                retained=cKDTree(X[~self.detached]) if (~self.detached).any() else None
                for j,z in enumerate(self.members):
                    valid[z]=okay[j];tet[z]=sample['tetra_id'][j];tau[z]=tt[j]
                    near[z]=np.min(global_state[z,12])<100e-6
                    uncertain[z]=(self.ref_tree.query(pos[j])[0]<=np.sqrt(3)*25e-6 or (retained is not None and retained.query(pos[j])[0]<=50e-6))
                    if okay[j]:
                        u[z]=sample['velocity_m_s'][j];drag[z]=F[j]/len(z)
                        first=z[self.drag_step[z]<0];self.drag_step[first]=step
            # Only actual detached valid groups receive flow; attached damping only.
            assert not np.any(drag[~self.detached]) and not np.any(numeric[self.detached])
            T=self.r['duration_s'];t=step*DT
            amp=.5*(1-np.cos(4*np.pi*t/T)) if t<T/4 else 1. if t<.75*T else .5*(1+np.cos(4*np.pi*(t/T-.75))) if t<T else 0.
            env=np.exp(-np.sum((X-self.r['load_centers_m'][0])**2,axis=1)/(2*self.r['gaussian_sigma_m']**2))
            load=self.g['areas'][:self.n,None]*self.load['scale_Pa']*amp*env[:,None]*self.reference_tangent[:self.n]
            data=np.zeros(self.n,DTYPE);data['step']=step;data['atom_id']=np.arange(1,self.n+1);data['component_id']=self.cid[self.lab];data['node_count']=self.sizes[self.lab];data['detach_step']=self.detach_step;data['valid_CFD']=valid;data['near_wall_uncertain']=near;data['clot_present_flow_uncertain']=uncertain
            for key,values in [('x_m',X),('v_force_half_m_s',V),('u_COM_m_s',u),('drag_N',drag),('numerical_N',numeric),('wall_N',global_state[:,6:9]),('load_N',load),('force_before_flow_N',global_state[:,9:12]),('tau_s',tau),('tetra_id',tet)]:data[key]=values
            data[self.detached].tofile(self.file);self.audit_count+=int(self.detached.sum())
            if step%20==0 or step==4660 or changed:
                data.tofile(self.fullfile);self.snapshot_count+=self.n
            # Exact first break/detach states for static frames, no 10us rounding.
            if changed:
                if self.topology_events and not (self.path/'first_break_state.npz').exists():self.save_event('first_break',step,global_state)
                if self.detached.any() and not (self.path/'first_detachment_state.npz').exists():self.save_event('first_detachment',step,global_state)
            self.history.append(dict(step=step,time_s=t,native_bonds=self.nbonds,topology_updated=changed,detached_nodes=int(self.detached.sum()),valid_drag_nodes=int(valid.sum()),invalid_detached_nodes=int((self.detached&~valid).sum()),drag_on_attached_nodes=0,numerical_damping_on_detached_nodes=0,force_law='COM mass-weighted distribution; tau100us for1-2; SN sphere for>=3'))
        self.comm.Bcast(drag,root=0)
        f[local]=drag[ii]+numeric[ii]
        # Separate diagnostic per-atom fix arrays; never touch coordinates/velocity.
        for prefix,values in [('drag',drag[ii]),('num',numeric[ii])]:
            for k,axis in enumerate('xyz'):
                arr=self.lmp.numpy.extract_atom('d_'+prefix+axis);arr[:nlocal]=0;arr[local]=values[:,k]
        self.seconds+=time.perf_counter()-started;self.calls+=1
    def save_event(self,name,step,state):
        X=self.g['X'].copy();X[:self.n]=state[:,:3];v=np.zeros_like(X);v[:self.n]=state[:,3:6]
        np.savez_compressed(self.path/(name+'_state.npz'),step=step,position_m=X,velocity_force_half_m_s=v,attached=~self.detached,component_id=self.cid[self.lab],alive=self.alive)
    def close(self):
        if self.comm.rank==0:
            self.file.close();self.fullfile.close()
            np.savez_compressed(self.path/'ONLINE_TOPOLOGY_EVENTS.npz',step_and_initial_bond_index=np.array(self.topology_events,int).reshape(-1,2),detach_step=self.detach_step,drag_step=self.drag_step)
            csv_out(self.path/'EARLY_FLOW_STEP_SUMMARY.csv',self.history)
            csv_out(OUT/'DETACHMENT_DRAG_EVENT_AUDIT.csv',[dict(atom_id=i+1,t_detach_s=int(s)*DT,detach_step=int(s),t_drag_on_s=int(self.drag_step[i])*DT if self.drag_step[i]>=0 else None,drag_on_step=int(self.drag_step[i]),delay_s=int(self.drag_step[i]-s)*DT if self.drag_step[i]>=0 else None) for i,s in enumerate(self.detach_step) if s>=0])
            dump(self.path/'EARLY_FLOW_CALLBACK.json',dict(calls=self.calls,seconds_rank0=self.seconds,audit_rows=self.audit_count,snapshot_rows=self.snapshot_count,maximum_detached_nodes=self.max_detached,maximum_detached_component_nodes=self.max_component,topology_event_count=len(self.topology_events),native_source_modified=False,no_velocity_writes=True,fluid_force_only_detached=True,numerical_damping_only_attached=True))

class AuditedWall(WallCallback):
    def compute(self,*args):
        super().compute(*args);self.last_force=args[-1].copy()

def driver():
    from mpi4py import MPI
    from lammps import lammps
    comm=MPI.COMM_WORLD;p=Path.cwd();L=lammps(comm=comm,cmdargs=['-log','log.lammps','-screen','screen.lammps']);wall=flow=None
    try:
        prefix,suffix=(p/'in.lammps').read_text().split('\nrun ',1);L.commands_string(prefix)
        wall=AuditedWall(L,comm,p,1927);flow=EarlyFlow(L,comm,p,wall)
        L.set_fix_external_callback('nativewall',wall.callback);L.set_fix_external_callback('earlyflow',flow.callback)
        L.commands_string('run '+suffix)
        # Full precision accepted terminal state, after native final_integrate.
        x=np.ctypeslib.as_array(L.gather_atoms('x',1,3)).reshape(-1,3).copy();v=np.ctypeslib.as_array(L.gather_atoms('v',1,3)).reshape(-1,3).copy()
        if comm.rank==0:np.savez_compressed(p/'terminal_full_precision.npz',position_m=x,velocity_m_s=v)
    finally:
        if flow:flow.close()
        if wall:wall.close()
        L.close()

def prepare(path,zero=False):
    path.mkdir(parents=True,exist_ok=False)
    for file in ['body.data','request.json','LOAD_PROVIDER.json']:shutil.copy2(OLD/file,path/file)
    s=(OLD/'in.lammps').read_text()
    s=re.sub(r'^fix drag clot viscous .*\n','',s,flags=re.M)
    s=s.replace('fix fixed anchors','fix flowprops all property/atom d_dragx d_dragy d_dragz d_numx d_numy d_numz\nfix earlyflow clot external pf/callback 1 1\nfix fixed anchors')
    s=s.replace('f_nativewall[3]\ndump_modify','f_nativewall[3] d_dragx d_dragy d_dragz d_numx d_numy d_numz\ndump_modify')
    # Native output precision is sufficient for strict handoff, no rounded reset.
    s=s.replace('format float %.15g','format float %.17g')
    if zero:s=s.replace('run 4660\n','run 0\n')
    (path/'in.lammps').write_text(s)

def launch(path):
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONPATH=f'{SOURCE}/python:{DEPS}:{ROOT}',LD_LIBRARY_PATH=str(LIBDIR)+os.pathsep+os.environ.get('LD_LIBRARY_PATH',''),MPI4PY_MPIABI='openmpi',MPI4PY_LIBMPI=os.environ.get('MPI4PY_LIBMPI','/usr/lib/x86_64-linux-gnu/libmpi.so.40'))
    cmd=['mpirun','--bind-to','core','-np','2',sys.executable,'-m','pd_clot.lammps_fragmentation.early_flow_coupling','driver'];start=time.perf_counter()
    dump(path/'RUN_COMMAND.json',dict(argv=cmd,cwd=path,environment={k:env[k] for k in ['PYTHONPATH','LD_LIBRARY_PATH','MPI4PY_MPIABI','MPI4PY_LIBMPI','OMP_NUM_THREADS']}))
    with (path/'launcher.log').open('w') as f:proc=subprocess.run(['/usr/bin/time','-f','%U %S %e %M %x','-o','resource.time',*cmd],cwd=path,env=env,stdout=f,stderr=subprocess.STDOUT)
    text=(path/'log.lammps').read_text() if (path/'log.lammps').exists() else '';m=re.search(r'Loop time of ([\d.e+-]+) on (\d+) procs for (\d+) steps',text)
    dump(path/'PERFORMANCE.json',dict(exit_code=proc.returncode,wall_seconds=time.perf_counter()-start,explicit_steps=int(m[3]) if m else None,loop_wall_s=float(m[1]) if m else None,mpi_ranks=2))
    assert proc.returncode==0 and 'ERROR:' not in text,(path/'launcher.log').read_text()

def preflight():
    path=OUT/'interface_check_run0';prepare(path,True);launch(path)
    q=json.loads((path/'EARLY_FLOW_CALLBACK.json').read_text());assert q['calls']==1 and q['topology_event_count']==0 and q['audit_rows']==0
    z=np.load(path/'terminal_full_precision.npz');g=load_geometry();assert np.max(abs(z['position_m']-g['X']))<1e-15 and not z['velocity_m_s'].any()
    rows=np.fromfile(path/'saved_frame_force_audit.bin',dtype=DTYPE);assert len(rows)==1927 and not rows['drag_N'].any() and not rows['numerical_N'].any()
    # Direct physical law check, no auxiliary dynamics or future topology.
    for count in [1,2,6]:
        v=np.array([[.4,.1,-.2]]);u=np.array([[.005,0,0]]);F,tau,re=drag_law(np.array([count]),v,u)
        assert np.dot(F[0],v[0]-u[0])<0 and np.allclose((F/count)*count,F)
        if count<=2:assert tau[0]==TAU
    dump(OUT/'INTERFACE_READINESS.json',dict(status='PASS',native_steps_advanced=0,mpi_ranks=2,all_initial_bonds_exact=True,all_initial_nodes_exact=True,callback_and_audit_write_verified=True,attached_flow_force_zero=True,source_SHA=subprocess.check_output(['git','rev-parse','HEAD'],cwd=SOURCE,text=True).strip(),library_sha256=digest(LIBDIR/'liblammps.so'),source_files={str(p.relative_to(SOURCE)):digest(p) for p in [SOURCE/'src/BPM/bond_bpm.cpp',SOURCE/'src/BPM/bond_bpm_peri.cpp',SOURCE/'src/library.cpp',SOURCE/'src/verlet.cpp']},native_hook_required=False))

def run():
    assert json.loads((OUT/'INTERFACE_READINESS.json').read_text())['status']=='PASS'
    assert (RESULTS/'animations/y_clot_point_cloud_only_preview.mp4').exists(),'Deliver old-data point-cloud preview first'
    assert not (OUT/'MECHANICAL_CASE_LAUNCHED.json').exists(),'Exactly one new dynamics case permitted'
    prepare(CASE);dump(OUT/'MECHANICAL_CASE_LAUNCHED.json',dict(case=CASE,epoch=time.time(),allowed_count=1,old_baseline_rerun=False));launch(CASE)
    from .y_v3_cases import analyze
    g=load_geometry();analyze(dict(path=CASE,request=json.loads((CASE/'request.json').read_text()),geometry=g,load=json.loads((CASE/'LOAD_PROVIDER.json').read_text()),steps=4660))

if __name__=='__main__':{'driver':driver,'preflight':preflight,'run':run}[sys.argv[1]]()
