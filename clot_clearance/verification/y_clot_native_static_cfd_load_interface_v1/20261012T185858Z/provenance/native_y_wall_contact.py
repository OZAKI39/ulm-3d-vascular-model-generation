"""Native open Y sidewall penalty, called inside LAMMPS post_force on each MPI rank."""
import csv,json,time,traceback
import numpy as np
from .native_y_wall_query import NativeYWallQuery

CONTACT=dict(model_id='NATIVE_Y_WALL_REPULSION_V1',effective_radius_m=20e-6,
             stiffness_N_m=50.,damping_ratio=.5,center_penetration_tolerance_m=2e-6,
             dt_s=.5e-6,node_mass_kg=1000*(50e-6)**3)
CONTACT['damping_N_s_m']=2*CONTACT['damping_ratio']*np.sqrt(CONTACT['stiffness_N_m']*CONTACT['node_mass_kg'])
CONTACT['dt_omega']=CONTACT['dt_s']*np.sqrt(CONTACT['stiffness_N_m']/CONTACT['node_mass_kg'])

def penalty(clearance,normal,velocity):
    overlap=np.maximum(CONTACT['effective_radius_m']-clearance,0)
    magnitude=np.where(overlap>0,CONTACT['stiffness_N_m']*overlap+CONTACT['damping_N_s_m']*np.maximum(-np.einsum('ij,ij->i',velocity,normal),0),0)
    return magnitude[:,None]*normal

class WallCallback:
    def __init__(self,lmp,comm,path,n):
        self.lmp=lmp;self.comm=comm;self.path=path;self.n=n;self.query=NativeYWallQuery()
        self.previous=np.full((n,3),np.nan);self.clearance=np.full(n,np.inf);self.normal=np.zeros((n,3));self.triangle=np.full(n,-1)
        self.exited=np.zeros(n,bool);self.min_clear=np.full(n,np.inf);self.out_count=np.zeros(n,int);self.max_force=np.zeros(n)
        self.events=[];self.calls=0;self.seconds=0
        self.file=(path/f'wall_rank{comm.rank}.csv').open('w');self.writer=csv.writer(self.file)
        self.writer.writerow(['step','rank','mobile_nodes','min_clearance_m','max_center_penetration_m','contact_nodes','outside_nodes','sum_normal_force_N','sum_fx_N','sum_fy_N','sum_fz_N','query_wall_seconds'])
    def callback(self,caller,step,nlocal,tags,x,f):
        try:self.compute(step,nlocal,tags,x,f)
        except BaseException:
            (self.path/f'CALLBACK_ERROR_rank{self.comm.rank}.txt').write_text(traceback.format_exc());self.comm.Abort(31)
    def compute(self,step,nlocal,tags,x,f):
        started=time.perf_counter();f[:]=0;local=np.flatnonzero(tags<=self.n);ids=tags[local].astype(int)-1;X=x[local]
        changed=np.any(X!=self.previous[ids],axis=1)&~self.exited[ids]
        jj=ids[changed];xx=X[changed];prev=self.previous[jj].copy();old_clear=self.clearance[jj].copy()
        if len(jj):
            q=self.query.query(xx);self.clearance[jj]=q['clearance_m'];self.normal[jj]=q['inward_normal'];self.triangle[jj]=q['wall_triangle_id']
            for k in range(len(jj)):
                if not np.isfinite(prev[k]).all():continue
                port=None
                if not q['inside_lumen'][k]:port=self.query.port_crossing(prev[k],xx[k])
                if port:
                    self.exited[jj[k]]=True;self.events.append(dict(step=int(step),atom_id=int(jj[k]+1),event='OPEN_PORT_EXIT',port=port,xyz_m=xx[k].tolist()))
                elif abs(old_clear[k])<=np.linalg.norm(xx[k]-prev[k])+1e-12:
                    hit=self.query.swept_wall_crossing(prev[k],xx[k])
                    if hit[0]:self.events.append(dict(step=int(step),atom_id=int(jj[k]+1),event='SWEPT_WALL_CROSSING',xyz_m=xx[k].tolist()))
        self.previous[ids]=X
        vel=self.lmp.numpy.extract_atom('v')[:nlocal][local]
        F=penalty(self.clearance[ids],self.normal[ids],vel);F[self.exited[ids]]=0;f[local]=F
        self.min_clear[ids]=np.minimum(self.min_clear[ids],self.clearance[ids]);outside=(self.clearance[ids]<0)&~self.exited[ids]
        self.out_count[ids]+=outside.astype(int);self.max_force[ids]=np.maximum(self.max_force[ids],np.linalg.norm(F,axis=1))
        for name,values in [('wclear',self.clearance[ids]),('wpen',np.maximum(-self.clearance[ids],0)),('wtri',self.triangle[ids])]:
            arr=self.lmp.numpy.extract_atom('d_'+name);arr[:nlocal]=0;arr[local]=values
        elapsed=time.perf_counter()-started;self.seconds+=elapsed;self.calls+=1
        self.writer.writerow([int(step),self.comm.rank,len(ids),float(self.clearance[ids].min(initial=np.inf)),float(np.maximum(-self.clearance[ids][~self.exited[ids]],0).max(initial=0)),int(np.count_nonzero(np.linalg.norm(F,axis=1)>0)),int(outside.sum()),float(np.linalg.norm(F,axis=1).sum()),*F.sum(0).tolist(),elapsed])
        if step%100==0:self.file.flush()
    def close(self):
        self.file.close()
        np.savez_compressed(self.path/f'wall_node_audit_rank{self.comm.rank}.npz',minimum_clearance_m=self.min_clear,outside_step_count=self.out_count,maximum_force_N=self.max_force)
        from .config import dump
        dump(self.path/f'wall_callback_rank{self.comm.rank}.json',dict(rank=self.comm.rank,calls=self.calls,callback_seconds=self.seconds,locator_query_seconds=self.query.query_seconds,queried_points=self.query.query_points,events=self.events,model=CONTACT))

def minimal_checks():
    q=NativeYWallQuery();i=len(q.tri)//2;p=q.tri[i].mean(0);normal=q.normals[i]
    x=np.array([p+normal*10e-6,p+normal*40e-6,p-normal*1e-6]);z=q.query(x);F=penalty(z['clearance_m'],z['inward_normal'],np.zeros_like(x))
    assert np.dot(F[0],normal)>0 and np.linalg.norm(F[1])==0 and np.dot(F[2],normal)>0
    checks=q.minimal_checks();checks.update(contact_force_probe_N=F.tolist(),contact=CONTACT)
    return checks
