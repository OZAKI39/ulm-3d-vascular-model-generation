"""One small geometry: original/fixed T0, real sort, 2-rank migration, one fracture/restart.
No integration fixes: controlled coordinate changes, native compute and true events.
"""
from pathlib import Path
import numpy as np,json,sys,ctypes as ct,time,resource
from scipy.spatial import cKDTree
from mpi4py import MPI
from lammps import lammps
C=json.loads(Path(__file__).with_name('config.json').read_text());ROOT=Path.cwd();O=ROOT/'verification'/C['task']/C['run'];mode=sys.argv[1];variant=sys.argv[2];out=O/(variant+'_'+mode);out.mkdir(exist_ok=True)
comm=MPI.COMM_WORLD;rank=comm.rank;DX=50e-6;V=DX**3;K=33333.333333333336;G=940.;DELTA=150.5e-6
lib=ct.CDLL(str(O/('probe_'+variant+'.so')));lib.observe.argtypes=[ct.c_void_p,ct.c_int];lib.observe.restype=ct.POINTER(ct.c_double);lib.grow_atom_storage.argtypes=[ct.c_void_p]
q=np.array(np.meshgrid(np.arange(7)-3,np.arange(7)-3,np.arange(7)-3,indexing='ij')).reshape(3,-1).T.astype(float);q=q[~((q[:,0]==3)&(q[:,1]==3)&(q[:,2]>1))];q[:,0]-=.1
q=q[np.random.default_rng(41011).permutation(len(q))];X=q*DX;n=len(X);b=np.array(sorted(cKDTree(X).query_pairs(DELTA)),int);a,z=b.T;r0=np.linalg.norm(X[z]-X[a],axis=1);h=np.where(abs(r0-DELTA)<=DX/2,(DELTA+DX/2-r0)/DX,1.);m=np.bincount(np.r_[a,z],weights=np.tile(h*r0*V,2),minlength=n);degree=np.bincount(b.ravel(),minlength=n)
if rank==0:
 np.savez_compressed(out/'reference.npz',id=np.arange(1,n+1),X=X,bonds=b,r0=r0,taper=h,V=np.full(n,V),m=m,degree=degree)
 s=f'Tiny nonuniform LPS identity specimen\n\n{n} atoms\n{len(b)} bonds\n\n1 atom types\n1 bond types\n\n-0.00025 0.00025 xlo xhi\n-0.00025 0.00025 ylo yhi\n-0.00025 0.00025 zlo zhi\n\nMasses\n\n1 {1000*V:.17g}\n\nAtoms # bond\n\n'
 s+=''.join(f'{i+1} 1 1 '+ ' '.join(f'{v:.17g}' for v in x)+'\n' for i,x in enumerate(X));s+='\nBonds\n\n'+''.join(f'{k+1} 1 {i+1} {j+1}\n' for k,(i,j) in enumerate(b));(out/'body.data').write_text(s)
comm.Barrier()
def independent(x,alive,theta_override=None):
 aa,zz=b[alive].T;rr=r0[alive];hh=h[alive];y=x[zz]-x[aa];length=np.linalg.norm(y,axis=1);e=length-rr;e[abs(e)<2.2204e-16]=0
 theta=3*V/m*np.bincount(np.r_[aa,zz],weights=np.tile(hh*e,2),minlength=n)
 U=.5*V*(K-5*G/3)*np.sum(theta**2)+7.5*G*V**2*np.sum(hh/rr*e**2*(1/m[aa]+1/m[zz]))
 tt=theta if theta_override is None else theta_override
 t=V**2*hh*((3*K-5*G)*(tt[aa]/m[aa]+tt[zz]/m[zz])+15*G*(1/m[aa]+1/m[zz])*e/rr)
 force=np.zeros_like(x);f=t[:,None]*y/length[:,None];np.add.at(force,aa,f);np.add.at(force,zz,-f)
 return theta,force,float(U)
def new(restart=None):
 L=lammps(comm=comm,cmdargs=['-log',str(out/('restart.log' if restart else 'log.lammps')),'-screen','none'])
 start='units si\ndimension 3\nboundary f f f\natom_style bond\natom_modify map array sort 0 0\nnewton on off\nprocessors '+ ('2 1 1' if comm.size==2 else '1 1 1')+'\n'
 if mode=='sort':start=start.replace('sort 0 0',f'sort 2 {DX:.17g}')
 text=start+(f'read_restart {restart}\n' if restart else f'read_data {out/"body.data"} extra/bond/per/atom 150 extra/special/per/atom 200\n')
 text+=f'''lattice sc {DX:.17g}
neighbor {DX:.17g} bin
neigh_modify every 1 delay 0 check no
comm_modify cutoff {3*DELTA:.17g}
fix props all property/atom d_vfrac d_theta ghost yes
'''
 if not restart:text+=f'set group all d_vfrac {V:.17g}\npair_style zero {DELTA:.17g}\npair_coeff * *\nbond_style bpm/peri\nbond_coeff 1 lps {K:.17g} {G} {DELTA} .15 0\n'
 text+='special_bonds lj 0 1 1 coul 1 1 1\ntimestep 5e-7\nthermo 1\nthermo_style custom step atoms bonds pe\nthermo_modify lost error lost/bond error\nrun 0\n'
 if rank==0:(out/('restart.in' if restart else 'in.lammps')).write_text(text)
 L.commands_string(text);return L
reports={};states={}
def snap(L,name):
 nl=L.extract_global('nlocal');ng=L.extract_global('nghost');ids=L.numpy.extract_atom('id')[:nl].copy();xx=L.numpy.extract_atom('x')[:nl].copy();ff=L.numpy.extract_atom('f')[:nl].copy()*V
 mm=np.ctypeslib.as_array(lib.observe(L.lmp,0),shape=(nl+ng,)).copy();tt=np.ctypeslib.as_array(lib.observe(L.lmp,1),shape=(nl+ng,)).copy();allid=L.numpy.extract_atom('id')[:nl+ng].copy()
 nb=L.extract_atom('num_bond');bt=L.extract_atom('bond_type');ba=L.extract_atom('bond_atom');edges=[]
 for i in range(nl):
  for j in range(nb[i]):
   if bt[i][j]>0 and ids[i]<ba[i][j]:edges.append((ids[i]-1,ba[i][j]-1))
 packets=comm.gather((ids,xx,ff,mm[:nl],tt[:nl],np.full(nl,rank),allid[nl:],mm[nl:],tt[nl:],edges),root=0)
 pe=L.get_thermo('pe')*V
 if rank==0:
  col=[np.concatenate([p[k] for p in packets]) for k in range(9)];order=np.argsort(col[0]);iid,x,f,nm,nt,owner=[c[order] for c in col[:6]];assert np.array_equal(iid,np.arange(1,n+1))
  live=set(e for p in packets for e in p[9]);alive=np.array([tuple(e) in live for e in b]);theta,ref,U=independent(x,alive);errm=nm-m;errf=f-ref;errt=nt-theta
  fdnode=0;eps=DX*1e-5;xplus=x.copy();xminus=x.copy();xplus[fdnode,0]+=eps;xminus[fdnode,0]-=eps;fd=-(independent(xplus,alive)[2]-independent(xminus,alive)[2])/(2*eps)
  scaleF=float(np.max(np.linalg.norm(ref,axis=1)));scaleT=float(np.max(abs(theta)));scaleM=float(m.max());absolute={'m':float(abs(errm).max()),'theta':float(abs(errt).max()),'force_N':float(np.linalg.norm(errf,axis=1).max())}
  ghost_id,ghost_m,ghost_theta=col[6:9];ge=float(np.max(abs(ghost_m-m[ghost_id-1]))) if len(ghost_id) else None
  report={'absolute_errors':absolute,'scales':{'m_m4':scaleM,'theta':scaleT,'force_N':scaleF},'relative_errors':{'m':absolute['m']/scaleM,'theta':absolute['theta']/scaleT if scaleT else None,'force':absolute['force_N']/scaleF if scaleF else None},'wrong_m_node_count':int((abs(errm)>2e-12*scaleM).sum()),'wrong_m_atom_ids':iid[abs(errm)>2e-12*scaleM].tolist(),'ghost_m_max_abs_error':ge,'ghost_count':len(ghost_id),'active_bonds':int(alive.sum()),'independent_energy_J':U,'native_PE_times_V_J':pe,'independent_FD_force_N':fd,'analytical_Fx_atom1_N':float(ref[0,0]),'force_sum_N':f.sum(0).tolist(),'torque_Nm':np.cross(x-x.mean(0),f).sum(0).tolist(),'local_id_order':[p[0].tolist() for p in packets],'owner_rank_by_id':owner.tolist()}
  reports[name]=report;states[name]=(x,f,nm,nt,alive)
  np.savez_compressed(out/(name+'.npz'),id=iid,x=x,force_N=f,m=nm,theta=nt,owner=owner,alive=alive,reference_force_N=ref,reference_theta=theta,ghost_id=ghost_id,ghost_m=ghost_m,ghost_theta=ghost_theta)
  (out/'report.json').write_text(json.dumps(reports,indent=2,allow_nan=False));print(variant,mode,name,report['relative_errors'],report['wrong_m_node_count'],flush=True)
def coords(L,x):
 nl=L.extract_global('nlocal');ids=L.numpy.extract_atom('id')[:nl];L.numpy.extract_atom('x')[:nl]=x[ids-1]
t0=time.perf_counter();L=new();snap(L,'initial')
D=np.array([[1.012,.004,0],[0,.997,.002],[0,0,1.002]]);x=X@D.T;coords(L,x);L.command('run 1 pre no post no');snap(L,'T0')
if mode=='sort':
 x=x@np.array([[0,-1,0],[1,0,0],[0,0,1]]).T;coords(L,x);L.command('run 1 pre no post no');snap(L,'T1_sort')
elif mode=='migration':
 x[:,0]+=.2*DX;coords(L,x);L.command('run 1 pre no post no');snap(L,'T2_migration')
elif mode in ['fracture','mpi_fracture']:
 corner=int(np.argmax(X.sum(1)));axis=0
 if mode=='mpi_fracture':
  corner=int(np.argmin(np.linalg.norm(q-np.array([.9,3,3]),axis=1)));axis=1
 x[corner,axis]+=.3*DX;coords(L,x);L.command('run 1 pre no post no');snap(L,'break')
 L.command('run 1 pre no post no');snap(L,'postbreak_same_coordinates')
 L.command('write_restart '+str(out/'tiny.restart'));lib.grow_atom_storage(L.lmp);L.command('run 1 pre no post no');snap(L,'growth_after_break')
 L.command('run 0');snap(L,'new_run_after_break');L.close();L=new(out/'tiny.restart');snap(L,'restart_after_break')
L.close()
if rank==0:
 (out/'PERFORMANCE.json').write_text(json.dumps({'wall_s':time.perf_counter()-t0,'cpu_user_s':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'max_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'mpi_ranks':comm.size,'coordinate_driven_no_integrator':True},indent=2))
