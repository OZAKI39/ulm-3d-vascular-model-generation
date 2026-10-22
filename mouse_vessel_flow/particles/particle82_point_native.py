"""Compiled diagnostic tracer with explicit tetra adjacency and original P1 data."""
from pathlib import Path
import ctypes,hashlib,subprocess,tempfile
import numpy as np


class NativePointTracer:
    def __init__(self,env):
        self.env=env;g=env.field.geometry;tet=g.tetra;n=len(g.points)
        faces=np.sort(np.stack([np.delete(tet,i,axis=1) for i in range(4)],axis=1),axis=2).reshape(-1,3)
        keys=(faces[:,0]*n+faces[:,1])*n+faces[:,2];order=np.argsort(keys)
        equal=np.flatnonzero(np.diff(keys[order])==0)
        nb=np.full(len(keys),-1,np.int64)
        nb[order[equal]]=order[equal+1]//4;nb[order[equal+1]]=order[equal]//4
        role=np.zeros(len(keys),np.int32)
        sorted_keys=keys[order]
        for name,s in env.boundaries.items():
            if not name.startswith('OUTLET_'):continue
            local=s.faces.reshape(-1,4)[:,1:];ids=np.sort(np.asarray(s.point_data['GlobalNodeID'],np.int64)[local]-1,axis=1)
            sk=(ids[:,0]*n+ids[:,1])*n+ids[:,2];idx=np.searchsorted(sorted_keys,sk)
            assert np.array_equal(sorted_keys[idx],sk)
            role[order[idx]]=int(name[-2:])
        self.arrays=[np.ascontiguousarray(a) for a in [g.origins,g.inverse,g.weight_tolerance,env.field.velocity_nodes_m_s,g.tetra,nb.reshape(-1,4),role.reshape(-1,4)]]
        source=Path(__file__).with_suffix('.cpp');key=hashlib.sha256(source.read_bytes()).hexdigest()[:20]
        library=Path(tempfile.gettempdir())/f'particle82_point_{key}.so'
        if not library.exists():
            temp=library.with_suffix('.tmp.so')
            subprocess.run(['g++','-std=c++17','-O3','-fno-fast-math','-ffp-contract=off','-shared','-fPIC',str(source),'-o',str(temp)],check=True)
            temp.replace(library)
        self.lib=ctypes.CDLL(str(library));P=ctypes.c_void_p;I=ctypes.c_int;L=ctypes.c_int64;D=ctypes.c_double
        self.lib.p82_init.argtypes=[P]*7+[L]
        self.lib.p82_sample.argtypes=[P,L,P];self.lib.p82_sample.restype=I
        self.lib.p82_trace.argtypes=[P,L,D,D,D,D,I,P,P,P];self.lib.p82_trace.restype=I
        self.lib.p82_init(*[a.ctypes.data for a in self.arrays],len(tet))
        self.source_sha256=hashlib.sha256(source.read_bytes()).hexdigest()

    def trace(self,position,*,step_m=.2e-6,error=1e-11,horizon_m=2e-3):
        p=np.ascontiguousarray(position,np.float64);cell,_=self.env.field.locate(p)
        if cell<0:return dict(outlet=None,end_reason='INITIAL_LOCATOR_OUTSIDE',path=np.array([[0.,*p]]),vtk_reason=None)
        out=np.empty((100000,4));count=ctypes.c_int();outlet=ctypes.c_int()
        status=self.lib.p82_trace(p.ctypes.data,cell,error,step_m,30.,horizon_m,len(out),out.ctypes.data,ctypes.byref(count),ctypes.byref(outlet))
        path=out[:count.value].copy();role=f'OUTLET_{outlet.value:02d}' if status==1 else None
        reason={1:role,2:'POINT_TIME_HORIZON_30S',3:'POINT_STEP_BUDGET',4:'POINT_DT_FLOOR',5:'POINT_LOCATOR_FAILURE',6:'POINT_PATH_LENGTH_BUDGET'}[status]
        if role:
            hit=self.env.classifier.first_event(path[-2,1:],path[-1,1:])
            if hit is None or hit.role!=role:raise ValueError('Native diagnostic cap crossing disagrees with original triangle classifier')
        return dict(outlet=role,end_reason=reason,path=path,vtk_reason=None,integrator='P1_RK23_DOUBLE_TETRA_ADJACENCY',
                    terminal_rule='LOCAL_EULER_ODE_CAP_EVENT_WITHIN_2_PERCENT_SPATIAL_STEP_NOT_ASSIGNED_OUTLET')
