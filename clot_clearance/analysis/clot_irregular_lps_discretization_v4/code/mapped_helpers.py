"""Actual mapped reference-cell quadrature and conservative visibility helpers."""
from common import *
from itertools import product
import ctypes,subprocess
from scipy.sparse import csr_matrix,save_npz,load_npz
from scipy.spatial import cKDTree

def cells(ref):
 bits=np.array(list(product([0,1],repeat=3)));n=int(ref['clot_count']);N=len(ref['X'])
 corners=np.empty((N,8,3));corners[:n]=ref['control_corners_m'][ref['control_corner_ids']]
 corners[n:]=ref['X'][n:,None]+(bits-.5)*DX
 faces=[]
 for axis in range(3):
  others=[k for k in range(3) if k!=axis]
  for side in range(2):
   ids=[]
   for a,b in [(0,0),(1,0),(1,1),(0,1)]:
    key=np.zeros(3,int);key[axis]=side;key[others]=[a,b];ids.append(int(np.flatnonzero((bits==key).all(1))[0]))
   q=bits[ids];normal=np.cross(q[1]-q[0],q[2]-q[0]);
   if normal[axis]*(2*side-1)<0:ids=[ids[k] for k in [0,3,2,1]]
   faces.append(ids)
 return corners,np.array(faces,int)

def initialize(tool):
 corners,faceids=cells(tool.ref);tool.cell_corner_offset=corners-tool.ref['X'][:,None];tool.cell_reference_radius=np.linalg.norm(tool.cell_corner_offset,axis=2).max(1)
 tool.face_quads=np.ascontiguousarray(corners[:tool.n,faceids].reshape(-1,4,3));q=tool.face_quads
 c1=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);c2=np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0]);a1=np.linalg.norm(c1,axis=1);a2=np.linalg.norm(c2,axis=1)
 tool.offset=((q[:,0]+q[:,1]+q[:,2])*a1[:,None]+(q[:,0]+q[:,2]+q[:,3])*a2[:,None])/(3*(a1+a2)[:,None])-tool.ref['X'][tool.owner]
 tool.reference_area_vectors=(c1+c2)/2;area=np.linalg.norm(tool.reference_area_vectors,axis=1);normal=tool.reference_area_vectors/area[:,None]
 tool.reference_normal=np.vstack([tool.base.reference_normal,normal]);tool.reference_area=np.r_[tool.base.reference_area,area]
 tool.face_corner_ids=faceids;tool.cell_pattern=tool.ref['control_tet_pattern']

def load_kernel(tool):
 target=O/'build/libmapped_geometry.so';source=A/'mapped_geometry.cpp';stamp=target.with_suffix('.source.sha256')
 if not target.exists() or not stamp.exists() or stamp.read_text().strip()!=sha(source):
  subprocess.run(['g++','-O3','-std=c++17','-fPIC','-shared',str(source),'-o',str(target)],check=True,capture_output=True);stamp.write_text(sha(source)+'\n')
 tool.raylib=ctypes.CDLL(str(target));lib=tool.raylib
 lib.ray_mapped_cells.argtypes=[ctypes.c_int,ctypes.c_int]+[ctypes.c_void_p]*8+[ctypes.c_double,ctypes.c_void_p,ctypes.c_void_p];lib.ray_mapped_cells.restype=None
 lib.ray_triangles.argtypes=[ctypes.c_int,ctypes.c_int]+[ctypes.c_void_p]*3+[ctypes.c_double,ctypes.c_void_p];lib.ray_triangles.restype=None
 lib.crossings_build.argtypes=[ctypes.c_int]+[ctypes.c_void_p]*5;lib.crossings_build.restype=ctypes.c_int64
 lib.crossings_copy.argtypes=[ctypes.c_void_p,ctypes.c_void_p];lib.crossings_copy.restype=None

def reference_crossings(tool):
 cache=O/'provenance/mapped_dual_face_crossings.npz';stamp=cache.with_suffix('.reference.sha256');fingerprint=sha(O/'provenance/reference.npz')+sha(A/'mapped_geometry.cpp')
 if cache.exists():
  assert stamp.read_text().strip()==fingerprint
  return load_npz(cache)
 load_kernel(tool);X=np.ascontiguousarray(tool.ref['X']);b=np.ascontiguousarray(tool.ref['bonds'],dtype=np.int64)
 mid=X[b].mean(1);radius=tool.ref['r0']/2+tool.cell_reference_radius[:tool.n].max()+1e-12
 nearby=cKDTree(X[:tool.n]).query_ball_point(mid,radius);ptr=np.r_[0,np.cumsum([len(v) for v in nearby])].astype(np.int64);cand=np.concatenate(nearby).astype(np.int64)
 length=tool.raylib.crossings_build(len(b),*[v.ctypes.data for v in [X,b,ptr,cand,tool.face_quads]])
 rows=np.empty(length,np.int64);cols=rows.copy();tool.raylib.crossings_copy(rows.ctypes.data,cols.ctypes.data)
 out=csr_matrix((np.ones(length,np.uint8),(rows,cols)),shape=(tool.M,len(b)));out.sum_duplicates();out.data[:]=1
 save_npz(cache,out);stamp.write_text(fingerprint+'\n');return out

def mapped_ray(tool,ray_ids,center,normal,owners,x,deformation,valid,radii):
 corners=x[:,None]+np.einsum('nij,nkj->nki',deformation,tool.cell_corner_offset)
 t=corners[:,tool.cell_pattern];orig=np.ascontiguousarray(t[:,:,0]);matrix=np.transpose(t[:,:,1:]-t[:,:,:1],(0,1,3,2));inv=np.zeros_like(matrix)
 inv[valid]=np.linalg.inv(matrix[valid]);inv=np.ascontiguousarray(inv)
 inputs=[np.ascontiguousarray(center[ray_ids]),np.ascontiguousarray(normal[ray_ids]),np.ascontiguousarray(owners[ray_ids],np.int64),np.ascontiguousarray(x),orig,inv,np.ascontiguousarray(radii),np.ascontiguousarray(valid,np.uint8)]
 cl=np.zeros(len(ray_ids),np.uint8);bl=cl.copy();tool.raylib.ray_mapped_cells(len(ray_ids),tool.N,*[v.ctypes.data for v in inputs],3e-6,cl.ctypes.data,bl.ctypes.data)
 return inputs,cl,bl
