"""Experimental optional Warp mesh queries. Never imported by default production."""
from time import perf_counter
import numpy as np


def mesh_audit(surface):
    faces=surface.faces.reshape(-1,4)[:,1:]
    points=np.asarray(surface.points,dtype=float)
    tri=points[faces]
    edge=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]])
    ordered=np.sort(edge,axis=1)
    _,inverse,counts=np.unique(ordered,axis=0,return_inverse=True,return_counts=True)
    orientation=np.bincount(inverse,weights=np.where(edge[:,0]<edge[:,1],1,-1))
    twice_area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    volume=float(np.einsum("ij,ij->i",tri[:,0],np.cross(tri[:,1],tri[:,2])).sum()/6)
    return dict(point_count=len(points),triangle_count=len(faces),boundary_edges=int((counts==1).sum()),
        nonmanifold_edges=int((counts>2).sum()),inconsistent_edge_winding=int((orientation!=0).sum()),
        zero_area_triangles=int((twice_area==0).sum()),signed_volume_m3=volume,
        closed_consistent_outward=bool(np.all(counts==2) and np.all(orientation==0) and volume>0 and np.all(twice_area>0)),
        limitation="Topological/winding and nondegeneracy audit; no new proof excluding every geometric self-intersection")


class WarpWallBackend:
    """Stable Warp mesh BVH and official sign query, with measured host transfers.

    GPU query accepts a contiguous CuPy float32 array via cuda-array interface.
    The current CPU tracker requires a host adapter; those transfers are recorded,
    not hidden in a GPU-only timing. Float32 mesh precision is tested without
    loosening the original subnanometre geometric/contact thresholds.
    """
    def __init__(self,surface,device="cuda:0"):
        import cupy as cp
        import warp as wp
        self.cp,self.wp,self.device=cp,wp,device
        begin=perf_counter()
        # Match the original VTK's exact duplicate-point cleaning/triangulation.
        self.surface=surface.clean(tolerance=0.).triangulate()
        self.source_audit=mesh_audit(self.surface)
        # Audit the coordinates actually uploaded, including float32 quantization.
        query_surface=self.surface.copy(deep=True)
        query_surface.points=np.asarray(self.surface.points,np.float32)
        self.audit=mesh_audit(query_surface)
        if not self.source_audit["closed_consistent_outward"] or not self.audit["closed_consistent_outward"]:
            raise ValueError("Warp signed wall query requires audited closed outward mesh")
        wp.init()
        self.points=wp.array(np.asarray(self.surface.points,np.float32),dtype=wp.vec3,device=device)
        self.indices=wp.array(self.surface.faces.reshape(-1,4)[:,1:].astype(np.int32).ravel(),dtype=wp.int32,device=device)
        self.mesh=wp.Mesh(points=self.points,indices=self.indices)
        self.build_seconds=perf_counter()-begin
        self.host_to_device_bytes=self.device_to_host_bytes=0
        self.calls=self.points_queried=0
        self.zero_copy_buffer_alias_checks=0
        self.launch_and_synchronize_seconds=self.transfer_seconds=0.

    def query_gpu(self,points):
        from mevo_mb.coupling.warp_wall_kernel import query_mesh
        cp,wp=self.cp,self.wp
        if points.dtype!=cp.float32 or not points.flags.c_contiguous:
            raise ValueError("Zero-copy Warp input requires contiguous CuPy float32 Nx3")
        n=len(points)
        distance=cp.empty(n,cp.float32)
        closest=cp.empty((n,3),cp.float32)
        normal=cp.empty((n,3),cp.float32)
        face=cp.empty(n,cp.int32)
        buffers=[(points,wp.vec3),(distance,wp.float32),(closest,wp.vec3),(normal,wp.vec3),(face,wp.int32)]
        arrays=[wp.array(array,dtype=dtype,copy=False,device=self.device) for array,dtype in buffers]
        assert all(wp_array.ptr==cp_array.data.ptr for wp_array,(cp_array,_) in zip(arrays,buffers))
        self.zero_copy_buffer_alias_checks+=len(arrays)
        # Synchronize the CuPy producer before consuming its borrowed allocation.
        cp.cuda.get_current_stream().synchronize()
        begin=perf_counter()
        wp.launch(query_mesh,dim=n,inputs=[self.mesh.id,arrays[0]],outputs=arrays[1:],device=self.device)
        wp.synchronize_device(self.device)
        # Host launch + synchronization includes the first JIT compilation;
        # this is deliberately not labelled device-only kernel time.
        self.launch_and_synchronize_seconds+=perf_counter()-begin
        self.calls+=1
        self.points_queried+=n
        return distance,closest,normal,face

    def query(self,points):
        cp=self.cp
        x=np.asarray(points,dtype=np.float32,order="C")
        begin=perf_counter()
        gpu=cp.asarray(x)
        cp.cuda.get_current_stream().synchronize()
        self.transfer_seconds+=perf_counter()-begin
        self.host_to_device_bytes+=x.nbytes
        values=self.query_gpu(gpu)
        begin=perf_counter()
        result=tuple(cp.asnumpy(value) for value in values)
        self.transfer_seconds+=perf_counter()-begin
        self.device_to_host_bytes+=sum(v.nbytes for v in result)
        if np.any(result[3]<0) or not all(np.isfinite(v).all() for v in result[:3]):
            raise RuntimeError("Warp closest mesh query failed; no substitute distance accepted")
        return tuple(v.astype(np.float64) for v in result[:3])+(result[3],)
