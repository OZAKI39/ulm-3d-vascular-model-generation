"""Deterministic >=10k wall-query and projection fixtures, including gold contacts."""
import argparse
import json
import time

import h5py
import numpy as np
import pyvista as pv
import vtk

from mevo_mb.coupling.bifurcation import BifurcationGeometry
from mevo_mb.coupling.wall_backends import PersistentVTKWallBackend
from scripts.stage4i5.common import BASE,GOLD,new_run,finish,write,sha


def wall_input():
    config=json.loads((GOLD/"cases/B2/RUN_IDENTITY.json").read_text())["config"]["geometry"]
    path=config["mesh_surface"].split("/workspace/mevo_acoustic_microbubbles/")[-1]
    return config,pv.read(GOLD/"inputs"/path)


def fixtures(surface):
    rng=np.random.default_rng(20261005)
    groups=[]
    for label in ("B2","B3"):
        with h5py.File(GOLD/"cases"/label/"particles/trajectories.h5") as h:
            x=h["positions_m"][:][h["active"][:]]
            groups.append(("gold_"+label,x[rng.choice(len(x),1000,replace=False)]))
    faces=surface.faces.reshape(-1,4)[:,1:]
    vertices=np.asarray(surface.points,dtype=float)
    tri=vertices[faces]
    normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    normal/=np.linalg.norm(normal,axis=1)[:,None]
    centers=tri.mean(axis=1)
    idx=rng.integers(len(tri),size=2000)
    # Around the actual projection threshold; no modification of solver tolerance.
    offsets=-1.78e-6+rng.choice([-2e-12,-1e-12,0.,1e-12,2e-12],2000)
    groups.append(("projection_threshold",centers[idx]+offsets[:,None]*normal[idx]))
    radius=1.3e-3*np.sqrt(rng.random(1500))*.8
    angle=rng.uniform(0,2*np.pi,1500)
    groups.append(("parent_lumen",np.column_stack([rng.uniform(-.012,-.002,1500),radius*np.cos(angle),radius*np.sin(angle)])))
    daughter=[]
    for sign in (-1,1):
        d=np.array([np.sqrt(3)/2,sign*.5,0.])
        tangent=np.array([-sign*.5,np.sqrt(3)/2,0.])
        r=.8e-3*np.sqrt(rng.random(750))
        a=rng.uniform(0,2*np.pi,750)
        daughter.append(rng.uniform(.003,.0105,750)[:,None]*d+r[:,None]*np.cos(a)[:,None]*tangent+
                        (r*np.sin(a))[:,None]*np.array([0.,0.,1.]))
    groups.append(("daughter_lumen",np.concatenate(daughter)))
    groups.append(("bifurcation_apex",rng.uniform([-.0015,-.001,-.0007],[.0015,.001,.0007],(1000,3))))
    idx=rng.integers(len(tri),size=1000)
    groups.append(("slightly_inside_outside",centers[idx]+rng.uniform(-5e-6,5e-6,1000)[:,None]*normal[idx]))
    idx=rng.integers(len(tri),size=1000)
    groups.append(("edges_vertices",(tri[idx,0]+tri[idx,1])/2+rng.normal(size=(1000,3))*1e-8))
    return np.concatenate([x for _,x in groups]),np.concatenate([np.full(len(x),label,dtype="S32") for label,x in groups])


def reference_query(implicit,points):
    distance=np.empty(len(points))
    closest=np.empty_like(points)
    normals=np.empty_like(points)
    cp=[0.,0.,0.]
    grad=[0.,0.,0.]
    for i,x in enumerate(points):
        distance[i]=implicit.EvaluateFunctionAndGetClosestPoint(x,cp)
        implicit.EvaluateGradient(x,grad)
        closest[i],normals[i]=cp,grad
    return distance,closest,normals


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--backend",choices=["vtk","warp"],default="vtk")
    args=parser.parse_args()
    out=new_run("p1_vtk_contact" if args.backend=="vtk" else "p2_warp",dict(test="wall_fixture",backend=args.backend,points=10000,seed=20261005))
    config,surface=wall_input()
    from mevo_mb.coupling.warp_wall import mesh_audit
    audit=mesh_audit(surface.clean(tolerance=0.).triangulate())
    write(out/"MESH_WATERTIGHTNESS_AUDIT.json",audit)
    fixture=BASE/"protocol"/"WALL_FIXTURES.npz"
    # Sealed protocol remains untouched: fixture is versioned with this benchmark.
    if fixture.exists():
        data=np.load(fixture)
        points,labels=data["points"],data["labels"]
    else:
        points,labels=fixtures(surface)
    np.savez_compressed(out/"fixtures.npz",points=points,labels=labels)
    implicit=vtk.vtkImplicitPolyDataDistance()
    implicit.SetInput(surface)
    if args.backend=="vtk":
        backend=PersistentVTKWallBackend(surface)
    else:
        from mevo_mb.coupling.warp_wall import WarpWallBackend
        backend=WarpWallBackend(surface)
        write(out/"WARP_MESH_REPRESENTATION_AUDIT.json",dict(source_original=backend.source_audit,
            uploaded_float32=backend.audit,source_coordinate_dtype=str(backend.surface.points.dtype),
            uploaded_coordinate_dtype="float32"))
    # Untimed warmup; initialization and kernel compilation excluded explicitly.
    reference_query(implicit,points[:4])
    backend.query(points[:4])
    begin=time.perf_counter()
    ref=reference_query(implicit,points)
    ref_s=time.perf_counter()-begin
    begin=time.perf_counter()
    result=backend.query(points)
    candidate_s=time.perf_counter()-begin
    distance,closest,normals,_=result
    normal_cos=np.sum(normals*ref[2],axis=1)/(np.linalg.norm(normals,axis=1)*np.linalg.norm(ref[2],axis=1))
    angles=np.rad2deg(np.arccos(np.clip(normal_cos,-1,1)))
    ref_excess=ref[0]+1.78e-6*np.linalg.norm(ref[2],axis=1)
    candidate_excess=distance+1.78e-6*np.linalg.norm(normals,axis=1)
    projection_mismatch=(ref_excess>1e-12)!=(candidate_excess>1e-12)
    constraint_mismatch=(ref_excess>=-1e-12)!=(candidate_excess>=-1e-12)
    groups=[]
    for label in np.unique(labels):
        mask=labels==label
        groups.append(dict(group=label.decode(),count=int(mask.sum()),
            max_distance_error_m=float(np.max(np.abs(distance[mask]-ref[0][mask]))),
            max_closest_point_error_m=float(np.max(np.linalg.norm(closest[mask]-ref[1][mask],axis=1))),
            max_normal_angle_deg=float(angles[mask].max()),
            inside_outside_mismatches=int(((distance[mask]<0)!=(ref[0][mask]<0)).sum()),
            projection_decision_mismatches=int(projection_mismatch[mask].sum()),
            velocity_constraint_decision_mismatches=int(constraint_mismatch[mask].sum())))
    original=BifurcationGeometry(config)
    original.bind_mesh_surface(surface)
    candidate=BifurcationGeometry(config)
    candidate.bind_mesh_surface(surface)

    def level(x):
        d,_,n,_=backend.query(x)
        return d,n

    candidate.level_and_gradient=level
    projection={}
    try:
        a,hit_a=original.project(points,1.78e-6)
        b,hit_b=candidate.project(points,1.78e-6)
        error=np.linalg.norm(a-b,axis=1)
        remaining=original.wall_distance(b)
        projection=dict(max_position_error_m=float(error.max()),rms_position_error_m=float(np.sqrt(np.mean(error**2))),
            hit_mismatches=int((hit_a!=hit_b).sum()),maximum_penetration_m=float(np.maximum(1.78e-6-remaining,0).max()),
            projection_failures=0)
    except RuntimeError as exc:
        projection=dict(projection_failures=1,error=repr(exc))
    np.savez_compressed(out/"query_results.npz",reference_distance=ref[0],reference_closest=ref[1],reference_normals=ref[2],
        distance=distance,closest=closest,normals=normals,cell_ids=result[3],
        projection_mismatch=projection_mismatch,velocity_constraint_mismatch=constraint_mismatch)
    write(out/"BOUNDARY_DIFFERENCES.json",[dict(index=int(i),group=labels[i].decode(),point=points[i].tolist(),
        reference_excess_m=float(ref_excess[i]),candidate_excess_m=float(candidate_excess[i]))
        for i in np.flatnonzero(projection_mismatch|constraint_mismatch)])
    backend_details={key:getattr(backend,key) for key in
        ("calls","points_queried","launch_and_synchronize_seconds","transfer_seconds","host_to_device_bytes","device_to_host_bytes","zero_copy_buffer_alias_checks") if hasattr(backend,key)}
    backend_details["timing_scope"]="Host launch/synchronize and transfer aggregates include fixture warmup and projection queries; timed query excludes warmup. No CUDA-event device-only time claimed."
    finish(out,dict(status="BENCHMARK_COMPLETE",backend=args.backend,point_count=len(points),groups=groups,
        reference_query_s=ref_s,candidate_query_s=candidate_s,speedup=ref_s/candidate_s,
        build_s=backend.build_seconds,projection=projection,fixture_sha256=sha(out/"fixtures.npz"),
        mesh_audit=backend.audit if args.backend=="warp" else audit,source_mesh_audit=audit,backend_details=backend_details,
        adoption="PENDING_STRICT_SHORT_TRAJECTORY_GATE",VTK_version=vtk.vtkVersion.GetVTKVersion(),
        source="https://vtk.org/doc/nightly/html/classvtkStaticCellLocator.html",
        original_audit="Original vtkImplicitPolyDataDistance already builds vtkCellLocator once and reuses thread-local generic cells. Separate value and gradient calls repeat closest searches; no per-query locator rebuild existed."))


if __name__=="__main__":
    main()
