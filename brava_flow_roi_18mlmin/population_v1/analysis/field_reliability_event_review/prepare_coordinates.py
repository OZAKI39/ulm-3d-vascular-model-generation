"""Geometry-only all-segment clipping and fixed common quadrature coordinates."""
from common import *
import gc
import faulthandler
import pyvista as pv
import vtk
from types import SimpleNamespace
from particle_3d.geometry import TetraGeometry
from particle_3d.population_v1.bvh import tetra_adjacency
from diagnostic_core import segment_cells
from bounds_index import BoundsIndex

def expand(seed,adj,layers):
    ids=np.unique(seed)
    for _ in range(layers):
        ids=np.unique(np.r_[ids,adj[ids].ravel()]);ids=ids[ids>=0]
    return ids

def cv_surface(points,tetra,ids,wallkeys):
    faces={}
    for cid in ids:
        cell=tetra[cid]
        for opposite in range(4):
            face=np.delete(cell,opposite);key=tuple(sorted(face))
            if key in faces:del faces[key]
            else: faces[key]=(face,int(cid),int(cell[opposite]))
    tri=[];normal=[];area=[];wall=[];nodes=[];parents=[]
    for key,(face,cid,opposite) in sorted(faces.items()):
        xyz=points[face];cross=np.cross(xyz[1]-xyz[0],xyz[2]-xyz[0])
        if np.dot(cross,points[opposite]-xyz[0])>0:
            face=face[[0,2,1]];xyz=points[face];cross=-cross
        a=np.linalg.norm(cross)/2
        tri.append(xyz);normal.append(cross/(2*a));area.append(a);wall.append(key in wallkeys);nodes.append(face);parents.append(cid)
    return dict(triangles_m=np.array(tri),normal_out=np.array(normal),area_m2=np.array(area),
                original_node_ids=np.array(nodes),parent_cells=np.array(parents),wall_face=np.array(wall),cells=np.array(ids))

def main():
    faulthandler.dump_traceback_later(60,repeat=True)
    protocol=read(HERE/'EVENT_FIELD_PROTOCOL.json');assert read(HERE/'FIELD_COMPARABILITY.json')['reference_payload_count']==2
    (HERE/'vtk').mkdir(exist_ok=True)
    parts={name:{} for name in NAMES};pads=[];configuration={}
    for name in NAMES:
        print('LOAD GEOMETRY',name,flush=True)
        grid=pv.read(meshpath(name));tet=grid.cells.reshape(-1,5)[:,1:];points=np.asarray(grid.points,float)
        # Match production's FP64 numerical coordinates without moving points.
        before=points.copy();grid.points=points.copy();assert np.array_equal(grid.points,before)
        geom=TetraGeometry(points,tet);loc=BoundsIndex(points,tet)
        obj=SimpleNamespace(geometry=geom,_locator=loc)
        print('GEOMETRY READY',name,grid.n_cells,geom.candidate_padding_m,flush=True)
        pads.append(geom.candidate_padding_m)
        if name=='M0':
            adj=tetra_adjacency(tet)
            wall=pv.read(paths('M0')[0]/'SV_MESH/mesh-surfaces/WALL.vtp')
            wallids=np.asarray(wall.point_data['GlobalNodeID'])[wall.faces.reshape(-1,4)[:,1:]]-1
            wallkeys={tuple(sorted(v)) for v in wallids}
            wall.save(HERE/'vtk/original_wall.vtp')
        for pid in PIDS:
            event=protocol['events'][str(pid)];raw=np.load(OLD/f'original_{pid}.npz');x=raw['x_m'];t=raw['physical_time_s']
            cache=HERE/f'partitions_{name}_{pid}.json'
            seg=read(cache) if cache.exists() else []
            if seg:
                # Bounded control against the preceding VTK-derived partitions.
                for j in np.linspace(0,len(x)-2,25,dtype=int):
                    new,_=segment_cells(obj,x[j],x[j+1])
                    assert np.allclose(np.asarray(new),np.asarray(seg[j]),rtol=0,atol=1e-12)
            if not seg:
                for j,(a,b) in enumerate(zip(x[:-1],x[1:])):
                    if j%200==0: print('SEGMENT',name,pid,j,flush=True)
                    spans,_=segment_cells(obj,a,b)
                    seg.append(spans)
            parts[name][pid]=seg
            put(f'partitions_{name}_{pid}.json',seg)
            print('PARTITION',name,pid,len(seg),sum(len(s) for s in seg),flush=True)
            if name=='M0':
                crossed=sorted(set(c for spans in seg for _,_,c in spans))
                roi=expand(crossed,adj,2)
                np.savez_compressed(HERE/f'roi_{pid}.npz',cells=roi,crossed=crossed,bounds_m=np.array([points[tet[roi]].min((0,1)),points[tet[roi]].max((0,1))]))
                selected=grid.extract_cells(roi);selected.cell_data['original_canonical_cell']=roi
                selected.save(HERE/f'vtk/original_roi_{pid}.vtu')
                keyseg=int(np.where(t==event['A']['time_s'])[0][0])
                seeds=[c for _,_,c in seg[keyseg]] if pid==822 else [event['A']['tetra']]
                for layer in [0,1,2]:
                    cv=cv_surface(points,tet,expand(seeds,adj,layer),wallkeys)
                    cv['volume_m3']=np.array(geom.volume_m3[cv['cells']].sum())
                    np.savez_compressed(HERE/f'cv_{pid}_L{layer}.npz',**cv)
                    poly=pv.PolyData(cv['triangles_m'].reshape(-1,3),np.c_[np.full(len(cv['triangles_m']),3),np.arange(3*len(cv['triangles_m'])).reshape(-1,3)])
                    poly.cell_data['wall_face']=cv['wall_face'].astype(np.uint8);poly.cell_data['normal_out']=cv['normal_out']
                    poly.save(HERE/f'vtk/cv_{pid}_L{layer}.vtp')
                plane=event['plane'];p=np.asarray(plane['origin_m']);e1=np.asarray(plane['e1']);e2=np.asarray(plane['e2'])
                projected=(points[tet[event['A']['tetra']]]-p)@np.stack([e1,e2],axis=1)
                lo=projected.min(0);hi=projected.max(0);span=hi-lo;lo-=.1*span;hi+=.1*span
                uu,vv=np.meshgrid(np.linspace(lo[0],hi[0],161),np.linspace(lo[1],hi[1],161))
                planar=p+uu.ravel()[:,None]*e1+vv.ravel()[:,None]*e2
                configuration[pid]=dict(plane={**plane,'bounds_m':[lo.tolist(),hi.tolist()]},roi_layers=2,ROI_cell_count=len(roi),crossed_cell_count=len(crossed),CV_seed_cells=seeds,
                  coordinates='Original world SI; plane plots use dot(x-origin,e1/e2); 3D plots translated explicitly in plotting metadata; no geometric deformation')
                np.savez_compressed(HERE/f'plane_{pid}.npz',x_m=planar,axis1_m=uu,axis2_m=vv)
        del geom,loc,obj,grid;gc.collect()
    put('PARTITIONS.json',finite_json(parts));put('SPATIAL_CONFIGURATION.json',configuration)
    put('GEOMETRY_QUERY_ENGINEERING.json',dict(original_failure='VTK FindCellsWithinBounds stalled in fine mesh, repeated after exact FP64 container conversion; precision-only hypothesis rejected.',
      replacement='scipy cKDTree conservative AABB candidates + unchanged exact barycentric clipping',production_code_changed=False,geometry_values_changed=False,
      controls='25 deterministic real segments per cached M0/medium event agree with previous candidate enumeration; brute-force AABB completeness test in targeted tests',
      failed_attempts_retained='monitoring/prepare_coordinates*, included in EXECUTION_RECEIPT'))
    for pid in PIDS:
        e=protocol['events'][str(pid)];raw=np.load(OLD/f'original_{pid}.npz');x=raw['x_m'];t=raw['physical_time_s']
        ds=norm(np.diff(x,axis=0));arc=np.r_[0,np.cumsum(ds)]
        positions=[];slist=[];times=[];spatial=[];temporal=[];stages=[];qorder=[];segmentid=[];spans_record=[]
        for j in range(len(x)-1):
            roots=sorted(set([0.,1.]+[max(0.,min(1.,r)) for name in NAMES for a,b,c in parts[name][pid][j] for r in [a,b]]))
            # Roots from different tessellations are retained; no original point is discarded.
            for a,b in zip(roots[:-1],roots[1:]):
                if b<=a:continue
                spans_record.append([j,a,b])
                for order in [1,2,4]:
                    nodes,weights=np.polynomial.legendre.leggauss(order)
                    for node,w in zip(nodes,weights):
                        alpha=(a+b)/2+(b-a)*node/2;tt=t[j]+alpha*(t[j+1]-t[j])
                        positions.append(x[j]+alpha*(x[j+1]-x[j]));slist.append(arc[j]+alpha*ds[j]);times.append(tt)
                        spatial.append(ds[j]*(b-a)*w/2);temporal.append((t[j+1]-t[j])*(b-a)*w/2)
                        stages.append(min(2,int(np.searchsorted(e['phase_boundaries_s'][1:3],tt,side='right'))));qorder.append(order);segmentid.append(j)
        np.savez_compressed(HERE/f'path_{pid}.npz',x_m=positions,arc_m=slist,time_s=times,arc_weight_m=spatial,time_weight_s=temporal,
          phase=stages,quadrature_order=qorder,original_segment=segmentid,original_x_m=x,original_time_s=t,original_arc_m=arc,union_spans=np.array(spans_record))
        a=np.asarray(e['A']['x_m']);b=np.asarray(e['B']['x_m'])
        np.savez_compressed(HERE/f'rise_{pid}.npz',x_m=a+np.linspace(0,1,1001)[:,None]*(b-a),arc_m=np.linspace(0,np.linalg.norm(b-a),1001))
        profiles=e['normal_profiles'];d=np.array(profiles['distance_m'])
        profile=np.concatenate([np.array(o)+d[:,None]*np.array(n) for o,n in zip(profiles['origins_m'],profiles['normals'])])
        np.savez_compressed(HERE/f'profiles_{pid}.npz',x_m=profile,d_m=np.tile(d,2),profile_index=np.repeat([0,1],len(d)))
    print('COMMON COORDINATES READY',flush=True)
    faulthandler.cancel_dump_traceback_later()

if __name__=='__main__':main()
