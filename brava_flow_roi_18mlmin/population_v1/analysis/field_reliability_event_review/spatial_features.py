"""Describe where the fixed historical path's background speed gain occurs."""
from common import *
from analyze import quantities,pair_summary
import pyvista as pv

def main():
    protocol=read(HERE/'EVENT_FIELD_PROTOCOL.json');out={}
    for pid in PIDS:
        e=protocol['events'][str(pid)];p=np.load(HERE/f'path_{pid}.npz');a=np.array(e['A']['x_m']);b=np.array(e['B']['x_m'])
        ta=e['A']['time_s'];tb=e['B']['time_s'];t=p['original_time_s'];x=p['original_x_m'];arc=p['original_arc_m']
        ia=int(np.flatnonzero(t==ta)[0]);ib=int(np.flatnonzero(t==tb)[0]);spans=p['union_spans'];mid=p['quadrature_order']==1
        keep=(spans[:,0]>=ia)&(spans[:,0]<ib);ss=spans[keep];seg=ss[:,0].astype(int)
        dvec=(ss[:,2]-ss[:,1])[:,None]*(x[seg+1]-x[seg]);length=norm(dvec)
        r=dict(event_path_length_um=float((arc[ib]-arc[ia])*1e6),fields={})
        for name in NAMES:
            z=np.load(HERE/f'{name}_path_{pid}.npz');f=z['flow'][mid][keep];cells=z['cell'][mid][keep];g=f[:,4:13].reshape(-1,3,3)
            du=np.einsum('nij,nj->ni',g,dvec);u0=f[:,:3]-du/2;u1=f[:,:3]+du/2;gain=(norm(u1)-norm(u0))*1000
            bycell=[]
            for cid in np.unique(cells):
                sel=cells==cid
                bycell.append(dict(cell=int(cid),signed_speed_gain_mm_s=float(gain[sel].sum()),path_length_um=float(length[sel].sum()*1e6),
                   first_arc_from_A_um=float((arc[seg[sel][0]]+ss[sel][0,1]*norm(x[seg[sel][0]+1]-x[seg[sel][0]])-arc[ia])*1e6)))
            bycell.sort(key=lambda v:v['signed_speed_gain_mm_s'],reverse=True)
            profile=np.load(HERE/f'{name}_profiles_{pid}.npz');coord=np.load(HERE/f'profiles_{pid}.npz');wall=np.load(HERE/f'wall_profiles_{pid}.npz')
            normal=[]
            for k in [0,1]:
                valid=(coord['profile_index']==k)&wall['original_normal_ray_valid']&profile['valid']&(coord['d_m']>0)
                normal.append(dict(ray=k,valid_positive_distance_count=int(valid.sum()),max_distance_um=float(coord['d_m'][valid].max()*1e6) if valid.any() else None,
                   cell_count=int(len(np.unique(profile['cell'][valid]))),cells=np.unique(profile['cell'][valid]).tolist()))
            r['fields'][name]=dict(all_cells_signed_gain=bycell,net_scalar_gain_mm_s=float(gain.sum()),negative_gain_mm_s=float(gain[gain<0].sum()),normal_profiles=normal)
        plane=np.load(HERE/f'plane_{pid}.npz');wall=np.load(HERE/f'wall_plane_{pid}.npz');n=wall['normal']
        if pid==822:n=np.broadcast_to(e['plane']['e1'],n.shape)
        raw={name:np.load(HERE/f'{name}_plane_{pid}.npz') for name in NAMES};q={name:quantities(v['flow'],n) for name,v in raw.items()}
        mask=raw['M0']['valid'];floor=max(read(HERE/f'{name}_SAMPLING.json')['velocity_roundoff_scale_m_s'] for name in NAMES)
        r['independent_plane']=dict(grid_points=len(mask),original_in_domain_points=int(mask.sum()),reference_comparisons={},
          weighting='Uniform fixed spatial raster; diagnostic of this plane, not full-volume average',outside_original_points_retained=int((~mask).sum()))
        for first,second in [(NAMES[0],NAMES[1]),(NAMES[0],NAMES[2]),(NAMES[1],NAMES[2])]:
            r['independent_plane']['reference_comparisons'][first+'__'+second]=pair_summary(q[first],q[second],np.ones(len(mask)),mask,floor)
        out[pid]=r
    inlet=pv.read(paths('M0')[0]/'SV_MESH/mesh-surfaces/INLET.vtp');tri=inlet.points[inlet.faces.reshape(-1,4)[0,1:]]
    n=unit(np.cross(tri[1]-tri[0],tri[2]-tri[0]));origin=inlet.points.mean(0)
    e=protocol['events']['822'];dist=[float(abs((np.array(e[k]['x_m'])-origin)@n)*1e6) for k in ['A','B']]
    out[822]['inlet_geometry_context']=dict(distance_to_inlet_plane_um_A_B=dist,
        cap_planarity_deviation_um=float(np.max(abs((inlet.points-origin)@n))*1e6),
        interpretation='A and B lie in the inlet-adjacent region of the declared flat-profile artificial port. Not an independent physiological acceleration measurement; BC was identical in all fields.')
    put('SPATIAL_FEATURES.json',finite_json(out));print('SPATIAL FEATURES COMPLETE')

if __name__=='__main__':main()
