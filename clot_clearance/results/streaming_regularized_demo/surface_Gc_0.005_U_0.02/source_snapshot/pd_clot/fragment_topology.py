"""Failure-driven graph lineage, anchor reachability and conservative clearance."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


class BondDamage:
    def __init__(self, count, config):
        self.D = np.zeros(count)
        self.active = np.ones(count, dtype=bool)
        self.g = np.ones(count)
        self.c = config
        if not 0 < config.get('D_break',1.) <= 1: raise ValueError('D_break must be in (0,1]')

    def advance(self, amplitude):
        c=self.c; old=self.D.copy(); old_active=self.active.copy()
        if c['mode']!='cyclic_accumulation' or c['driver']!='bond_stretch_amplitude':
            raise ValueError('Fragmentation demo requires measured cyclic stretch amplitude')
        increment=c['DeltaN']*c['C_damage']*np.maximum(amplitude/c['Q_ref']-1,0)**c['m_damage']
        if not c['enabled']:increment[:]=0
        self.D=np.where(self.active,np.clip(self.D+increment,0,1),self.D)
        self.active &= self.D<c.get('D_break',1.)
        self.g=np.where(self.active,1-self.D,0.)
        broken=np.flatnonzero(old_active & ~self.active)
        return broken, old, increment


def surface_particles(cloud, integrity, config):
    h=cloud.spacing;R=int(np.ceil(cloud.horizon/h))
    offsets=np.array([[i,j,k] for i in range(-R,R+1) for j in range(-R,R+1) for k in range(-R,R+1)],float)*h
    length=np.linalg.norm(offsets,axis=1)
    bulk=np.sum((1-length[(length>0)&(length<cloud.horizon)]/cloud.horizon)**2)*h**3
    weights=cloud.weight*integrity
    degree=np.bincount(cloud.pairs.ravel(),weights=np.repeat(weights,2),minlength=len(cloud.X))*h**3
    ratio=degree/bulk
    exposed=(ratio<config['weighted_bulk_fraction']) & ~cloud.fixed
    return exposed,ratio


class FragmentTracker:
    def __init__(self, cloud, clearance_x):
        self.cloud=cloud; self.clearance_x=clearance_x
        self.labels=np.zeros(len(cloud.X),dtype=np.int64); self.next_id=1
        self.detached_at=np.full(len(cloud.X),-1,dtype=np.int64)
        self.detach_position=np.full_like(cloud.X,np.nan)
        self.cleared=np.zeros(len(cloud.X),dtype=bool)
        self.lineage=[]; self.crossings=[]; self.previous_com={}
        self.previous_positions=None

    def update(self, integrity, x, v, N):
        g=self.cloud;p=g.pairs[integrity>0]
        graph=coo_matrix((np.ones(2*len(p)),(p[:,[0,1]].ravel(),p[:,[1,0]].ravel())),shape=(len(x),len(x))).tocsr()
        number,raw=connected_components(graph,directed=False)
        groups=[np.flatnonzero(raw==k) for k in range(number)]
        parent_groups={}
        for ids in groups:
            parents=np.unique(self.labels[ids])
            if len(parents)!=1:raise AssertionError('Unexpected healing/merge in failure-only bond graph')
            parent_groups.setdefault(int(parents[0]),[]).append(ids)
        labels=np.empty(len(x),dtype=np.int64)
        for parent,children in parent_groups.items():
            # Anchor-bearing child keeps the old ID; otherwise largest child keeps it.
            children.sort(key=lambda ids:(not bool(g.fixed[ids].any()),-len(ids),int(ids.min())))
            for k,ids in enumerate(children):
                fid=parent if k==0 else self.next_id
                if k:self.next_id+=1;self.lineage.append(dict(cycles=int(N),parent=parent,child=fid,particle_ids=ids.tolist()))
                labels[ids]=fid
        self.labels=labels
        anchored=np.unique(labels[g.fixed]);attached=np.isin(labels,anchored)
        new_detached=(~attached)&(self.detached_at<0)
        self.detached_at[new_detached]=N;self.detach_position[new_detached]=x[new_detached]
        components=[]
        for fid in np.unique(labels):
            ids=np.flatnonzero(labels==fid);volume=float(g.volume[ids].sum())
            com=np.average(x[ids],axis=0,weights=g.volume[ids]);vel=np.average(v[ids],axis=0,weights=g.volume[ids])
            is_attached=bool(attached[ids[0]])
            travel=np.zeros(3) if is_attached else com-np.average(self.detach_position[ids],axis=0,weights=g.volume[ids])
            previous=None if self.previous_positions is None else np.average(self.previous_positions[ids],axis=0,weights=g.volume[ids])
            # Unique material ledger prevents double counting when cleared fragments split.
            if not is_attached and previous is not None and previous[0]<self.clearance_x<=com[0]:
                if np.any(~self.cleared[ids]):
                    fresh=ids[~self.cleared[ids]]
                    if len(fresh):
                        self.crossings.append(dict(cycles=int(N),fragment_id=int(fid),particle_ids=fresh.tolist(),
                            newly_cleared_volume_m3=float(g.volume[fresh].sum()),center_of_mass_m=com.tolist(),
                            inherited_or_new_component=previous is None))
                    self.cleared[ids]=True
            components.append(dict(fragment_id=int(fid),particle_ids=ids.tolist(),attached_to_base=is_attached,
                particle_count=len(ids),volume_m3=volume,center_of_mass_m=com.tolist(),velocity_of_center_of_mass_m_s=vel.tolist(),
                bounding_box_m=[x[ids].min(axis=0).tolist(),x[ids].max(axis=0).tolist()],
                travel_vector_since_detachment_m=travel.tolist(),travel_distance_m=float(np.linalg.norm(travel))))
            self.previous_com[int(fid)]=com.copy()
        total=float(g.volume.sum());free=[row for row in components if not row['attached_to_base']]
        stats=dict(total_components=len(components),attached_components=len(components)-len(free),detached_components=len(free),
            attached_volume_fraction=float(g.volume[attached].sum()/total),detached_volume_fraction=float(g.volume[~attached].sum()/total),
            cleared_volume_fraction=float(g.volume[self.cleared].sum()/total),largest_component_volume_m3=max(row['volume_m3'] for row in components),
            detached_fragment_volumes_m3=[row['volume_m3'] for row in free],
            maximum_fragment_travel_m=max([row['travel_distance_m'] for row in free],default=0.))
        self.previous_positions=x.copy()
        return attached,components,stats
