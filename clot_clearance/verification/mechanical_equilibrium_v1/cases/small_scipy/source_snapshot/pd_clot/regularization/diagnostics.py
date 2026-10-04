"""Diagnostic-only component resolution, localization and unique clearance."""
import numpy as np

ATTACHED = 0
RESOLVED_FRAGMENT = 1
UNDER_RESOLVED_DEBRIS = 2
SINGLETON_NUMERICAL_COMPONENT = 3
NAMES = ['ATTACHED', 'RESOLVED_FRAGMENT', 'UNDER_RESOLVED_DEBRIS', 'SINGLETON_NUMERICAL_COMPONENT']


def classify(cloud, integrity, ranks, attached, labels, components, config):
    category = np.zeros(len(cloud.X), np.uint8)
    rows = []
    active_pairs = cloud.pairs[integrity > 0]
    if np.any(labels[active_pairs[:, 0]] != labels[active_pairs[:, 1]]):
        raise AssertionError('Active bond crosses connected components')
    for component in components:
        row = dict(component)
        ids = np.asarray(row['particle_ids'])
        fid = row['fragment_id']
        internal = int(np.sum(labels[active_pairs[:, 0]] == fid))
        supported = float(np.sum(cloud.volume[ids]*(ranks[ids] >= config['minimum_support_rank']))/cloud.volume[ids].sum())
        if row['attached_to_base']:
            kind = ATTACHED
        elif len(ids) == 1:
            kind = SINGLETON_NUMERICAL_COMPONENT
        elif (len(ids) >= config['minimum_particles'] and row['volume_m3'] >= config['minimum_volume_m3']
              and internal >= config['minimum_internal_bonds'] and supported >= config['minimum_supported_fraction']):
            kind = RESOLVED_FRAGMENT
        else:
            kind = UNDER_RESOLVED_DEBRIS
        category[ids] = kind
        row.update(resolution_class=NAMES[kind], internal_active_bonds=internal,
                   supported_volume_fraction=supported, equivalent_spherical_diameter_m=(6*row['volume_m3']/np.pi)**(1/3))
        rows.append(row)
    total = float(cloud.volume.sum())
    detached = float(cloud.volume[~attached].sum())
    volumes = [float(cloud.volume[category == k].sum()) for k in range(4)]
    lowrank = float(cloud.volume[(~attached) & (ranks < 2)].sum())
    singleton_ratio = volumes[3]/detached if detached else None
    lowrank_ratio = lowrank/detached if detached else None
    inadequate = bool(detached and (singleton_ratio > config['singleton_flag_threshold'] or lowrank_ratio > config['lowrank_flag_threshold']))
    stats = dict(detached_total_volume_m3=detached, resolved_fragment_volume_m3=volumes[1],
                 under_resolved_debris_volume_m3=volumes[2], singleton_volume_m3=volumes[3],
                 resolved_fragment_volume_fraction=volumes[1]/total, under_resolved_debris_fraction=volumes[2]/total,
                 singleton_volume_fraction=volumes[3]/total, singleton_detached_volume_fraction=singleton_ratio,
                 lowrank_detached_volume_m3=lowrank, P_singleton=singleton_ratio, P_lowrank=lowrank_ratio,
                 largest_resolved_fragment_volume_m3=max([r['volume_m3'] for r in rows if r['resolution_class']==NAMES[1]], default=0.),
                 fragment_resolution_status=('NOT_ASSESSED_NO_DETACHMENT' if not detached else
                                             'FRAGMENT_RESOLUTION_INADEQUATE' if inadequate else 'FRAGMENT_RESOLUTION_ADEQUATE'),
                 rank_particle_fractions=[float(np.mean(ranks == k)) for k in range(4)],
                 detached_rank_volume_fractions=[float(cloud.volume[(~attached)&(ranks==k)].sum()/detached) if detached else None for k in range(4)])
    assert np.isclose(sum(volumes), total, rtol=1e-12)
    assert np.isclose(sum(volumes[1:]), detached, rtol=1e-12, atol=1e-30)
    return category, rows, stats


class ClassifiedClearance:
    def __init__(self, count):
        self.crossing_class = np.full(count, -1, np.int8)

    def update(self, cleared, category, volume):
        fresh = cleared & (self.crossing_class < 0)
        self.crossing_class[fresh] = category[fresh]
        if np.any(self.crossing_class[cleared] == ATTACHED):
            raise AssertionError('Attached material cannot be recorded as cleared')
        total = float(volume.sum())
        resolved = volume[self.crossing_class == RESOLVED_FRAGMENT].sum()/total
        under = volume[self.crossing_class >= UNDER_RESOLVED_DEBRIS].sum()/total
        all_clear = volume[self.crossing_class >= 0].sum()/total
        assert np.isclose(resolved+under, all_clear)
        return dict(total_clearance_volume_fraction=float(all_clear),
                    resolved_fragment_clearance_fraction=float(resolved),
                    under_resolved_debris_clearance_fraction=float(under))


def localization(cloud, damage, integrity, attached, center, radius):
    distance = np.linalg.norm(cloud.X-np.asarray(center), axis=1)
    inside = distance <= radius
    bond_midpoint = cloud.X[cloud.pairs].mean(axis=1)
    bond_inside = np.linalg.norm(bond_midpoint-center, axis=1) <= radius
    Dvol = damage*cloud.volume
    broken = integrity == 0
    detached = cloud.volume[~attached].sum()
    return dict(fraction_of_damage_inside_influence_region=float(Dvol[inside].sum()/Dvol.sum()) if Dvol.sum() else None,
                fraction_of_broken_bonds_inside_influence_region=float(np.mean(bond_inside[broken])) if np.any(broken) else None,
                fraction_of_detached_volume_originating_in_influence_region=float(cloud.volume[inside & ~attached].sum()/detached) if detached else None,
                damage_weighted_distance_from_source_m=float(np.dot(Dvol, distance)/Dvol.sum()) if Dvol.sum() else None,
                initial_volume_in_influence_region_fraction=float(cloud.volume[inside].sum()/cloud.volume.sum()))
