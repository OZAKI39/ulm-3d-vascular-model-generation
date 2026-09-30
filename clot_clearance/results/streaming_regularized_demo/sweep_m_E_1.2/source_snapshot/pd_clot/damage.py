import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def instantaneous(integrity, stretch, s1, s2):
    if not 1 <= s1 < s2:
        raise ValueError('Require 1 <= s1 < s2 (stretch is a length ratio)')
    z = np.clip((stretch-s1)/(s2-s1), 0, 1)
    candidate = 1 - z*z*(3-2*z)
    return np.minimum(integrity, candidate)


def accumulate(integrity, amplitude, c):
    if c['driver'] != 'bond_stretch_amplitude':
        raise ValueError('Only measured bond_stretch_amplitude is implemented')
    if min(c['C_damage'], c['DeltaN']) < 0 or c['Q_ref'] <= 0 or c['m_damage'] <= 0:
        raise ValueError('Invalid cyclic damage parameters')
    inc = c['DeltaN'] * c['C_damage'] * np.maximum(amplitude/c['Q_ref'] - 1, 0)**c['m_damage']
    return np.maximum(0, integrity-inc)


def particle_damage(cloud, integrity):
    ids = cloud.pairs.ravel()
    degree = np.bincount(ids, minlength=len(cloud.X))
    sums = np.bincount(ids, weights=np.repeat(1-integrity, 2), minlength=len(cloud.X))
    maximum = np.zeros(len(cloud.X))
    np.maximum.at(maximum, ids, np.repeat(1-integrity, 2))
    return sums/degree, maximum


def fragments(cloud, integrity, threshold):
    if not 0 <= threshold < 1:
        raise ValueError('connectivity_threshold must be in [0,1)')
    p = cloud.pairs[integrity > threshold]
    graph = coo_matrix((np.ones(2*len(p)), (p[:, [0,1]].ravel(), p[:, [1,0]].ravel())),
                       shape=(len(cloud.X), len(cloud.X))).tocsr()
    count, labels = connected_components(graph, directed=False)
    volumes = np.bincount(labels, weights=cloud.volume)
    attached = np.unique(labels[cloud.fixed])
    sizes = np.bincount(labels)
    active = np.bincount(p.ravel(), minlength=len(cloud.X))
    return labels, active, dict(number_of_fragments=int(count),
        fragment_sizes=sizes.tolist(), fragment_volumes_m3=volumes.tolist(),
        largest_fragment_volume_m3=float(volumes.max()),
        largest_fragment_fraction=float(volumes.max()/volumes.sum()),
        detached_volume_fraction=float(volumes[~np.isin(np.arange(count),attached)].sum()/volumes.sum()))
