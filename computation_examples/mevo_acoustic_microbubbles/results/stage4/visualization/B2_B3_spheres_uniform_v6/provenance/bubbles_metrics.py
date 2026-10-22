"""Spatial clustering metrics; clustering never implies fusion or coalescence."""
import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist, squareform


def aggregation_metrics(positions_m, radii_m, config):
    pairs = pdist(positions_m)
    distance = squareform(pairs)
    labels_count, labels = connected_components(distance <= config["cluster_threshold_m"], directed=False)
    np.fill_diagonal(distance, np.inf)
    nn = distance.min(axis=1)
    gap = distance-radii_m[:, None]-radii_m[None, :]
    sizes = np.bincount(labels)
    hist, edges = np.histogram(pairs, bins=config["pair_histogram_bins"], range=(0, config["pair_histogram_max_m"]))
    return {"mean_nearest_neighbor_m": float(nn.mean()), "median_nearest_neighbor_m": float(np.median(nn)),
            "cluster_count": int(labels_count), "largest_cluster_size": int(sizes.max()),
            "largest_cluster_fraction": float(sizes.max()/len(positions_m)),
            "near_contact_fraction": float(np.mean(np.min(gap, axis=1) <= config["near_contact_gap_m"])),
            "minimum_surface_gap_m": float(gap.min()), "pair_distance_mean_m": float(pairs.mean())}, labels, hist, edges
