"""Deterministic KMeans and exploratory cluster metrics without sklearn."""

from __future__ import annotations

import numpy as np

from .sampling_types import ClusteringResult


def _squared_distances(values: np.ndarray, centers: np.ndarray) -> np.ndarray:
    return np.sum((values[:, None, :] - centers[None, :, :]) ** 2, axis=2)


def _initial_centers(values: np.ndarray, n_clusters: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    chosen = [int(rng.integers(0, len(values)))]
    nearest = np.sum((values - values[chosen[0]]) ** 2, axis=1)
    while len(chosen) < n_clusters:
        nearest[chosen] = 0.0
        total = float(np.sum(nearest))
        if total <= 1.0e-15:
            candidate = next(index for index in range(len(values)) if index not in chosen)
        else:
            candidate = int(rng.choice(len(values), p=nearest / total))
            if candidate in chosen:
                candidate = int(np.argmax(nearest))
        chosen.append(candidate)
        nearest = np.minimum(nearest, np.sum((values - values[candidate]) ** 2, axis=1))
    return values[chosen].copy()


def _silhouette(values: np.ndarray, labels: np.ndarray, n_clusters: int) -> float | None:
    if n_clusters <= 1 or len(values) <= n_clusters:
        return None
    distances = np.linalg.norm(values[:, None, :] - values[None, :, :], axis=2)
    scores: list[float] = []
    for index, label in enumerate(labels.tolist()):
        same = np.flatnonzero(labels == label)
        same = same[same != index]
        if not len(same):
            scores.append(0.0)
            continue
        within = float(np.mean(distances[index, same]))
        outside = [
            float(np.mean(distances[index, labels == other]))
            for other in range(n_clusters)
            if other != label and np.any(labels == other)
        ]
        nearest_other = min(outside)
        denominator = max(within, nearest_other)
        scores.append((nearest_other - within) / denominator if denominator > 0 else 0.0)
    return float(np.mean(scores))


def deterministic_kmeans(
    values: np.ndarray,
    *,
    n_clusters: int,
    feature_names: tuple[str, ...],
    seed: int,
    max_iter: int,
) -> ClusteringResult:
    matrix = np.asarray(values, dtype=float)
    if matrix.ndim != 2 or not len(matrix):
        raise ValueError("KMeans requires a non-empty 2-D feature matrix")
    cluster_count = min(max(1, int(n_clusters)), len(matrix))
    centers = _initial_centers(matrix, cluster_count, seed)
    assignments = np.full(len(matrix), -1, dtype=np.int64)
    for _ in range(max_iter):
        squared = _squared_distances(matrix, centers)
        updated_assignments = np.argmin(squared, axis=1).astype(np.int64)
        if np.array_equal(updated_assignments, assignments):
            break
        assignments = updated_assignments
        updated_centers = centers.copy()
        for cluster_id in range(cluster_count):
            members = matrix[assignments == cluster_id]
            if len(members):
                updated_centers[cluster_id] = np.mean(members, axis=0)
            else:
                own_distance = squared[np.arange(len(matrix)), assignments]
                replacement = int(np.argmax(own_distance))
                updated_centers[cluster_id] = matrix[replacement]
                assignments[replacement] = cluster_id
        if np.allclose(updated_centers, centers, rtol=0.0, atol=1.0e-10):
            centers = updated_centers
            break
        centers = updated_centers
    squared = _squared_distances(matrix, centers)
    assignments = np.argmin(squared, axis=1).astype(np.int64)
    distances = np.sqrt(squared[np.arange(len(matrix)), assignments])
    sizes = tuple(int(np.count_nonzero(assignments == index)) for index in range(cluster_count))
    return ClusteringResult(
        method="kmeans",
        n_clusters=cluster_count,
        feature_names=feature_names,
        scaled_features=matrix.copy(),
        assignments=assignments,
        centers=centers,
        distances_to_center=distances,
        inertia=float(np.sum(distances ** 2)),
        silhouette_score=_silhouette(matrix, assignments, cluster_count),
        cluster_sizes=sizes,
    )


def exploratory_cluster_scan(
    values: np.ndarray,
    *,
    feature_names: tuple[str, ...],
    candidate_k: tuple[int, ...],
    seed: int,
    max_iter: int,
) -> list[dict[str, float | int | None]]:
    results: list[dict[str, float | int | None]] = []
    for cluster_count in sorted(set(candidate_k)):
        if cluster_count < 1 or cluster_count > len(values):
            continue
        result = deterministic_kmeans(
            values,
            n_clusters=cluster_count,
            feature_names=feature_names,
            seed=seed,
            max_iter=max_iter,
        )
        results.append(
            {
                "n_clusters": result.n_clusters,
                "inertia": result.inertia,
                "silhouette_score": result.silhouette_score,
            }
        )
    return results

