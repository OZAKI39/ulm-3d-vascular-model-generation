"""Original-ID connected-component lineage and volume-weighted morphology."""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from .config import RHO, V


class Lineage:

    def __init__(self, ref, data=None):
        self.ref = ref
        self.registry = (
            {} if data is None else {int(k): v for k, v in data["registry"].items()}
        )
        self.splits = [] if data is None else data["splits"]
        self.lookup = {tuple(q["atom_ids"]): k for k, q in self.registry.items()}
        self.current = []

    def update(self, top, s):
        lab = top["labels"]
        labels = np.flatnonzero(~top["attached"])
        order = np.argsort(lab, kind="stable")
        starts = np.r_[0, np.cumsum(np.bincount(lab))]
        old = self.current
        self.current = []
        for label in labels:
            ids = np.sort(order[starts[label] : starts[label + 1]])
            key = tuple((ids + 1).tolist())
            if key not in self.lookup:
                cid = len(self.registry) + 1
                self.lookup[key] = cid
                vol = float(self.ref["V"][ids].sum())
                self.registry[cid] = dict(
                    atom_ids=list(key),
                    nodes=len(ids),
                    reference_volume_m3=vol,
                    mass_kg=RHO * vol,
                    equivalent_nominal_nodes=vol / V,
                    birth_time_s=s["time"],
                    birth_step=s["step"],
                    last_seen_time_s=s["time"],
                    split_time_s=None,
                )
            cid = self.lookup[key]
            self.registry[cid]["last_seen_time_s"] = s["time"]
            self.current.append((int(label), cid, ids))
        cidmap = np.zeros(len(top["attached"]), int)
        for label, cid, ids in self.current:
            cidmap[label] = cid
        for _, cid, ids in old:
            children = cidmap[np.unique(lab[ids])].tolist()
            if len(children) > 1:
                self.splits.append(
                    dict(
                        time_s=s["time"],
                        step=s["step"],
                        parent_component_id=cid,
                        children_component_ids=children,
                    )
                )
                self.registry[cid]["split_time_s"] = s["time"]

    def data(self):
        return dict(registry=self.registry, splits=self.splits)

    def records(self, s, top, nodework):
        lab = top["labels"]
        b = self.ref["bonds"]
        live = s["d"] < 1
        livecounts = np.bincount(lab[b[live, 0]], minlength=len(top["attached"]))
        out = []
        for label, cid, ids in self.current:
            vol = self.ref["V"][ids]
            weight = vol / vol.sum()
            x = s["x"][ids]
            v = s["velocity"][ids]
            center = weight @ x
            vc = weight @ v
            dx = x - center
            cov = (dx * weight[:, None]).T @ dx
            eig, axes = np.linalg.eigh(cov)
            eig = eig[::-1]
            axes = axes[:, ::-1]
            lengths = np.sqrt(np.maximum(eig, 0))
            ratio = lengths / max(lengths[0], 1e-300)
            typ = (
                "SMALL_DEBRIS"
                if len(ids) < 20
                else (
                    "FILAMENT"
                    if ratio[1] < 0.2
                    else "THIN_SHEET" if ratio[2] < 0.2 else "COHESIVE_BLOCK"
                )
            )
            out.append(
                dict(
                    component_id=cid,
                    **self.registry[cid],
                    centroid_m=center,
                    centroid_velocity_m_s=vc,
                    principal_covariance_axes=axes,
                    principal_variances_m2=eig,
                    principal_radii_m=lengths,
                    aspect_ratios=ratio,
                    radius_of_gyration_m=float(np.sqrt(eig.sum())),
                    equivalent_sphere_diameter_m=float(
                        (6 * vol.sum() / np.pi) ** (1 / 3)
                    ),
                    internal_alive_bonds=int(livecounts[label]),
                    internal_graph_density=(
                        float(2 * livecounts[label] / (len(ids) * (len(ids) - 1)))
                        if len(ids) > 1
                        else None
                    ),
                    internal_relative_kinetic_J=float(
                        0.5 * RHO * np.sum(vol[:, None] * (v - vc) ** 2)
                    ),
                    cumulative_member_external_work_J=float(nodework[ids].sum()),
                    morphology=typ
                )
            )
        return out


def components(ref, d):
    bonds = ref["bonds"][d < 1]
    N = len(ref["X"])
    n = int(ref["clot_count"])
    graph = csr_matrix(
        (
            np.ones(2 * len(bonds)),
            (np.r_[bonds[:, 0], bonds[:, 1]], np.r_[bonds[:, 1], bonds[:, 0]]),
        ),
        shape=(N, N),
    )
    count, labels = connected_components(graph, directed=False)
    attached = np.zeros(count, bool)
    attached[labels[n:]] = True
    return dict(labels=labels, attached=attached, free=~attached[labels])
