"""Check all 149 source states of the latest fracture video, without dynamics."""
import json
from pathlib import Path
import numpy as np
from verify_export import ROOT, META, records, sha


def main():
    exports = {r["source"]: r for r in records() if "export" in r}
    source = "verification/y_clot_irregular_current_method_v4/20261020T041100Z/"
    hashes = json.loads((ROOT / "verification/y_clot_fracture_reference_style_v2/"
                         "20261020T094240Z/SOURCE_HASHES.json").read_text())
    ref = np.load(ROOT / source / "provenance/reference.npz", allow_pickle=False)
    X, bonds, ids = ref["X"], ref["bonds"], ref["atom_ids"]
    n = int(ref["clot_count"])
    states = []
    for rel, want in hashes.items():
        path = ROOT / rel
        if path.exists():
            assert sha(path) == want, rel
        else:
            row = exports[rel]
            assert row["source_sha256"] == want, rel
            path = ROOT / row["export"]
            assert sha(path) == row["export_sha256"], rel
        if "/keys/" in rel or rel.endswith("/last_accepted.npz"):
            with np.load(path, allow_pickle=False) as z:
                states.append((float(z["time"]), z["x"].copy(), z["d"].copy(), z["atom_ids"].copy()))
    states.sort(key=lambda s: s[0])
    assert len(states) == 149 and abs(states[-1][0] - .0036) < 1e-14
    last = np.zeros(len(bonds))
    for t, x, d, atom_ids in states:
        assert np.array_equal(ids, atom_ids) and np.isfinite(x).all()
        assert np.array_equal(x[n:], X[n:]) and np.all(d >= last - 1e-12)
        assert np.all((d >= 0) & (d <= 1))
        last = d
    parent = list(range(len(X)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in bonds[last < 1]:
        a, b = find(int(i)), find(int(j))
        if a != b:
            parent[a] = b
    labels = [find(i) for i in range(len(X))]
    anchored = set(labels[n:])
    from collections import Counter
    free = Counter(label for label in labels[:n] if label not in anchored)
    result = dict(status="PASS", states=len(states), sources=len(hashes),
                  physical_end_s=states[-1][0], nodes=len(X), clot_nodes=n,
                  reference_bonds=len(bonds), terminal_free_nodes=sum(free.values()),
                  terminal_free_components=len(free), terminal_largest_free=max(free.values()),
                  terminal_target_20_to_200=sum(20 <= x <= 200 for x in free.values()),
                  exact_selected_arrays=True, new_physics_steps=0)
    assert (result["terminal_free_nodes"], result["terminal_free_components"],
            result["terminal_largest_free"], result["terminal_target_20_to_200"]) == (1743, 912, 461, 1)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
