"""Same full B64 arrays: gzip1 versus lossless byte shuffle plus gzip1."""

import argparse
from time import perf_counter
import h5py
import numpy as np
from .common import BASE, F10C, identity, sha, write


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    source = F10C / "exclusive_v1_batch/B64_batch/trajectory.h5"
    begin = perf_counter()
    nodes = []
    seen = {}

    def read(group, prefix=""):
        for key in group:
            obj = group[key]
            name = prefix + key
            address = h5py.h5o.get_info(obj.id).addr
            if address in seen:
                nodes.append(dict(path=name, link=seen[address]))
                continue
            seen[address] = name
            node = dict(path=name, attrs=dict(obj.attrs))
            nodes.append(node)
            if isinstance(obj, h5py.Group):
                node["group"] = True
                read(obj, name + "/")
            else:
                node.update(data=obj[()], dtype=obj.dtype, compressed=obj.compression is not None)

    with h5py.File(source) as h:
        attrs = dict(h.attrs)
        read(h)
    read_s = perf_counter() - begin
    rows = []
    for mode in ["legacy", "optimized"]:
        start = perf_counter()
        path = out / (mode + ".h5")
        with h5py.File(path, "x") as h:
            h.attrs.update(attrs)
            for node in nodes:
                name = node["path"]
                if "link" in node:
                    h[name] = h[node["link"]]
                    continue
                if node.get("group"):
                    obj = h.require_group(name)
                else:
                    opts = (
                        dict(compression="gzip", compression_opts=1, shuffle=mode == "optimized")
                        if node["compressed"]
                        else {}
                    )
                    obj = h.create_dataset(name, data=node["data"], dtype=node["dtype"], **opts)
                obj.attrs.update(node["attrs"])
        wall = perf_counter() - start
        start = perf_counter()
        count = 0
        with h5py.File(path) as h:
            for node in nodes:
                if "data" not in node:
                    continue
                aa = np.asarray(node["data"])
                bb = np.asarray(h[node["path"]][()])
                assert aa.dtype == bb.dtype
                assert (
                    np.array_equal(aa, bb, equal_nan=True)
                    if aa.dtype.kind in "fc"
                    else np.array_equal(aa, bb)
                )
                if aa.dtype.kind not in "OUS":
                    assert aa.tobytes() == bb.tobytes()
                count += 1
        rows.append(
            dict(
                mode=mode,
                hdf5_write_s=wall,
                read_verify_s=perf_counter() - start,
                bytes=path.stat().st_size,
                datasets_bitwise=count,
                sha256=sha(path),
            )
        )
    write(
        out / "RESULTS.json",
        dict(
            status="PASS",
            source_sha256=sha(source),
            input_read_s=read_s,
            rows=rows,
            scope="Full original arrays, dtype, attributes and hard links preserved. Imports excluded; filesystem cache not flushed. Compression is lossless; no output omitted.",
        ),
    )
    print(rows, flush=True)


if __name__ == "__main__":
    main()
