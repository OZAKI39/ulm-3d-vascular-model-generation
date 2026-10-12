"""Bounded offline CPU profile of existing B64 data; no physics integration."""

import argparse
import cProfile
import io
import json
import platform
import pstats
from time import perf_counter
import h5py
import numpy as np
from .common import BASE, F10C, identity, sha, write
from .postprocess_reference import process


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    path = F10C / "exclusive_v1_batch/B64_batch/trajectory.h5"
    t = perf_counter()
    with h5py.File(path) as h:
        x = h["positions_m"][:]
        times = h["time_s"][:]
        ei, et = h["edge_events/batch_i_j_linked"][:], h["edge_events/time_s"][:]
        radius = float(h["radii_m"][0])
    read_s = perf_counter() - t
    rows = []
    for repeat in ["first_after_import", "warm"]:
        metrics, lineage, morph, top, timing = process(x, radius, times, ei, et)
        rows.append(dict(repeat=repeat, **timing))
        print(rows[-1], flush=True)
    pr = cProfile.Profile()
    pr.enable()
    process(x, radius, times, ei, et)
    pr.disable()
    pr.dump_stats(str(out / "legacy_cpu.prof"))
    stream = io.StringIO()
    stats = pstats.Stats(pr, stream=stream).sort_stats("cumulative")
    stats.print_stats(45)
    (out / "PROFILE.txt").write_text(stream.getvalue())
    funcs = [
        dict(
            file=k[0], line=k[1], name=k[2], primitive_calls=v[0], calls=v[1], self_s=v[2], cumulative_s=v[3]
        )
        for k, v in stats.stats.items()
    ]
    funcs.sort(key=lambda row: row["cumulative_s"], reverse=True)
    t = perf_counter()
    serialized = json.dumps(dict(lineage=lineage, topology=top), allow_nan=False)
    json_s = perf_counter() - t
    t = perf_counter()
    target = out / "reference_postprocessing.h5"
    with h5py.File(target, "x") as h:
        for group, arrays in [("metrics", metrics), ("morphology", morph)]:
            for key, value in arrays.items():
                h.create_dataset(group + "/" + key, data=value, compression="gzip", compression_opts=1)
        h.attrs["events_json"] = serialized
    write_s = perf_counter() - t
    t = perf_counter()
    with h5py.File(target) as h:
        for group, arrays in [("metrics", metrics), ("morphology", morph)]:
            for key, value in arrays.items():
                np.testing.assert_array_equal(h[group + "/" + key][:], value)
    verify_s = perf_counter() - t
    write(
        out / "RESULTS.json",
        dict(
            status="PASS",
            source=str(path),
            input_sha256=sha(path),
            host=platform.platform(),
            numpy=np.__version__,
            input_read_s=read_s,
            shape=x.shape,
            repeats=rows,
            profile=funcs[:100],
            json_serialize_s=json_s,
            hdf5_write_s=write_s,
            output_verify_s=verify_s,
            bytes=target.stat().st_size,
            not_separately_measured=[
                "GPU state readback (offline input)",
                "centroid/Rg/tensor inside original shape/frame_metrics; see inclusive profile, not additive wall components",
                "original device quality checks (no GPU work)",
            ],
            timing_scope="Filesystem cache not flushed; imports excluded from repeated CPU walls; instrumented cProfile cost is separate",
        ),
    )
    print(stream.getvalue(), flush=True)


if __name__ == "__main__":
    main()
