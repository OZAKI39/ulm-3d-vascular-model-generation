"""Lossless bounded-memory compression of Stage2.5 cache only."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import time
from mevo_mb.bubbles.adaptive_table import compress_raw
from common import atomic_json


def main():
    start = time.perf_counter()
    paths = sorted(Path("results/stage2_5/response_table/node_cache").glob("*/attempt_*"))
    original = 0
    files = 0
    with ProcessPoolExecutor(max_workers=4) as pool:
        for rows in pool.map(compress_raw, paths, chunksize=8):
            original += sum(r["original_bytes"] for r in rows)
            files += len(rows)
    atomic_json(
        "results/stage2_5/response_table/LOSSLESS_COMPRESSION.json",
        dict(
            raw_files=files,
            original_bytes=original,
            elapsed_s=time.perf_counter() - start,
            note="All decompressed bytes verified before replacement; Stage1/2 files untouched",
        ),
    )
    print(files, original, flush=True)


if __name__ == "__main__":
    main()
