"""Adaptive APECSS database with durable per-node cache and rejected leaves.

No radial equations are implemented here. Cell tests are empirical evidence;
fixed shell manifold and calibration assumptions are inherited from Stage2.
"""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import fcntl
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import time

import h5py
import numpy as np

from .apecss import APECSS_COMMIT, ApecssBubbleResponseModel, RadialCase
from .parameters import SonoVueParameterDistribution
from .response_table import ResponseTableBoundsError, ResponseTableBubbleModel

PROFILES = {"strict": 0.02, "production": 0.05, "exploratory": 0.10}
CORNERS = np.array(list(itertools.product([0.0, 1.0], repeat=3)))
CHECKS = np.array(
    [
        [0.5, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.75, 0.75, 0.75],
        [0.25, 0.5, 0.5],
        [0.75, 0.5, 0.5],
        [0.5, 0.25, 0.5],
        [0.5, 0.75, 0.5],
        [0.5, 0.5, 0.25],
        [0.5, 0.5, 0.75],
    ]
)


def durable_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with tmp.open("w") as stream:
        json.dump(data, stream, sort_keys=True, allow_nan=False)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)
    fd = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def node_key(case, executable_sha256):
    values = asdict(case)
    values.pop("provenance")
    # Equivalent int/float physical inputs have one key. Existing cache values
    # already use these declared types, so normalization preserves their keys.
    values = {k: int(v) if k in ("cycles", "samples_per_cycle") else float(v) for k, v in values.items()}
    identity = dict(
        model="RPAR reversible Marmottant",
        case=values,
        apecss_commit=APECSS_COMMIT,
        binary=executable_sha256,
        analysis_schema=1,
    )
    return hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def compress_raw(attempt):
    """Losslessly compress ONLY new cache raw files; byte hash verifies recovery."""
    records = []
    for raw in Path(attempt).glob("RPAR_*.txt"):
        content = raw.read_bytes()
        target = raw.with_suffix(raw.suffix + ".gz")
        tmp = target.with_suffix(".tmp")
        with gzip.open(tmp, "wb", compresslevel=1) as stream:
            stream.write(content)
        if gzip.decompress(tmp.read_bytes()) != content:
            raise RuntimeError("Lossless raw compression verification failed")
        tmp.replace(target)
        records.append(
            dict(
                path=target.name,
                original_sha256=hashlib.sha256(content).hexdigest(),
                original_bytes=len(content),
            )
        )
        raw.unlink()
    if records:
        durable_json(Path(attempt) / "RAW_COMPRESSION.json", records)
    return records


class ResponseNodeCache:
    """flock serializes identical keys; completed attempts recover after a crash.

    Failed nodes remain explicit and cached. A partial attempt can be resumed by
    a fresh attempt, while a successfully finished APECSS case is never repeated.
    """

    def __init__(self, root, executable):
        self.root, self.executable = Path(root), Path(executable).resolve()
        self.binary_hash = hashlib.sha256(self.executable.read_bytes()).hexdigest()
        self.root.mkdir(parents=True, exist_ok=True)

    def get(self, case):
        key = node_key(case, self.binary_hash)
        folder = self.root / key
        folder.mkdir(exist_ok=True)
        with (folder / "LOCK").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            record = folder / "NODE.json"
            if record.exists():
                meta = json.loads(record.read_text())
                if meta["status"] == "COMPLETE":
                    payload = folder / meta["payload"]
                    if (
                        not payload.is_file()
                        or hashlib.sha256(payload.read_bytes()).hexdigest() != meta["payload_sha256"]
                    ):
                        raise ValueError("Incomplete/corrupt table cache payload; repair explicitly")
                return dict(meta, cache_hit=True)
            # Recover finished upstream output if interrupted before atomic publication.
            attempts = sorted(folder.glob("attempt_*"))
            recovered = None
            for attempt in attempts:
                result = attempt / "RESULT.json"
                if result.is_file():
                    m = json.loads(result.read_text())
                    if m.get("status") == "COMPLETE" and (attempt / "cycles.npz").is_file():
                        recovered = (attempt, m)
                        break
            if recovered is None:
                attempt = folder / f"attempt_{len(attempts):04d}"
                try:
                    meta, _, _ = ApecssBubbleResponseModel(self.executable).solve(case, attempt)
                except Exception as error:
                    meta = dict(status="FAILED", error=repr(error), key=key, case=asdict(case))
                    durable_json(record, meta)
                    return dict(meta, cache_hit=False)
            else:
                attempt, meta = recovered
            payload = attempt / "cycles.npz"
            with np.load(payload) as data:
                h = data["harmonics_cosine_m3"]
                r = data["radius_m"]
            rb = case.radius_m / np.sqrt(1 + case.sigma_initial_N_m / case.elasticity_N_m)
            rr = rb * np.sqrt(1 + case.sigma_water_N_m / case.elasticity_N_m)
            metrics = dict(
                meta["metrics"],
                cycles=case.cycles,
                transient_cycles_discarded=case.cycles - 8,
                entered_buckling_retained_cycles=bool(np.any(r < rb)),
                entered_rupture_retained_cycles=bool(np.any(r > rr)),
                solver_termination=meta["returncode"],
            )
            compress_raw(attempt)
            meta = dict(
                status="COMPLETE",
                key=key,
                case=asdict(case),
                metrics=metrics,
                harmonics_real=h.real.tolist(),
                harmonics_imag=h.imag.tolist(),
                payload=str(payload.relative_to(folder)),
                payload_sha256=hashlib.sha256(payload.read_bytes()).hexdigest(),
                recovered_completed_attempt=recovered is not None,
                elapsed_s=meta["elapsed_s"],
                apecss_commit=APECSS_COMMIT,
            )
            durable_json(record, meta)
            return dict(meta, cache_hit=recovered is not None)


def cached_task(task):
    cache, exe, case = task
    return ResponseNodeCache(cache, exe).get(case)


def harmonics(node):
    return np.array(node["harmonics_real"]) + 1j * np.array(node["harmonics_imag"])


def interpolate(weights, corner_values):
    w = np.prod(np.where(CORNERS, weights, 1 - weights), axis=1)
    return np.tensordot(w, corner_values, axes=1)


def response_errors(pred, actual, f, max_radius_pred, min_radius_pred):
    h = harmonics(actual)
    m = actual["metrics"]
    scale = max(abs(h[0]), 1e-30)
    return dict(
        complex_fundamental=abs(pred[0] - h[0]) / scale,
        amplitude=abs(abs(pred[0]) - abs(h[0])) / scale,
        phase_rad=abs(np.angle(pred[0] * np.conj(h[0]))),
        max_R_ratio=abs(max_radius_pred - m["max_R_over_R0"]) / m["max_R_over_R0"],
        min_R_ratio=abs(min_radius_pred - m["min_R_over_R0"]) / m["min_R_over_R0"],
        derivative_weighted_harmonics=np.linalg.norm(np.arange(1, len(h) + 1) * (pred - h))
        / max(np.linalg.norm(np.arange(1, len(h) + 1) * h), 1e-30),
        dot_volume_energy=abs(
            0.5 * np.sum(abs(2 * np.pi * f * np.arange(1, len(h) + 1) * pred) ** 2)
            - m["mean_dot_volume_squared_m6_s2"]
        )
        / max(m["mean_dot_volume_squared_m6_s2"], 1e-60),
    )


class AdaptiveResponseTableBuilder:
    def __init__(
        self,
        directory,
        executable,
        parameter_file,
        *,
        profile="production",
        max_depth=16,
        max_nodes=12000,
        workers=6,
        bounds=((1.6e-6, 1.96e-6), (2.3e6, 2.7e6), (5e3, 40e3)),
        cache=None,
    ):
        if profile not in PROFILES:
            raise ValueError("Unknown tolerance profile")
        if workers < 1 or max_nodes < 17 or max_depth < 0:
            raise ValueError("Positive worker/node budget and nonnegative depth required")
        self.out = Path(directory)
        self.out.mkdir(parents=True, exist_ok=True)
        self.exe = str(Path(executable).resolve())
        self.cache = str(Path(cache or self.out / "cache").resolve())
        self.parameters = SonoVueParameterDistribution(parameter_file)
        self.bounds = np.array(bounds, float)
        if (
            self.bounds.shape != (3, 2)
            or not np.isfinite(self.bounds).all()
            or np.any(self.bounds[:, 1] <= self.bounds[:, 0])
            or np.any(self.bounds < 0)
        ):
            raise ValueError("Invalid adaptive domain")
        self.config = dict(
            profile=profile,
            tolerance=PROFILES[profile],
            max_depth=max_depth,
            max_nodes=max_nodes,
            bounds=self.bounds.tolist(),
            parameter_sha256=hashlib.sha256(Path(parameter_file).read_bytes()).hexdigest(),
            sigma_initial_N_m=0.036,
            apecss_commit=APECSS_COMMIT,
            executable_sha256=hashlib.sha256(Path(executable).read_bytes()).hexdigest(),
            checks=CHECKS.tolist(),
            refinement="measured_axis_curvature_v1",
        )
        self.workers = workers
        self.nodes = {}
        self.leaves = []
        self.queue = []
        checkpoint = self.out / "CHECKPOINT.json"
        if checkpoint.exists():
            data = json.loads(checkpoint.read_text())
            if data["config"] != self.config:
                raise ValueError("Resume configuration mismatch")
            self.nodes, self.leaves, self.queue = data["nodes"], data["leaves"], data["queue"]
        else:
            self.queue = [dict(lo=self.bounds[:, 0].tolist(), hi=self.bounds[:, 1].tolist(), depth=0)]
        self.cache_hits = 0
        self.new_calls = 0

    def case(self, point):
        r, f, p = point
        s = self.parameters.parameter_set(r)
        return RadialCase(
            r,
            f,
            p,
            s.shell_elasticity_N_m.value,
            s.shell_viscosity_kg_s.value,
            0.036,
            provenance=json.dumps(s.to_dict(), sort_keys=True) + "; sigma_init .036 hypothesis",
        )

    @staticmethod
    def point_key(point):
        return ",".join(float(x).hex() for x in point)

    def fetch(self, points, pool):
        unique = {self.point_key(p): p for p in points if self.point_key(p) not in self.nodes}
        if len(self.nodes) + len(unique) > self.config["max_nodes"]:
            return False
        tasks = [(self.cache, self.exe, self.case(p)) for p in unique.values()]
        for key, node in zip(unique, pool.map(cached_task, tasks)):
            self.nodes[key] = node
            self.cache_hits += node["cache_hit"]
            self.new_calls += not node["cache_hit"]
        return True

    def checkpoint(self):
        durable_json(
            self.out / "CHECKPOINT.json",
            dict(config=self.config, nodes=self.nodes, leaves=self.leaves, queue=self.queue),
        )

    def build(self):
        start = time.perf_counter()
        with ProcessPoolExecutor(max_workers=self.workers) as pool:
            while self.queue:
                # Evaluate a breadth-first batch to keep the pool supplied.
                batch = self.queue[:32]
                self.queue = self.queue[32:]
                for cell in batch:
                    lo, hi = np.array(cell["lo"]), np.array(cell["hi"])
                    points = lo + np.vstack([CORNERS, CHECKS]) * (hi - lo)
                    if not self.fetch(points, pool):
                        cell.update(status="INVALID_CELL", reason="node_budget")
                        self.leaves.append(cell)
                        continue
                    nodes = [self.nodes[self.point_key(p)] for p in points]
                    if any(
                        n["status"] != "COMPLETE"
                        or not n["metrics"]["valid"]
                        or not n["metrics"]["harmonic_energy_resolved"]
                        for n in nodes
                    ):
                        cell.update(status="INVALID_CELL", reason="physical_or_numerical_node_failure")
                        self.leaves.append(cell)
                        continue
                    corners = np.array([harmonics(n) for n in nodes[:8]])
                    maxr = np.array([n["metrics"]["max_R_over_R0"] for n in nodes[:8]])
                    minr = np.array([n["metrics"]["min_R_over_R0"] for n in nodes[:8]])
                    errors = []
                    for w, point, node in zip(CHECKS, points[8:], nodes[8:]):
                        errors.append(
                            response_errors(
                                interpolate(w, corners),
                                node,
                                point[1],
                                interpolate(w, maxr),
                                interpolate(w, minr),
                            )
                        )
                    worst = {k: float(max(e[k] for e in errors)) for k in errors[0]}
                    cell["errors"] = worst
                    if max(worst.values()) <= self.config["tolerance"]:
                        cell.update(status="VALID_CELL", reason="nine_interior_checks")
                        self.leaves.append(cell)
                    elif cell["depth"] >= self.config["max_depth"]:
                        cell.update(status="INVALID_CELL", reason="max_depth")
                        self.leaves.append(cell)
                    else:
                        center = harmonics(nodes[8])
                        scores = []
                        for axis_index in range(3):
                            hlow = harmonics(nodes[8 + 3 + 2 * axis_index])
                            hhigh = harmonics(nodes[8 + 4 + 2 * axis_index])
                            scores.append(
                                float(np.linalg.norm(np.arange(1, 9) * ((hlow + hhigh) / 2 - center)))
                            )
                        axis = (
                            int(np.argmax(scores))
                            if max(scores) > 1e-30
                            else int(np.argmax((hi - lo) / (self.bounds[:, 1] - self.bounds[:, 0])))
                        )
                        mid = (lo[axis] + hi[axis]) / 2
                        for side in (0, 1):
                            aa, bb = lo.copy(), hi.copy()
                            if side:
                                aa[axis] = mid
                            else:
                                bb[axis] = mid
                            self.queue.append(dict(lo=aa.tolist(), hi=bb.tolist(), depth=cell["depth"] + 1))
                self.checkpoint()
                print(
                    f"leaves={len(self.leaves)} queued={len(self.queue)} nodes={len(self.nodes)}", flush=True
                )
        self.export()
        total = float(np.prod(self.bounds[:, 1] - self.bounds[:, 0]))
        coverage = (
            sum(np.prod(np.array(c["hi"]) - c["lo"]) for c in self.leaves if c["status"] == "VALID_CELL")
            / total
        )
        return dict(
            status="PASS_EMPIRICAL_DOMAIN" if coverage > 1 - 1e-12 else "PARTIAL_DOMAIN",
            coverage_fraction=float(coverage),
            nodes=len(self.nodes),
            leaves=len(self.leaves),
            valid_leaves=sum(c["status"] == "VALID_CELL" for c in self.leaves),
            new_apecss_calls=self.new_calls,
            cache_hits=self.cache_hits,
            runtime_s=time.perf_counter() - start,
            profile=self.config["profile"],
            workers=self.workers,
        )

    def export(self):
        temp = self.out / "adaptive_response.h5.tmp"
        n = len(self.leaves)
        h = np.zeros((n, 8, 8), complex)
        maxr = np.zeros((n, 8))
        minr = np.zeros((n, 8))
        for i, c in enumerate(self.leaves):
            if c["status"] != "VALID_CELL":
                continue
            for j, p in enumerate(np.array(c["lo"]) + CORNERS * (np.array(c["hi"]) - c["lo"])):
                node = self.nodes[self.point_key(p)]
                h[i, j] = harmonics(node)
                maxr[i, j] = node["metrics"]["max_R_over_R0"]
                minr[i, j] = node["metrics"]["min_R_over_R0"]
        with h5py.File(temp, "w") as f:
            f.attrs.update(
                schema_version=2,
                complete=True,
                metadata_json=json.dumps(
                    dict(
                        self.config,
                        medium={"density_kg_m3": 1000.0, "sound_speed_m_s": 1500.0, "viscosity_Pa_s": 0.002},
                        physics_status="partially validated",
                        parameter_calibration_validated=False,
                    )
                ),
            )
            for k, v in dict(
                lower=[c["lo"] for c in self.leaves],
                upper=[c["hi"] for c in self.leaves],
                valid=[c["status"] == "VALID_CELL" for c in self.leaves],
                harmonics=h,
                max_R=maxr,
                min_R=minr,
            ).items():
                f[k] = v
            f.flush()
        with temp.open("rb") as stream:
            os.fsync(stream.fileno())
        temp.replace(self.out / "adaptive_response.h5")


class AdaptiveResponseTableBubbleModel(ResponseTableBubbleModel):
    def __init__(self, path, *, backend="numpy"):
        if backend == "cupy":
            import cupy as xp
        elif backend == "numpy":
            xp = np
        else:
            raise ValueError("Unknown table backend")
        self.xp = xp
        self.path = Path(path)
        if (self.path.parent / "TABLE_REJECTION.json").exists():
            raise ValueError("Adaptive table explicitly rejected by independent validation")
        with h5py.File(path) as h:
            if h.attrs.get("schema_version") != 2 or not h.attrs.get("complete", False):
                raise ValueError("Incomplete adaptive table")
            self.metadata = json.loads(h.attrs["metadata_json"])
            self.lower = xp.asarray(h["lower"][:])
            self.upper = xp.asarray(h["upper"][:])
            self.valid = xp.asarray(h["valid"][:], bool)
            self.harmonics = xp.asarray(h["harmonics"][:])
            self.max_R = xp.asarray(h["max_R"][:])
            self.min_R = xp.asarray(h["min_R"][:])
        if (
            self.lower.shape != self.upper.shape
            or self.lower.shape[1:] != (3,)
            or self.harmonics.shape != (len(self.lower), 8, 8)
            or self.valid.shape != (len(self.lower),)
            or self.max_R.shape != (len(self.lower), 8)
            or self.min_R.shape != (len(self.lower), 8)
        ):
            raise ValueError("Invalid adaptive table shape")
        if (
            bool(xp.any(~xp.isfinite(self.harmonics[self.valid])))
            or bool(xp.any(~xp.isfinite(self.lower)))
            or bool(xp.any(~xp.isfinite(self.upper)))
            or bool(xp.any(self.upper <= self.lower))
            or bool(xp.any(~xp.isfinite(self.max_R[self.valid])))
            or bool(xp.any(~xp.isfinite(self.min_R[self.valid])))
            or bool(xp.any(self.min_R[self.valid] <= 0))
            or bool(xp.any(self.max_R[self.valid] < self.min_R[self.valid]))
        ):
            raise ValueError("Invalid adaptive cell values")
        self.axes = tuple(xp.asarray(v) for v in self.metadata["bounds"])

    def query(self, R0_m, frequency_Hz, pressure_amplitude_Pa):
        xp = self.xp
        q = xp.broadcast_arrays(
            xp.asarray(R0_m, float), xp.asarray(frequency_Hz, float), xp.asarray(pressure_amplitude_Pa, float)
        )
        shape = q[0].shape
        points = xp.stack([x.ravel() for x in q], axis=1)
        if bool(xp.any(~xp.isfinite(points))):
            raise ResponseTableBoundsError("Nonfinite adaptive query")
        selected = xp.full(len(points), -1, dtype=int)
        # Chunk leaves and queries: bounded memory and vector kernels on GPU.
        for begin in range(0, len(points), 4096):
            block = points[begin : begin + 4096, None, :]
            chosen = selected[begin : begin + 4096]
            for start in range(0, len(self.lower), 128):
                stop = min(start + 128, len(self.lower))
                inside = (
                    xp.all(
                        (block >= self.lower[None, start:stop]) & (block <= self.upper[None, start:stop]),
                        axis=2,
                    )
                    & self.valid[None, start:stop]
                )
                chosen = xp.maximum(
                    chosen, xp.max(xp.where(inside, xp.arange(start, stop)[None, :], -1), axis=1)
                )
            selected[begin : begin + 4096] = chosen
        if bool(xp.any(selected < 0)):
            raise ResponseTableBoundsError("Outside domain or INVALID_CELL; no extrapolation")
        weights = (points - self.lower[selected]) / (self.upper[selected] - self.lower[selected])
        coeff = xp.prod(
            xp.where(xp.asarray(CORNERS)[None, :, :], weights[:, None, :], 1 - weights[:, None, :]), axis=2
        )
        values = xp.sum(coeff[:, :, None] * self.harmonics[selected], axis=1).reshape(shape + (8,))
        return dict(
            harmonics_m3=values,
            max_R_over_R0=xp.sum(coeff * self.max_R[selected], axis=1).reshape(shape),
            min_R_over_R0=xp.sum(coeff * self.min_R[selected], axis=1).reshape(shape),
            corner_envelope_max_R_over_R0=xp.max(self.max_R[selected], axis=1).reshape(shape),
            corner_envelope_min_R_over_R0=xp.min(self.min_R[selected], axis=1).reshape(shape),
            harmonic_resolved=xp.ones(shape, bool),
            within_response_table=xp.ones(shape, bool),
            interpolation_validated=xp.ones(shape, bool),
            parameter_calibration_validated=False,
        )
