#!/usr/bin/env python3
"""Continuous number-weighted SonoVue diameters from the frozen histogram.

This is not a parametric fit: each original bin has a constant PDF, hence
a linear CDF. Empirical bin masses, zero bins and source edges are retained.
One vector of U ~ Uniform[0, 1) is inverted analytically. There is no second
random draw, smoothing, fitted distribution family or mean matching.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np

SAMPLER_VERSION = "SONOVUE_CONTINUOUS_EMPIRICAL_INVERSE_CDF_V0"
DEFAULT_HISTOGRAM = Path(__file__).resolve().parents[1] / "input/FROZEN_SONOVUE_HISTOGRAM.csv"
REQUIRED_COLUMNS = ("diameter_center_um", "diameter_low_um", "diameter_high_um", "probability")
EXPECTED_BIN_WIDTH_UM = 0.1  # Source-defined bin width; support is never hard-coded.
WIDTH_ATOL_UM = 1e-10
FLOAT_TOLERANCE = 32 * np.finfo(np.float64).eps


class HistogramError(ValueError):
    """Invalid source histogram; never silently repair scientific inputs."""


def sha256_file(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _integer(value: int, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an explicitly supplied non-negative integer")
    if value < 0:
        raise ValueError(f"{name} must be non-negative")
    return int(value)


class SonoVueDistribution:
    """Read-only bin table with vectorized CDF, PDF and analytic inverse.

    CDF plateaus use a right-side inverse convention: for 0 <= u < 1,
    select the first positive-mass interval with cdf_high > u. Thus an
    exact plateau probability maps to the next positive bin's low edge.
    inverse(0) is the first positive bin's low edge; inverse(1) is the last
    positive bin's high edge. Uniform sampling never draws 1.

    Samples belong to half-open positive intervals [low, high). If floating
    arithmetic rounds a u<1 result to high, move it to nextafter(high, low).
    This one-ULP boundary guard prevents assignment to a neighboring zero
    bin; it does not truncate the distribution's support or remove tails.
    Gaps between original bins, if present, have zero density.
    """

    def __init__(self, histogram: str | Path = DEFAULT_HISTOGRAM):
        self.histogram_path = Path(histogram).resolve()
        raw = self.histogram_path.read_bytes()  # FileNotFoundError remains explicit.
        self.histogram_sha256 = hashlib.sha256(raw).hexdigest()
        try:
            reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
            if reader.fieldnames is None or not set(REQUIRED_COLUMNS) <= set(reader.fieldnames):
                raise HistogramError(f"Required CSV columns: {', '.join(REQUIRED_COLUMNS)}")
            if len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise HistogramError("Duplicate CSV column names")
            rows = list(reader)
            if not rows or any(None in row for row in rows):
                raise HistogramError("Empty or malformed CSV rows")
            table = np.array([[float(row[name]) for name in REQUIRED_COLUMNS] for row in rows], dtype=np.float64)
        except (UnicodeError, csv.Error, TypeError, KeyError, ValueError) as exc:
            raise HistogramError(f"CSV parse failed: {exc}") from exc
        if not np.isfinite(table).all():
            raise HistogramError("Required numeric values must be finite")
        self.center, self.low, self.high, raw_probability = table.T.copy()
        self.width = self.high - self.low
        if np.any(self.low <= 0) or np.any(self.width <= 0):
            raise HistogramError("Diameter edges must be positive and every bin width > 0")
        if np.any(np.diff(self.center) <= 0) or np.any(np.diff(self.low) <= 0):
            raise HistogramError("Bins must already be ordered; automatic sorting is forbidden")
        if np.any(self.low[1:] < self.high[:-1]):
            raise HistogramError("Overlapping source intervals are forbidden")
        if not np.allclose(self.width, EXPECTED_BIN_WIDTH_UM, rtol=0, atol=WIDTH_ATOL_UM):
            raise HistogramError("Original bin widths must be approximately 0.1 um")
        if not np.allclose(self.center, (self.low + self.high) / 2, rtol=0, atol=WIDTH_ATOL_UM):
            raise HistogramError("Diameter centers must agree with their source edges")
        if np.any(raw_probability < 0):
            raise HistogramError("Negative probability is forbidden")
        try:
            self.input_probability_sum = math.fsum(raw_probability)
        except OverflowError as exc:
            raise HistogramError("Probability total overflows") from exc
        if not math.isfinite(self.input_probability_sum) or self.input_probability_sum <= 0:
            raise HistogramError("Probability sum must be finite and positive")
        # Only one common normalization factor; no per-bin reweighting or tail edits.
        self.probability = raw_probability / self.input_probability_sum
        self.normalized_probability_sum = math.fsum(self.probability)
        if np.any((raw_probability > 0) & (self.probability == 0)):
            raise HistogramError("Normalization underflows a positive bin; no tail mass may be dropped")
        if abs(self.normalized_probability_sum - 1) > FLOAT_TOLERANCE:
            raise HistogramError("Normalized probability sum fails floating-point tolerance")
        self.positive = np.flatnonzero(self.probability > 0)
        self.cdf_high = np.cumsum(self.probability, dtype=np.float64)
        # Exact terminal value, including any trailing zero bins. Only FP closure.
        self.cdf_high[self.positive[-1]:] = 1.0
        self.cdf_low = np.concatenate(([0.0], self.cdf_high[:-1]))
        self.cdf_mass = self.cdf_high - self.cdf_low
        if np.any(self.cdf_mass[self.positive] <= 0):
            raise HistogramError("A positive bin is unresolvable in float64 cumulative probability")
        if np.any(self.cdf_mass[self.probability == 0] != 0):
            raise HistogramError("Float64 closure would assign mass to a zero bin")
        if not np.all((self.cdf_high >= 0) & (self.cdf_high <= 1)):
            raise HistogramError("Cumulative probabilities outside [0,1]")
        self.pdf_per_um = self.probability / self.width
        self.support_min_um = float(self.low[0])
        self.support_max_um = float(self.high[-1])
        self.bin_width_um = float(np.median(self.width))
        self.cdf_mass_roundoff_max = float(np.max(np.abs(self.cdf_mass - self.probability)))
        if self.cdf_mass_roundoff_max > FLOAT_TOLERANCE:
            raise HistogramError("Cumulative probability closure exceeds fixed floating tolerance")
        for name in ("center", "low", "high", "width", "probability", "positive", "cdf_low", "cdf_high", "cdf_mass", "pdf_per_um"):
            getattr(self, name).flags.writeable = False

    @staticmethod
    def _query(values):
        x = np.asarray(values, dtype=np.float64)
        if np.isnan(x).any():
            raise ValueError("NaN is not a valid diameter query")
        return x

    def cdf(self, diameter_um):
        x = self._query(diameter_um)
        idx = np.clip(np.searchsorted(self.low, x, side="right") - 1, 0, len(self.low) - 1)
        fraction = np.clip((x - self.low[idx]) / self.width[idx], 0, 1)
        value = self.cdf_low[idx] + fraction * self.cdf_mass[idx]
        value = np.where(x < self.low[0], 0.0, value)
        value = np.where(x >= self.high[-1], 1.0, value)
        return float(value) if x.ndim == 0 else value

    def pdf(self, diameter_um):
        x = self._query(diameter_um)
        idx = np.clip(np.searchsorted(self.low, x, side="right") - 1, 0, len(self.low) - 1)
        inside = (x >= self.low[idx]) & (x < self.high[idx])
        value = np.where(inside, self.pdf_per_um[idx], 0.0)
        return float(value) if x.ndim == 0 else value

    def _inverse_with_indices(self, u_values):
        u = np.asarray(u_values, dtype=np.float64)
        if not np.isfinite(u).all() or np.any((u < 0) | (u > 1)):
            raise ValueError("Inverse-CDF probabilities must be finite and in [0,1]")
        slot = np.searchsorted(self.cdf_high[self.positive], u, side="right")
        idx = self.positive[np.minimum(slot, len(self.positive) - 1)]
        fraction = (u - self.cdf_low[idx]) / self.cdf_mass[idx]
        diameter = self.low[idx] + fraction * self.width[idx]
        upper = np.where(u == 1, self.high[idx], np.nextafter(self.high[idx], self.low[idx]))
        diameter = np.maximum(self.low[idx], np.minimum(diameter, upper))
        return diameter, idx

    def inverse_cdf(self, u_values):
        diameter, _ = self._inverse_with_indices(u_values)
        return float(diameter) if diameter.ndim == 0 else diameter

    def sample_with_provenance(self, n: int, seed: int):
        n, seed = _integer(n, "n"), _integer(seed, "seed")
        rng = np.random.default_rng(seed)
        u = rng.random(n, dtype=np.float64)
        diameter, idx = self._inverse_with_indices(u)
        return diameter, u, idx

    def sample_diameters(self, n: int, seed: int) -> np.ndarray:
        return self.sample_with_provenance(n, seed)[0]

    def target_mean_um(self) -> float:
        return float(np.sum(self.probability * (self.low + self.high) / 2))

    def export_cdf(self, output: str | Path) -> None:
        """Keep both edges of every bin, including duplicate shared edges and zeros."""
        output = Path(output)
        if output.exists():
            raise FileExistsError(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        rows = np.empty((2 * len(self.low), 9), dtype=np.float64)
        rows[0::2, 0], rows[1::2, 0] = self.low, self.high
        rows[0::2, 1], rows[1::2, 1] = self.cdf_low, self.cdf_high
        rows[:, 2:] = np.repeat(np.column_stack((np.arange(len(self.low)), self.low, self.high, self.probability, self.pdf_per_um, self.cdf_low, self.cdf_high)), 2, axis=0)
        np.savetxt(output, rows, delimiter=",", fmt=["%.17g", "%.17g", "%d"] + ["%.17g"] * 6,
                   header="diameter_um,cdf,bin_index,bin_low_um,bin_high_um,probability_mass,pdf_per_um,cdf_low,cdf_high", comments="")


def sample_diameters(n: int, seed: int, histogram: str | Path = DEFAULT_HISTOGRAM) -> np.ndarray:
    return SonoVueDistribution(histogram).sample_diameters(n, seed)


def write_population(distribution: SonoVueDistribution, n: int, seed: int, output: str | Path,
                     role: str = "SAMPLER_DEMO_ONLY") -> dict:
    if role not in ("SAMPLER_DEMO_ONLY", "SAMPLER_VALIDATION_ONLY"):
        raise ValueError("V0 CLI does not issue formal simulation populations")
    n, seed = _integer(n, "n"), _integer(seed, "seed")
    output = Path(output)
    metadata_path = output.with_suffix(".metadata.json")
    if output.exists() or metadata_path.exists():
        raise FileExistsError("Output/metadata already exists; immutable cases are not overwritten")
    output.parent.mkdir(parents=True, exist_ok=True)
    diameter, u, idx = distribution.sample_with_provenance(n, seed)
    # Sampling is vectorized; np.savetxt only serializes the already generated population.
    data = np.column_stack((np.arange(n, dtype=np.int64), diameter, distribution.low[idx], distribution.high[idx], u))
    np.savetxt(output, data, delimiter=",", fmt=["%d"] + ["%.17g"] * 4,
               header="bubble_id,diameter_um,source_bin_low_um,source_bin_high_um,u_sample", comments="")
    metadata = {
        "schema": "SONOVUE_POPULATION_PROVENANCE_V0", "role": role,
        "formal_simulation_population": False, "N": n, "seed": seed,
        "seed_policy": "FROZEN_PER_CASE", "histogram_file": distribution.histogram_path.name,
        "histogram_sha256": distribution.histogram_sha256, "sampler_version": SAMPLER_VERSION,
        "sampler_source_sha256": sha256_file(__file__), "random_generator": "numpy.random.default_rng",
        "bit_generator": type(np.random.default_rng(seed).bit_generator).__name__,
        "numpy_version": np.__version__, "python_version": platform.python_version(),
        "platform": platform.platform(), "dtype": "float64", "byte_order": sys.byteorder,
        "sampling_method": "inverse-CDF", "uniform_domain": "[0,1)",
        "support_min_um": distribution.support_min_um, "support_max_um": distribution.support_max_um,
        "diameter_sequence_sha256_float64_le": hashlib.sha256(diameter.astype("<f8", copy=False).tobytes()).hexdigest(),
        "population_file": output.name, "population_sha256": sha256_file(output),
        "source_bin_columns_role": "Provenance of analytic inverse interval; not a second random sampling stage",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histogram", type=Path, default=DEFAULT_HISTOGRAM)
    parser.add_argument("--n", required=True, type=int)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--role", choices=("SAMPLER_DEMO_ONLY", "SAMPLER_VALIDATION_ONLY"), default="SAMPLER_DEMO_ONLY")
    args = parser.parse_args(argv)
    if not args.histogram.is_file():
        print("STATUS = BLOCKED_INPUT_MISSING", file=sys.stderr)
        return 2
    try:
        distribution = SonoVueDistribution(args.histogram)
        meta = write_population(distribution, args.n, args.seed, args.output, args.role)
    except (ValueError, OSError) as exc:
        print(f"STATUS = FAIL; {exc}", file=sys.stderr)
        return 2
    print(json.dumps({"status": "PASS", "N": meta["N"], "seed": meta["seed"], "role": meta["role"],
                      "histogram_sha256": meta["histogram_sha256"], "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
