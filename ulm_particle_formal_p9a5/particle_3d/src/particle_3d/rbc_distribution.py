"""C57BL/6 V0 geometry distribution; explicit seed, paired candidates, no clipping.

Latent D/V are independent. Guard and shape conditioning change the accepted
distribution. This layer knows no flow, timestep or orientation integrator.
"""
from dataclasses import dataclass
import csv
import hashlib
import json
from pathlib import Path
import platform
import numpy as np

DEFAULT_CONTRACT = Path(__file__).resolve().parents[2] / "contracts/C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json"
DISTRIBUTION_ID = "C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0"
COLUMNS = ["rbc_id", "candidate_id", "D_um", "V_fL", "a_um", "b_um", "c_um",
           "a_m", "b_m", "c_m", "volume_m3", "r", "jeffery_lambda"]
DTYPE = np.dtype([(name, "<i8" if name in COLUMNS[:2] else "<f8") for name in COLUMNS])
LEDGER_DTYPE = np.dtype([("candidate_id", "<i8"), ("D_raw_um", "<f8"), ("V_raw_fL", "<f8"), ("status", "U16")])


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_digest(array):
    return hashlib.sha256(np.asarray(array, dtype="<f8").tobytes(order="C")).hexdigest()


def load_contract(path=DEFAULT_CONTRACT):
    contract = json.loads(Path(path).read_text())
    if contract["contract_name"] != DISTRIBUTION_ID or contract["version"] != "V0":
        raise ValueError("unsupported RBC distribution contract")
    return contract


def explicit_integer(value, name, minimum=0):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an explicit integer >= {minimum}")
    return int(value)


def classify_candidate(diameter_um, volume_fL, contract):
    """Ordered disjoint rejection; original D/V never altered."""
    d, v = float(diameter_um), float(volume_fL)
    low, high = contract["diameter"]["guard_um"]
    if not np.isfinite(d) or not low <= d <= high:
        return "D_GUARD", None
    low, high = contract["volume"]["guard_fL"]
    if not np.isfinite(v) or not low <= v <= high:
        return "V_GUARD", None
    a = d / 2
    c = 3 * v / (np.pi * d * d)
    if not np.isfinite(c) or not 0 < c < a:
        return "SHAPE", None
    r = c / a
    return "ACCEPT", (a, a, c, r, (r * r - 1) / (r * r + 1))


@dataclass(frozen=True)
class RBCPopulation:
    samples: np.ndarray
    candidates: np.ndarray
    metadata_json: str

    def __post_init__(self):
        for name in ["samples", "candidates"]:
            array = np.asarray(getattr(self, name))
            object.__setattr__(self, name, np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape))

    @property
    def metadata(self):
        return json.loads(self.metadata_json)


def sample_rbc_geometries(n, seed, contract=DEFAULT_CONTRACT):
    n = explicit_integer(n, "n", 1)
    seed = explicit_integer(seed, "seed")
    path = Path(contract)
    config = load_contract(path)
    rng = np.random.default_rng(seed)
    chunk_size = config["rng"]["chunk_size"]
    means = [config["diameter"]["mean_um"], config["volume"]["mean_fL"]]
    scales = [config["diameter"]["sd_um"], config["volume"]["sd_fL"]]
    accepted, candidates = [], []
    counts = {"D_GUARD": 0, "V_GUARD": 0, "SHAPE": 0}
    generated = 0
    while len(accepted) < n:
        raw = rng.normal(loc=means, scale=scales, size=(chunk_size, 2))
        generated += chunk_size
        for d, v in raw:
            candidate_id = len(candidates)
            status, derived = classify_candidate(d, v, config)
            candidates.append((candidate_id, d, v, status))
            if status != "ACCEPT":
                counts[status] += 1
                continue
            a, b, c, r, lam = derived
            accepted.append((len(accepted), candidate_id, d, v, a, b, c,
                             a * 1e-6, b * 1e-6, c * 1e-6, v * 1e-18, r, lam))
            if len(accepted) == n:
                break
    samples = np.array(accepted, dtype=DTYPE)
    ledger = np.array(candidates, dtype=LEDGER_DTYPE)
    metadata = dict(
        N=n, seed=seed, contract_name=DISTRIBUTION_ID, contract_version=config["version"],
        contract_sha256=digest(path), source_literature_ids=config["sources"],
        role="RBC_GEOMETRY_VALIDATION_ONLY", formal_hematocrit_population=False,
        production_particle_population=False, numpy_version=np.__version__,
        bit_generator=type(rng.bit_generator).__name__, python_version=platform.python_version(),
        platform=platform.platform(), chunk_size=chunk_size, generated_candidate_count=generated,
        candidate_count=len(candidates), unused_generated_tail_count=generated - len(candidates),
        diameter_guard_rejections=counts["D_GUARD"], volume_guard_rejections=counts["V_GUARD"],
        guard_rejection_count=counts["D_GUARD"] + counts["V_GUARD"],
        shape_rejection_count=counts["SHAPE"], acceptance_rate=n / len(candidates),
        count_semantics=config["rng"]["candidate_count_definition"],
        diameter_array_sha256=array_digest(samples["D_um"]), volume_array_sha256=array_digest(samples["V_fL"]),
        axes_array_sha256=array_digest(np.column_stack([samples[name] for name in ["a_m", "b_m", "c_m"]])),
        sample_structured_array_sha256=hashlib.sha256(samples.tobytes()).hexdigest(),
        array_hash_convention="float64 little-endian C order; axes rows [a_m,b_m,c_m]",
        cross_numpy_version_byte_identity_claimed=False)
    return RBCPopulation(samples, ledger, json.dumps(metadata, ensure_ascii=False, sort_keys=True))


def write_population(population, csv_path):
    """Stable UTF-8/LF CSV and adjacent CSV-bound provenance; no timestamps."""
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(COLUMNS)
        for row in population.samples:
            writer.writerow([str(int(row[n])) if n in COLUMNS[:2] else format(float(row[n]), ".17g") for n in COLUMNS])
    metadata = dict(population.metadata, csv_file=csv_path.name, csv_sha256=digest(csv_path))
    csv_path.with_suffix(".metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n")
    return metadata


def statistics(values):
    values = np.asarray(values, dtype=np.float64)
    quantiles = [1, 5, 25, 50, 75, 95, 99]
    return dict(N=len(values), mean=float(np.mean(values)), sd=float(np.std(values, ddof=1)),
                sd_definition="sample SD (ddof=1)", median=float(np.median(values)),
                quantiles={str(q): float(np.percentile(values, q)) for q in quantiles},
                role="FINAL_ACCEPTED_MODEL_DERIVED_NOT_DIRECT_LITERATURE_MEASUREMENT")


def quantile_indices(samples, fractions=(.05, .25, .5, .75, .95)):
    order = np.lexsort((samples["candidate_id"], samples["r"]))
    # Nearest observed order statistic, ties by source order; never invent axes.
    return np.array([order[int(np.floor(f * (len(order) - 1) + .5))] for f in fractions], dtype=np.int64)


def stratified_indices(samples, seed, n=64):
    seed = explicit_integer(seed, "selection seed")
    order = np.lexsort((samples["candidate_id"], samples["r"]))
    rng = np.random.default_rng(seed)
    return np.array([stratum[rng.integers(len(stratum))] for stratum in np.array_split(order, n)])
