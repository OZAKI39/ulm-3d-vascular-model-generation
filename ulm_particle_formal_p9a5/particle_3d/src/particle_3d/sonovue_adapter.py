"""Thin read-only adapter to the separately frozen SonoVue sampler; no copied sampler."""
from dataclasses import dataclass, asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import numpy as np
from .audit import check_hash
from .microbubble import positive_scalar


@dataclass(frozen=True)
class SingleMBSize:
    seed: int
    diameter_um: float
    diameter_m: float
    radius_m: float
    histogram_sha256: str
    sampler_version: str
    sampler_source_sha256: str
    diameter_sequence_sha256_float64_le: str
    numpy_version: str
    python_version: str
    N: int = 1
    formal_simulation_population: bool = False
    role: str = "DEMO / VALIDATION ONLY"

    def to_dict(self):
        return asdict(self)


def diameter_um_to_radius_m(diameter_um):
    return np.float64(0.5) * positive_scalar(diameter_um, "diameter_um") * np.float64(1e-6)


def read_sonovue(root):
    root = Path(root).resolve()
    receipt = {}
    for line in (root / "SHA256SUMS").read_text().splitlines():
        digest, relative = line.split("  ", 1)
        path = (root / relative).resolve()
        if not path.is_relative_to(root):
            raise ValueError("SonoVue manifest path escapes root")
        check_hash(path, digest)
        receipt[relative] = digest
    contract = json.loads((root / "contracts/SONOVUE_SAMPLER_CONTRACT_V0.json").read_text())
    if contract["status"] != "FROZEN_VALIDATED_PASS":
        raise ValueError("SonoVue contract is not frozen/validated")
    source = root / contract["sampler_source_file"]
    check_hash(source, contract["sampler_source_sha256"])
    check_hash(root / contract["histogram_file"], contract["histogram_sha256"])
    spec = importlib.util.spec_from_file_location("particle1_external_frozen_sonovue", source)
    module = importlib.util.module_from_spec(spec)
    # Executes the original module in memory; never writes a copied module/pyc.
    exec(compile(source.read_bytes(), str(source), "exec"), module.__dict__)
    if module.SAMPLER_VERSION != contract["sampler_version"]:
        raise ValueError("SonoVue version differs from frozen contract")
    distribution = module.SonoVueDistribution(root / contract["histogram_file"])
    return contract, distribution, receipt


def sample_single_validation_size(root, *, seed):
    contract, distribution, _ = read_sonovue(root)
    diameters = distribution.sample_diameters(n=1, seed=seed)
    diameter = float(diameters[0])
    return SingleMBSize(int(seed), diameter, diameter * 1e-6, float(diameter_um_to_radius_m(diameter)),
                        contract["histogram_sha256"], contract["sampler_version"], contract["sampler_source_sha256"],
                        hashlib.sha256(diameters.astype("<f8").tobytes()).hexdigest(), np.__version__, platform.python_version())
