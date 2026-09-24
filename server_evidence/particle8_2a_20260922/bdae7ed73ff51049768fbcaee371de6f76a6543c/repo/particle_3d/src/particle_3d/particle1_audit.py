"""P1 provenance gate: accepted P0 source remains byte-identical."""
import json
from pathlib import Path
from .audit import check_hash, sha256
from .sonovue_adapter import read_sonovue

P0_COMMIT = "7cfe5141382600e28582f05bff712a6f09c38a39"
FROZEN_COMMIT = "c83dccea0f1fe5cd9fd3f5939e2eba883aeba9c2"
P1_BRANCH = "dev/particle-1-single-mb-20260920"
VALIDATION_SEED = 20260920


def check_dependency(repo):
    repo = Path(repo)
    record = json.loads((repo / "particle_3d/reports/particle0/PARTICLE0_VALIDATION.json").read_text())
    if record["git_commit"] != P0_COMMIT or record["manual_visual_review"] != "PASS":
        raise ValueError("Particle-0 must match the user-accepted implementation")
    for name, digest in record["source_sha256"].items():
        check_hash(repo / name, digest)
    if record["tests"]["passed"] != 46 or record["tests"]["failed"] != 0:
        raise ValueError("P0 baseline test record changed")
    return record


def scope_record(repo, sonovue_root, frozen_summary):
    dependency = check_dependency(repo)
    contract, _, receipt = read_sonovue(sonovue_root)
    return dict(stage="Particle-1", particle0_dependency_commit=P0_COMMIT,
                particle0_manual_review=dependency["manual_review_evidence"], frozen_base_commit=FROZEN_COMMIT,
                branch=P1_BRANCH, mesh_sha256=frozen_summary["mesh_sha256"], flow_sha256=frozen_summary["flow_sha256"],
                particle_count=1, shape="sphere", core_units="SI", model="overdamped; no inertia",
                translation="V = u_inf", rotation="Omega = 0.5 * curl(u_inf)", position_update="explicit Euler; dt required",
                production_particle_timestep_frozen=False, formal_simulation_population=False,
                no_cfd_executed=True, particle2_started=False, particle0_enhanced_tests_added=False,
                sonovue_root=str(Path(sonovue_root).resolve()), sonovue_sampler_version=contract["sampler_version"],
                sonovue_histogram_sha256=contract["histogram_sha256"], sonovue_inventory=receipt,
                sonovue_sha256sums_sha256=sha256(Path(sonovue_root)/"SHA256SUMS"),
                excluded_physics=["RBC", "other bubbles", "inertia", "wall force", "lubrication", "contact", "adhesion", "Brownian", "gravity", "buoyancy", "acoustic force", "LAMMPS", "continuous injection"])
