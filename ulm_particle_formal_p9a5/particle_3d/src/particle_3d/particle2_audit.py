"""Particle-2 accepted dependency lock; old code, tests, images and data are read-only."""
import json
from pathlib import Path
from .audit import check_hash

P0_COMMIT="7cfe5141382600e28582f05bff712a6f09c38a39"
P1_COMMIT="6e59c605cb8251b9fa2e80d2dbed0cfa6daaf03c"
FROZEN_COMMIT="c83dccea0f1fe5cd9fd3f5939e2eba883aeba9c2"
P2_BRANCH="dev/particle-2-rbc-orientation-20260920"


def check_dependencies(repo):
    repo=Path(repo)
    records=[]
    for stage,commit in [("particle0",P0_COMMIT),("particle1",P1_COMMIT)]:
        record=json.loads((repo/f"particle_3d/reports/{stage}/{stage.upper()}_VALIDATION.json").read_text())
        if record["git_commit"]!=commit or record["manual_visual_review"]!="PASS":
            raise ValueError(f"{stage} must be the accepted implementation")
        for group in ["source_sha256","figure_sha256","data_sha256"]:
            for name,expected in record[group].items():check_hash(repo/name,expected)
        records.append(record)
    if records[1]["visual_step_jump_review"]!="PASS" or not records[1]["manual_review_evidence"]["particle2_authorized_to_start"]:
        raise ValueError("Particle-1 user acceptance required")
    return records
