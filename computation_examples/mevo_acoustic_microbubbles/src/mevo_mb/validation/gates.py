"""Downstream gates require named, hashed evidence and matching configuration scope."""

import hashlib
import json
from pathlib import Path


class ValidationGateError(RuntimeError):
    def __init__(self, metadata):
        self.metadata = metadata
        super().__init__("NOT_RUN_GATE: " + ", ".join(metadata["failed_gates"]))


BASIC_GATES = ["ACOUSTIC_PASS", "BUBBLE_RESPONSE_PASS", "PRIMARY_PASS", "SECONDARY_BASIC_PASS"]
REQUIRES = {
    "run-lazarus": BASIC_GATES + ["CLUSTERING_TREND_PASS"],
    "run-sonovue-clustering": BASIC_GATES + ["CLUSTERING_TREND_PASS"],
    "run-cylinder": BASIC_GATES + ["CLUSTERING_TREND_PASS", "FLOW_COUPLING_PASS"],
    "run-frequency-sweep": BASIC_GATES + ["CLUSTERING_TREND_PASS", "FLOW_COUPLING_PASS", "CYLINDER_PASS"],
}


class ValidationGateManager:
    def __init__(self, path):
        self.path = Path(path)
        self.registry = json.loads(self.path.read_text()) if self.path.exists() else {"gates": {}}

    def check(self, command, *, scope, force_demo=False):
        if command not in REQUIRES:
            raise ValueError("Unknown gated command")
        failed = []
        for name in REQUIRES[command]:
            gate = self.registry.get("gates", {}).get(name, {})
            valid = gate.get("status") == "PASS" and scope in gate.get("scopes", [])
            evidence = gate.get("evidence", [])
            valid &= bool(evidence)
            for item in evidence:
                p = Path(item["path"])
                valid &= p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest() == item["sha256"]
            if not valid:
                failed.append(name)
        meta = dict(
            command=command,
            scope=scope,
            failed_gates=failed,
            UNVALIDATED_OVERRIDE=bool(failed and force_demo),
            gate_status="OVERRIDDEN" if failed and force_demo else ("NOT_RUN_GATE" if failed else "PASS"),
        )
        # Admission is not a certificate for the subsequent trajectory. Its
        # measured validity and figure watermark are assigned after simulation.
        if failed and force_demo:
            meta.update(physics_status="demo only", validation_watermark="DEMO / NOT VALIDATED")
        if failed and not force_demo:
            raise ValidationGateError(meta)
        return meta


def validation_watermark(metadata):
    if metadata.get("UNVALIDATED_OVERRIDE"):
        return "DEMO / NOT VALIDATED"
    if metadata.get("paper_reproduction"):
        return "PAPER REPRODUCTION"
    if metadata.get("physics_status") == "validated":
        return "VALIDATED AGAINST " + metadata.get("reference", "DECLARED GATES")
    if metadata.get("physics_status") == "partially validated":
        return "PARTIALLY VALIDATED"
    return "NOT PHYSICALLY VALIDATED"
