"""Read-only, explicit input adaptation from the audited 3D project into SI."""
import math
from pathlib import Path


def resolve_source_path(value, root, historical_prefix):
    text = str(value).replace("\\", "/")
    root = Path(root).resolve()
    if text.startswith(historical_prefix):
        path = root / text[len(historical_prefix):]
    elif ":" in text:
        raise ValueError(f"Unknown historical path prefix: {value}")
    else:
        path = Path(text) if Path(text).is_absolute() else root/text
    path = path.resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"Source path escapes reference root: {value}")
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def length_to_m(value, units):
    if units not in ("m", "um"):
        raise ValueError(f"Unknown length unit: {units}")
    return value * {"m": 1.0, "um": 1e-6}[units]


def validate_contract(contract):
    if contract.get("status") != "PASS":
        raise ValueError("Source audit did not pass")
    if contract["geometry_units"] not in ("m", "um"):
        raise ValueError("Geometry units are not explicit")
    if contract["solver_units"] != "SI":
        raise ValueError("Solver must use SI")
    if not contract["inlets"] or not contract["outlets"]:
        raise ValueError("Both inlet and outlet metadata are required")
    ports = contract["inlets"] + contract["outlets"]
    ids = [p["port_id"] for p in ports]
    tags = [p["surface_entity_id"] for p in ports]
    if len(ids) != len(set(ids)) or len(tags) != len(set(tags)):
        raise ValueError("Port IDs and entity tags must be unique")
    for p in ports:
        if not math.isfinite(p["area_m2"]) or p["area_m2"] <= 0:
            raise ValueError("Port area must be finite and positive")
        normal = p["outward_normal"]
        if len(normal) != 3 or not all(math.isfinite(x) for x in normal) or sum(x*x for x in normal) <= 0:
            raise ValueError("Port normal must be a finite nonzero 3-vector")
        if not math.isclose(sum(x*x for x in normal), 1., rel_tol=1e-9):
            raise ValueError("Port normal must be unit length")
        if len(p["plane_origin_m"]) != 3 or not all(math.isfinite(x) for x in p["plane_origin_m"]):
            raise ValueError("Invalid plane origin")
    return True
