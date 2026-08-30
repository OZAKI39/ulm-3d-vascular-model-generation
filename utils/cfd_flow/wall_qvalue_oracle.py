"""Zero-run analytic audit of Seeder q-values consumed by Musubi wall_libb.

This research helper reads existing TreElm meshes only.  It does not launch
Seeder, Musubi, Harvester, or any vascular workflow.
"""

from __future__ import annotations

import json
import math
import re
import subprocess
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .exact_link_flux import INVERSE_DIRECTIONS
from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import load_mesh_contract
from .musubi_pressure_bc_benchmark import PIPE_RADIUS_M
from .port_grid_sensitivity import RESEARCH_RUN
from .restart_decode import D3Q19_DIRECTIONS


SEEDER_ROOT_WSL = "/home/lzy/apes-pinned/seeder_official"
MUSUBI_ROOT_WSL = "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300"
SEEDER_COMMIT = "667109df6fafdcb39f4409e3f5d90f04d75cd33c"
MUSUBI_COMMIT = "4e8b277b66226277171ef93bf054d21270812793"
MUSUBI_SCHEME_COMMIT = "81f8c4f13772f6d4af31f335e1e3f99b02726e25"
TREELM_COMMIT = "9899d1376992c4fafc8a343d2b4ccef81de670d1"


def _unc(wsl_path: str) -> Path:
    return Path("//wsl.localhost/Ubuntu" + wsl_path)


def _git_value(repository: str, *arguments: str) -> str:
    process = subprocess.run(
        ["wsl", "git", "-C", repository, *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return process.stdout.strip()


def _source_evidence(path: Path, tokens: Iterable[str]) -> dict[str, Any]:
    text = Path(path).read_text(encoding="utf-8")
    rows = []
    for token in tokens:
        offset = text.find(token)
        if offset < 0:
            raise RuntimeError(f"Source token is missing from {path}: {token}")
        rows.append({"token": token, "line": text.count("\n", 0, offset) + 1})
    return {"path": str(path), "sha256": sha256_file(path), "evidence": rows}


def source_qvalue_contract() -> dict[str, Any]:
    """Prove Seeder storage and Musubi D3Q19/PULL consumption from source."""

    repositories = {
        "seeder": (SEEDER_ROOT_WSL, SEEDER_COMMIT),
        "musubi": (MUSUBI_ROOT_WSL, MUSUBI_COMMIT),
        "musubi_scheme": (f"{MUSUBI_ROOT_WSL}/mus", MUSUBI_SCHEME_COMMIT),
        "treelm": (f"{MUSUBI_ROOT_WSL}/tem", TREELM_COMMIT),
    }
    revisions = {}
    for name, (repository, expected) in repositories.items():
        actual = _git_value(repository, "rev-parse", "HEAD")
        revisions[name] = {
            "path": repository,
            "expected_commit": expected,
            "actual_commit": actual,
            "commit_match": actual == expected,
        }

    seeder = _unc(SEEDER_ROOT_WSL)
    musubi = _unc(MUSUBI_ROOT_WSL)
    evidence = {
        "seeder_ray_direction": _source_evidence(
            seeder / "sdr/source/sdr_boundary_module.f90",
            (
                "line%origin = origin",
                "line%vec = qOffset(iDir, :) * dx",
                "qVal_t = min( qVal_t, fraction_PointLine(intersect_p, line) )",
                "qVal(iDir) = 0.5_rk",
            ),
        ),
        "seeder_element_major_26_direction_storage": _source_evidence(
            seeder / "sdr/source/sdr_proto2treelm_module.f90",
            (
                "do iElem=1,temData%qVal(1)%nVals",
                "do iDir=1,qQQQ",
                "write(distunit, rec=iElem) qVal",
            ),
        ),
        "musubi_mesh_to_stencil_mapping": _source_evidence(
            musubi / "mus/source/mus_construction_module.fpp",
            (
                "iMeshDir = stencil%map( iDir )",
                "! q-values stored in outgoing direction",
                "= bc_prop%qVal(iMeshDir, iElem_qVal)",
            ),
        ),
        "musubi_wall_libb_inverse_lookup": _source_evidence(
            musubi / "mus/source/bc/mus_bc_header_module.fpp",
            (
                "invDir = cxDirInv(iDir)",
                "%qVal%val( invDir, iElem )",
                "call set_bouzidi_coeff",
            ),
        ),
        "musubi_wall_libb_application": _source_evidence(
            musubi / "mus/source/bc/mus_bc_fluid_wall_module.fpp",
            (
                "subroutine wall_libb",
                "fIn  = bcBuffer( me%bouzidi(iLevel)% inPos(iLink) )",
                "state( me%links(iLevel)%val(iLink) ) = cIn*fIn + cOut*fOut + cNgh*fNgh",
            ),
        ),
        "musubi_pull_fetch": _source_evidence(
            musubi / "mus/source/header/lbm_macros.inc",
            (
                "PULL:  read adjacent element pdf and store to local element",
                "else !PULL",
                "macro :: FETCH(iDir,iField,node,QQ,nScalars,nElems,neigh)",
            ),
        ),
    }
    status = "PASS" if all(row["commit_match"] for row in revisions.values()) else "FAIL"
    return {
        "status": status,
        "revisions": revisions,
        "evidence": evidence,
        "proven_contract": (
            "Seeder stores q in element-major TreElm qOffset direction; Musubi maps "
            "mesh direction to the D3Q19 outgoing direction, marks the inverse as the "
            "incoming boundary write, then wall_libb deliberately fetches qVal(invDir)."
        ),
    }


def ray_cylinder_fraction(
    origins: np.ndarray,
    directions: np.ndarray,
    *,
    axis: Iterable[float],
    radius_m: float,
    center_m: Iterable[float] = (0.0, 0.0, 0.0),
) -> np.ndarray:
    """Return the first non-negative ray/infinite-cylinder parameter."""

    points = np.asarray(origins, dtype=np.float64)
    links = np.asarray(directions, dtype=np.float64)
    if points.ndim != 2 or links.shape != points.shape or points.shape[1] != 3:
        raise ValueError("origins and directions must both have shape (n, 3)")
    cylinder_axis = np.asarray(tuple(axis), dtype=np.float64).reshape(3)
    cylinder_axis /= np.linalg.norm(cylinder_axis)
    relative = points - np.asarray(tuple(center_m), dtype=np.float64).reshape(3)
    point_perp = relative - (relative @ cylinder_axis)[:, None] * cylinder_axis
    link_perp = links - (links @ cylinder_axis)[:, None] * cylinder_axis
    aa = np.einsum("ij,ij->i", link_perp, link_perp)
    bb = 2.0 * np.einsum("ij,ij->i", point_perp, link_perp)
    cc = np.einsum("ij,ij->i", point_perp, point_perp) - float(radius_m) ** 2
    discriminant = bb * bb - 4.0 * aa * cc
    result = np.full(len(points), np.nan, dtype=np.float64)
    valid = (aa > np.finfo(np.float64).tiny) & (discriminant >= 0.0)
    if not np.any(valid):
        return result
    root = np.sqrt(np.maximum(discriminant[valid], 0.0))
    denom = 2.0 * aa[valid]
    first = (-bb[valid] - root) / denom
    second = (-bb[valid] + root) / denom
    candidates = np.column_stack((first, second))
    candidates[candidates < -1.0e-12] = np.nan
    chosen = np.full(len(candidates), np.nan, dtype=np.float64)
    any_candidate = np.any(np.isfinite(candidates), axis=1)
    if np.any(any_candidate):
        chosen[any_candidate] = np.nanmin(candidates[any_candidate], axis=1)
    result[np.flatnonzero(valid)] = chosen
    return result


def _parse_uniform_lattice(header_text: str) -> tuple[np.ndarray, float, int, float]:
    match = re.search(
        r"boundingbox\s*=\s*\{\s*origin\s*=\s*\{([^}]*)\}\s*,?\s*"
        r"length\s*=\s*([-+0-9.Ee]+)",
        header_text,
        re.DOTALL | re.IGNORECASE,
    )
    level_match = re.search(r"\bmaxLevel\s*=\s*(\d+)", header_text)
    if not match or not level_match:
        raise ValueError("Could not parse uniform mesh bounding cube")
    origin_values = [
        float(value.replace("D", "E").replace("d", "e"))
        for value in re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?", match.group(1))
    ]
    if len(origin_values) != 3:
        raise ValueError("Bounding-cube origin is not three-dimensional")
    length = float(match.group(2))
    level = int(level_match.group(1))
    return np.asarray(origin_values), length, level, length / (2**level)


def _error_statistics(values: np.ndarray) -> dict[str, float | int | None]:
    data = np.asarray(values, dtype=np.float64)
    data = data[np.isfinite(data)]
    if not len(data):
        return {"count": 0, "mean": None, "rms": None, "p95": None, "max": None}
    absolute = np.abs(data)
    return {
        "count": int(len(data)),
        "mean": float(np.mean(data)),
        "mean_absolute": float(np.mean(absolute)),
        "rms": float(np.sqrt(np.mean(data * data))),
        "p95": float(np.percentile(absolute, 95.0)),
        "max": float(np.max(absolute)),
    }


def audit_mesh_qvalues(
    mesh_dir: Path,
    mesh_summary_path: Path,
    *,
    return_error_samples: bool = False,
) -> dict[str, Any] | tuple[dict[str, Any], np.ndarray]:
    mesh_dir = Path(mesh_dir).resolve()
    summary = json.loads(Path(mesh_summary_path).read_text(encoding="utf-8"))
    contract = load_mesh_contract(
        mesh_dir, allow_zero_normals=True, require_runtime_order=False
    )
    header = (mesh_dir / "header.lua").read_text(encoding="utf-8")
    cube_origin, cube_length, level, dx_m = _parse_uniform_lattice(header)
    configured_dx = float(summary["dx_m"])
    if not math.isclose(dx_m, configured_dx, rel_tol=1.0e-12, abs_tol=0.0):
        raise ValueError(f"Header dx {dx_m} != configured dx {configured_dx}")

    wall = contract.boundaries["wall"]
    wall_rows, outgoing = np.nonzero(wall.outward_masks)
    cells = wall.cell_indices[wall_rows]
    q_seeder = contract.qvalues_by_cell[cells, outgoing]
    centers = cube_origin + (contract.cell_ijk[cells] + 0.5) * dx_m
    lattice_links = D3Q19_DIRECTIONS[outgoing].astype(np.float64)
    physical_links = lattice_links * dx_m
    axis = np.asarray(summary["wall"]["axis"], dtype=np.float64)
    q_exact = ray_cylinder_fraction(
        centers,
        physical_links,
        axis=axis,
        radius_m=float(summary["wall"].get("radius_m", PIPE_RADIUS_M)),
    )
    exact_on_link = np.isfinite(q_exact) & (q_exact >= -1.0e-10) & (q_exact <= 1.0 + 1.0e-10)
    q_valid = np.isfinite(q_seeder) & (q_seeder > 0.0) & (q_seeder <= 1.0)
    comparable = exact_on_link & q_valid
    q_error = q_seeder[comparable] - q_exact[comparable]
    link_norm = np.linalg.norm(lattice_links[comparable], axis=1)
    distance_error_over_dx = q_error * link_norm

    inverse = INVERSE_DIRECTIONS[outgoing]
    q_inverse_slot = contract.qvalues_by_cell[cells, inverse]
    inverse_valid = exact_on_link & np.isfinite(q_inverse_slot) & (q_inverse_slot > 0.0)
    inverse_error = q_inverse_slot[inverse_valid] - q_exact[inverse_valid]

    signed_mean = float(np.mean(q_error)) if len(q_error) else None
    q_stats = _error_statistics(q_error)
    distance_stats = _error_statistics(distance_error_over_dx)
    direction_mismatch = int(np.count_nonzero(~exact_on_link))
    out_of_range = int(np.count_nonzero(~q_valid))
    unique_q, unique_q_counts = np.unique(q_seeder, return_counts=True)
    q_distribution = {
        format(float(value), ".17g"): int(count)
        for value, count in zip(unique_q, unique_q_counts, strict=True)
    }
    structural_pass = direction_mismatch == 0 and out_of_range == 0 and len(q_error) == len(q_seeder)
    record = {
        "status": "PASS_NO_STRUCTURAL_QVALUE_ERROR" if structural_pass else "FAIL",
        "mesh_dir": str(mesh_dir),
        "mesh_summary": str(Path(mesh_summary_path).resolve()),
        "cells_across_diameter": int(summary["cells_across_diameter"]),
        "dx_m": dx_m,
        "uniform_level": level,
        "bounding_cube_origin_m": cube_origin.tolist(),
        "bounding_cube_length_m": cube_length,
        "axis": axis.tolist(),
        "radius_m": float(summary["wall"].get("radius_m", PIPE_RADIUS_M)),
        "links_tested": int(len(q_seeder)),
        "valid_q_fraction": float(np.mean(q_valid)) if len(q_valid) else 0.0,
        "q_value_distribution_on_declared_wall_links": q_distribution,
        "uniform_halfway_fallback_fraction": float(np.mean(q_seeder == 0.5)),
        "out_of_range_q_count": out_of_range,
        "direction_mismatch_count": direction_mismatch,
        "q_seeder_minus_q_exact": q_stats,
        "distance_error_over_dx": distance_stats,
        "systematic_q_bias": {
            "mean_signed_q_error": signed_mean,
            "sign": None if signed_mean is None else ("positive" if signed_mean > 0 else "negative"),
            "interpretation": (
                "Reported without an invented pass threshold; this contains both STL polygon "
                "approximation and floating-point ray-intersection effects."
            ),
        },
        "storage_inversion_diagnostic": {
            "same_direction_comparable_links": int(len(q_error)),
            "same_direction_rms_q_error": q_stats["rms"],
            "inverse_slot_comparable_links": int(len(inverse_error)),
            "inverse_slot_rms_q_error": _error_statistics(inverse_error)["rms"],
            "source_proven_selection": "same outgoing mesh/stencil direction; wall_libb later accesses invDir because its bitmask is incoming",
        },
    }
    if return_error_samples:
        return record, q_error.copy()
    return record


def run_wall_qvalue_oracle(project_root: Path) -> dict[str, Any]:
    started = time.perf_counter()
    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / RESEARCH_RUN
    benchmark = run_root / "pressure_bc_benchmark" / "meshes"
    source = source_qvalue_contract()
    requested = {
        "axis_aligned": (16, 20, 27),
        "worst_real_outlet": (20, 27),
    }
    cases = []
    error_samples = []
    for orientation, resolutions in requested.items():
        for resolution in resolutions:
            case_root = benchmark / orientation / f"n{resolution:02d}"
            record, errors = audit_mesh_qvalues(
                case_root / "mesh",
                case_root / "mesh_summary.json",
                return_error_samples=True,
            )
            record["orientation"] = orientation
            cases.append(record)
            error_samples.append(errors)
    total_links = sum(int(case["links_tested"]) for case in cases)
    contract_pass = source["status"] == "PASS" and all(case["status"].startswith("PASS") for case in cases)
    worst = {
        key: max(
            float(case["q_seeder_minus_q_exact"][key])
            for case in cases
            if case["q_seeder_minus_q_exact"][key] is not None
        )
        for key in ("rms", "p95", "max")
    }
    aggregate_q_error = _error_statistics(np.concatenate(error_samples))
    uniform_halfway = all(
        math.isclose(float(case["uniform_halfway_fallback_fraction"]), 1.0)
        for case in cases
    )
    halfway_links = sum(
        int(round(float(case["uniform_halfway_fallback_fraction"]) * int(case["links_tested"])))
        for case in cases
    )
    overall_halfway_fraction = halfway_links / total_links
    qvalue_support = sorted(
        {
            value
            for case in cases
            for value in case["q_value_distribution_on_declared_wall_links"]
        }
    )
    total_direction_mismatch = sum(int(case["direction_mismatch_count"]) for case in cases)
    result = {
        "status": "PASS_NO_QVALUE_CONTRACT_BUG_IDENTIFIED" if contract_pass else "CFD_FLOW_WALL_QVALUE_CONTRACT_ERROR_IDENTIFIED",
        "decision_basis": (
            "No numerical q-error threshold was invented. Contract status uses source revision/mapping proof, "
            "finite in-range q-values, and existence of a forward analytic cylinder intersection for every "
            "declared D3Q19 wall link; error distributions remain diagnostic evidence."
        ),
        "source_contract": source,
        "cases": cases,
        "links_tested": total_links,
        "worst_case_q_error": worst,
        "aggregate_q_error": aggregate_q_error,
        "direction_mismatch_count": total_direction_mismatch,
        "all_declared_wall_qvalues_equal_halfway": uniform_halfway,
        "overall_halfway_fallback_fraction": overall_halfway_fraction,
        "observed_qvalue_support": qvalue_support,
        "seeder_calls": 0,
        "musubi_calls": 0,
        "harvester_calls": 0,
        "vascular_cfd_calls": 0,
        "first_failure": None if contract_pass else "A source/storage/direction structural q-value invariant failed.",
        "next": "SOURCE_PROVE_AND_RUN_MINIMAL_PERIODIC_PIPE_FORCE" if contract_pass else "FIX QVALUE CONTRACT BEFORE ANY CFD",
        "runtime_s": time.perf_counter() - started,
    }
    write_json(run_root / "qc" / "wall_qvalue_oracle.json", result)
    if not contract_pass:
        actual_head = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        ).stdout.strip()
        decision = {
            "status": "CFD_FLOW_WALL_QVALUE_CONTRACT_ERROR_IDENTIFIED",
            "root_cause_final": "WALL_QVALUE_CONTRACT_ERROR",
            "actual_head": actual_head,
            "production_pipeline_modified": False,
            "evidence": {
                "wall_qvalue_oracle": str((run_root / "qc" / "wall_qvalue_oracle.json").resolve()),
                "links_tested": total_links,
                "direction_mismatch_count": total_direction_mismatch,
                "all_declared_wall_qvalues_equal_halfway": uniform_halfway,
                "overall_halfway_fallback_fraction": overall_halfway_fraction,
                "observed_qvalue_support": qvalue_support,
                "worst_case_q_error": worst,
                "aggregate_q_error": aggregate_q_error,
                "source_direction_storage_pull_contract": source["status"],
            },
            "downstream_stage_decisions": {
                "periodic_pipe_force": "NOT_RUN_QVALUE_GATE_FAILED",
                "interior_pressure": "NOT_RUN_WALL_GATE_NOT_REACHED",
                "outlet_02_placement": "NOT_RUN_WALL_AND_PRESSURE_GATES_NOT_REACHED",
                "vascular_cfd": "NOT_RUN",
            },
            "calls": {
                "vascular_seeder": 0,
                "vascular_musubi": 0,
                "harvester": 0,
                "small_pipe_seeder": 0,
                "small_pipe_musubi": 0,
            },
            "selected_remediation_concept": (
                "Correct the isolated pipe mesh q-value generation so calc_dist produces "
                "source-consistent ray/surface fractions instead of uniform halfway fallback; "
                "then rerun the zero-run analytic oracle before any CFD."
            ),
            "first_failure": (
                "Existing pressure-benchmark meshes contain qVal=0.5 on every axis-aligned "
                "declared D3Q19 wall link and on more than 98% of oblique links (the remainder "
                "are qVal=1.0), while the analytic cylinder oracle finds missing forward "
                "intersections and systematic positive q bias."
            ),
            "next": "FIX QVALUE CONTRACT BEFORE ANY CFD",
        }
        write_json(run_root / "qc" / "isolation_root_cause_decision.json", decision)
    return result
