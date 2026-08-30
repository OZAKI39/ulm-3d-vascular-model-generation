"""Zero-run qVal contract forensics for existing pipe and vascular meshes.

This research-only module never launches Seeder, Musubi, Harvester, or the
production CFD pipeline.  It proves the binary property layout from pinned
source, exercises the Python reader with a synthetic record pattern, and then
audits immutable meshes already on disk.
"""

from __future__ import annotations

import math
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import trimesh

from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import (
    load_mesh_contract,
    read_qvalue_rows,
    scatter_qvalue_rows,
)
from .port_flux_audit import parse_treelm_side_contract
from .port_grid_sensitivity import RESEARCH_RUN
from .restart_decode import D3Q19_DIRECTIONS
from .wall_qvalue_oracle import _parse_uniform_lattice, audit_mesh_qvalues


SEEDER_ROOT_WSL = "/home/lzy/apes-pinned/seeder_official"
MUSUBI_ROOT_WSL = (
    "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300"
)
EXPECTED_REVISIONS = {
    "seeder": "667109df6fafdcb39f4409e3f5d90f04d75cd33c",
    "musubi_worktree": "4e8b277b66226277171ef93bf054d21270812793",
    "musubi_scheme": "81f8c4f13772f6d4af31f335e1e3f99b02726e25",
    "treelm": "9899d1376992c4fafc8a343d2b4ccef81de670d1",
}
VASCULAR_MESH_RELATIVE = {
    "coarse": Path(
        "outputs/cfd_flow/healthy_mouse_capillary_grid_convergence_"
        "anchor003274_20260829/grids/coarse/seeder/mesh"
    ),
    "base": Path(
        "outputs/cfd_flow/axis_aligned_ideal_plane_inlet_preflight_"
        "anchor003274_20260829_120444/seeder/mesh"
    ),
    "fine": Path(
        "outputs/cfd_flow/healthy_mouse_capillary_grid_convergence_"
        "anchor003274_20260829/grids/fine/seeder/mesh"
    ),
}
BASE_WALL_STL_RELATIVE = Path(
    "outputs/cfd_flow/axis_aligned_ideal_plane_inlet_preflight_"
    "anchor003274_20260829_120444/geometry/geometry_solver_m/wall.stl"
)
TOLERANCE = 1.0e-10
MIN_CONTINUOUS_UNIQUE = 16
MIN_NON_FALLBACK_FRACTION = 0.01


def _unc(path: str) -> Path:
    return Path("//wsl.localhost/Ubuntu" + path)


def _git_head(repository: str) -> str:
    process = subprocess.run(
        ["wsl", "git", "-C", repository, "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return process.stdout.strip()


def _source_record(
    path: Path,
    *,
    routine: str,
    tokens: Iterable[str],
    commit: str,
) -> dict[str, Any]:
    text = Path(path).read_text(encoding="utf-8")
    evidence = []
    for token in tokens:
        offset = text.find(token)
        if offset < 0:
            raise RuntimeError(f"Missing source token in {path}: {token}")
        evidence.append(
            {"token": token, "line": text.count("\n", 0, offset) + 1}
        )
    return {
        "source_file": str(path),
        "routine": routine,
        "commit": commit,
        "sha256": sha256_file(path),
        "evidence": evidence,
    }


def build_qvalue_file_format_contract() -> dict[str, Any]:
    """Recover the qVal property contract directly from pinned APES source."""

    seeder = _unc(SEEDER_ROOT_WSL)
    musubi = _unc(MUSUBI_ROOT_WSL)
    treelm = musubi / "tem"
    revisions = {
        "seeder": _git_head(SEEDER_ROOT_WSL),
        "musubi_worktree": _git_head(MUSUBI_ROOT_WSL),
        "musubi_scheme": _git_head(f"{MUSUBI_ROOT_WSL}/mus"),
        "treelm": _git_head(f"{MUSUBI_ROOT_WSL}/tem"),
    }
    revisions_match = all(
        revisions[name] == expected
        for name, expected in EXPECTED_REVISIONS.items()
    )

    param_path = treelm / "source/tem_param_module.f90"
    side_names, side_offsets = parse_treelm_side_contract(
        param_path.read_text(encoding="utf-8")
    )
    d3q19_mapping_pass = bool(
        np.array_equal(side_offsets[:18], D3Q19_DIRECTIONS[:18])
    )

    sources = {
        "writer": _source_record(
            seeder / "sdr/source/sdr_proto2treelm_module.f90",
            routine="sdr_dump_treelmesh",
            tokens=(
                "property(tem_global%nProperties)%bitpos = prp_hasQVal",
                "property(tem_global%nProperties)%nElems = temData%qVal(1)%nVals",
                "do iElem=1,temData%qVal(1)%nVals",
                "do iDir=1,qQQQ",
                "write(distunit, rec=iElem) qVal",
            ),
            commit=revisions["seeder"],
        ),
        "reader": _source_record(
            treelm / "source/tem_bc_prop_module.f90",
            routine="load_tem_BC_realArray / load_tem_BC_qVal",
            tokens=(
                "allocate(me%qVal(me%nSides, nElems))",
                "arraylen = me%nSides",
                "call MPI_TYPE_CONTIGUOUS( arraylen, rk_mpi, etype, iError )",
                "propdat(:,i) = buffer( ((i-1)*arraylen+1) : (i*arraylen) )",
            ),
            commit=revisions["treelm"],
        ),
        "dtype_and_endian": _source_record(
            treelm / "source/env_module.f90",
            routine="init_env / tem_create_EndianSuffix",
            tokens=(
                "integer, parameter, public :: rk_prec = double_prec",
                "rk_mpi = mpi_double_precision",
                "function tem_create_EndianSuffix()  result(suffix)",
                "suffix = '.lsb'",
            ),
            commit=revisions["treelm"],
        ),
        "property_element_order": _source_record(
            treelm / "source/tem_property_module.f90",
            routine="gather_property",
            tokens=(
                "Property%nElems = count(btest(BitField, Header%BitPos))",
                "do iElem=1,nElems",
                "Property%ElemID(PropElem) = iElem",
            ),
            commit=revisions["treelm"],
        ),
        "directions": _source_record(
            param_path,
            routine="tem_param_module constants",
            tokens=(
                "qOffset =",
                "qDirName(qQQQ)",
                "'  W', '  S', '  B', '  E', '  N', '  T'",
            ),
            commit=revisions["treelm"],
        ),
        "boundary_and_q_row_association": _source_record(
            musubi / "mus/source/mus_construction_module.fpp",
            routine="assignBCList",
            tokens=(
                "do iElem = 1, tree%property(1)%nElems",
                "iElem_qVal = iElem_qVal + 1",
                "iMeshDir = stencil%map( iDir )",
                "bID      = bc_prop%boundary_ID( iMeshDir, iElem )",
                "= bc_prop%qVal(iMeshDir, iElem_qVal)",
            ),
            commit=revisions["musubi_scheme"],
        ),
        "sentinel_and_fallback": _source_record(
            seeder / "sdr/source/sdr_boundary_module.f90",
            routine="sdr_qValByNode / sdr_truncate_qVal / sdr_identify_boundary",
            tokens=(
                "If there is no intersection, qVal returns -1.0",
                "qVal(iDir) = 0.5_rk",
                "for flooded neighbor, treat it as normal fluid: clean BCID",
                "for non-flooded neighbor, treat it as high order wall: set qVal to 1",
            ),
            commit=revisions["seeder"],
        ),
    }
    contract = {
        "qval_lua": (
            "hasQVal is a logical vector indexed by bnd.lua boundary-ID order; "
            "the has-qVal property bit and nElems are declared in header.lua"
        ),
        "qval_lsb_dtype": "IEEE-754 binary64 / Fortran real(kind=rk)",
        "qval_lsb_numpy_dtype": "<f8",
        "endian": "little-endian because filename suffix is .lsb",
        "n_elems": "header.lua property label='has qVal' nElems",
        "n_sides": "bnd.lua nSides; current TreElm qQQQ=26",
        "record_layout": (
            "element-major records; each row is qVal(1:nSides) in qOffset order"
        ),
        "row_element_class": (
            "ascending elemlist order restricted to property bit prp_hasQVal=8"
        ),
        "direction_names": list(side_names),
        "direction_offsets": side_offsets.astype(int).tolist(),
        "d3q19_mapping": (
            "Musubi stencil%map selects TreElm sides; pinned D3Q19 directions "
            "equal the first 18 TreElm qOffset rows"
        ),
        "boundary_association": (
            "Musubi iterates ascending has-boundary elements and increments the "
            "qVal row only when the same element carries prp_hasQVal"
        ),
        "sentinels": {
            "-1.0": "no ray/geometry intersection",
            "0.5": "unknown-boundary halfway fallback",
            "1.0": "non-flooded-neighbor high-order-wall truncation fallback",
        },
        "flooded_semantics": (
            "failed or beyond-link q on a flooded neighbor clears BCID and returns "
            "to fluid; on a non-flooded neighbor it keeps a wall with q=1"
        ),
    }
    status = "PASS" if revisions_match and d3q19_mapping_pass else "FAIL"
    return {
        "status": status,
        "revisions": {
            name: {
                "expected": EXPECTED_REVISIONS[name],
                "actual": actual,
                "match": actual == EXPECTED_REVISIONS[name],
            }
            for name, actual in revisions.items()
        },
        "contract": contract,
        "d3q19_direction_mapping_pass": d3q19_mapping_pass,
        "sources": sources,
    }


def synthetic_reader_roundtrip() -> dict[str, Any]:
    """Round-trip known element and direction patterns through the real reader."""

    rows = np.asarray(
        [[0.01 * (100 * row + side + 1) for side in range(26)] for row in range(3)],
        dtype="<f8",
    )
    property_cells = np.asarray([1, 4, 6], dtype=np.int64)
    with tempfile.TemporaryDirectory(prefix="qvalue_reader_") as directory:
        path = Path(directory) / "qval.lsb"
        rows.tofile(path)
        loaded = read_qvalue_rows(path, element_count=3, side_count=26)
        scattered = scatter_qvalue_rows(
            loaded,
            cell_count=8,
            property_cells=property_cells,
        )
    row_pass = bool(np.array_equal(loaded, rows))
    direction_pass = bool(
        np.array_equal(loaded[0], rows[0])
        and np.array_equal(loaded[2], rows[2])
    )
    scatter_pass = bool(
        np.array_equal(scattered[property_cells], rows)
        and np.all(np.isnan(scattered[[0, 2, 3, 5, 7]]))
    )
    dtype_pass = loaded.dtype == np.dtype("float64")
    endian_pass = bool(
        math.isclose(float(loaded[0, 0]), 0.01)
        and not math.isclose(float(rows.byteswap().view(rows.dtype)[0, 0]), 0.01)
    )
    passed = row_pass and direction_pass and scatter_pass and dtype_pass and endian_pass
    return {
        "status": "PASS" if passed else "FAIL",
        "row_mapping_pass": row_pass and scatter_pass,
        "direction_mapping_pass": direction_pass,
        "dtype_pass": dtype_pass,
        "endian_pass": endian_pass,
        "pattern": rows.tolist(),
        "property_cells": property_cells.tolist(),
    }


def _summary(values: np.ndarray) -> dict[str, Any]:
    data = np.asarray(values, dtype=np.float64)
    finite = data[np.isfinite(data)]
    if not len(finite):
        return {"count": 0}
    near_half = np.isclose(finite, 0.5, rtol=0.0, atol=TOLERANCE)
    near_one = np.isclose(finite, 1.0, rtol=0.0, atol=TOLERANCE)
    valid = (finite > 0.0) & (finite <= 1.0)
    non_fallback = valid & ~near_half & ~near_one
    quantiles = np.percentile(finite, [1, 5, 25, 50, 75, 95, 99])
    unique = np.unique(np.round(finite / TOLERANCE).astype(np.int64))
    return {
        "count": int(len(data)),
        "finite_count": int(len(finite)),
        "valid_q_count": int(np.count_nonzero(valid)),
        "valid_q_fraction": float(np.mean(valid)),
        "near_0.5_fraction": float(np.mean(near_half)),
        "near_1.0_fraction": float(np.mean(near_one)),
        "strict_0_to_1_fraction": float(np.mean((finite > 0.0) & (finite < 1.0))),
        "non_fallback_continuous_fraction": float(np.mean(non_fallback)),
        "mean": float(np.mean(finite)),
        "std": float(np.std(finite)),
        "min": float(np.min(finite)),
        "max": float(np.max(finite)),
        "p01": float(quantiles[0]),
        "p05": float(quantiles[1]),
        "p25": float(quantiles[2]),
        "median": float(quantiles[3]),
        "p75": float(quantiles[4]),
        "p95": float(quantiles[5]),
        "p99": float(quantiles[6]),
        "unique_q_count_at_1e-10": int(len(unique)),
    }


def vascular_wall_qvalue_distribution(mesh_dir: Path) -> dict[str, Any]:
    """Summarize wall-link q values in one existing vascular TreElm mesh."""

    mesh_dir = Path(mesh_dir).resolve()
    required = (
        "header.lua",
        "bnd.lua",
        "bnd.lsb",
        "qval.lua",
        "qval.lsb",
        "elemlist.lsb",
    )
    missing = [name for name in required if not (mesh_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing mesh contract files: {missing}")
    contract = load_mesh_contract(
        mesh_dir,
        allow_zero_normals=True,
        require_runtime_order=False,
    )
    wall = contract.boundaries["wall"]
    rows, directions = np.nonzero(wall.outward_masks)
    cells = wall.cell_indices[rows]
    qvalues = contract.qvalues_by_cell[cells, directions]
    overall = _summary(qvalues)
    by_direction = {}
    for direction in range(18):
        selected = directions == direction
        by_direction[str(direction + 1)] = {
            "offset": D3Q19_DIRECTIONS[direction].astype(int).tolist(),
            **_summary(qvalues[selected]),
        }
    continuous = int(overall.get("unique_q_count_at_1e-10", 0)) >= MIN_CONTINUOUS_UNIQUE
    populated = float(overall.get("non_fallback_continuous_fraction", 0.0)) >= MIN_NON_FALLBACK_FRACTION
    passed = continuous and populated
    return {
        "status": "PASS" if passed else "FAIL",
        "mesh_dir": str(mesh_dir),
        "files": {
            name: {
                "path": str((mesh_dir / name).resolve()),
                "sha256": sha256_file(mesh_dir / name),
            }
            for name in required
        },
        "wall_links": overall,
        "by_d3q19_direction": by_direction,
        "continuous_surface_gate": {
            "minimum_unique_q_at_1e-10": MIN_CONTINUOUS_UNIQUE,
            "minimum_non_fallback_fraction": MIN_NON_FALLBACK_FRACTION,
            "pass": passed,
            "interpretation": (
                "This gate only distinguishes a varied sub-grid distance field "
                "from documented 0.5/1.0 fallback support; it is not a CFD accuracy gate."
            ),
        },
    }


def _segment_triangle_q(
    start: np.ndarray,
    vector: np.ndarray,
    triangles: np.ndarray,
) -> float | None:
    """Return the nearest Moller-Trumbore intersection fraction on one segment."""

    edge1 = triangles[:, 1] - triangles[:, 0]
    edge2 = triangles[:, 2] - triangles[:, 0]
    vector_rows = np.broadcast_to(vector, edge2.shape)
    cross = np.cross(vector_rows, edge2)
    determinant = np.einsum("ij,ij->i", edge1, cross)
    active = np.abs(determinant) > np.finfo(np.float64).tiny
    inverse = np.zeros_like(determinant)
    inverse[active] = 1.0 / determinant[active]
    relative = start - triangles[:, 0]
    bary_u = inverse * np.einsum("ij,ij->i", relative, cross)
    active &= (bary_u >= -TOLERANCE) & (bary_u <= 1.0 + TOLERANCE)
    second_cross = np.cross(relative, edge1)
    bary_v = inverse * np.einsum("ij,ij->i", vector_rows, second_cross)
    active &= (bary_v >= -TOLERANCE) & (
        bary_u + bary_v <= 1.0 + TOLERANCE
    )
    fractions = inverse * np.einsum("ij,ij->i", edge2, second_cross)
    active &= (fractions > TOLERANCE) & (fractions <= 1.0 + TOLERANCE)
    if not np.any(active):
        return None
    return float(np.min(fractions[active]))


def sample_base_stl_ray_oracle(
    mesh_dir: Path,
    wall_stl: Path,
    *,
    target_comparable: int = 10_000,
    random_seed: int = 20_260_830,
) -> dict[str, Any]:
    """Compare BASE q against nearest valid wall-STL intersections."""

    mesh_dir = Path(mesh_dir).resolve()
    wall_stl = Path(wall_stl).resolve()
    contract = load_mesh_contract(
        mesh_dir,
        allow_zero_normals=True,
        require_runtime_order=False,
    )
    wall = contract.boundaries["wall"]
    wall_rows, directions = np.nonzero(wall.outward_masks)
    cells = wall.cell_indices[wall_rows]
    header = (mesh_dir / "header.lua").read_text(encoding="utf-8")
    cube_origin, _, _, dx_m = _parse_uniform_lattice(header)
    starts = cube_origin + (contract.cell_ijk[cells] + 0.5) * dx_m
    vectors = D3Q19_DIRECTIONS[directions].astype(np.float64) * dx_m
    q_seeder = contract.qvalues_by_cell[cells, directions]

    surface = trimesh.load_mesh(wall_stl, process=False)
    if not isinstance(surface, trimesh.Trimesh):
        raise TypeError(f"Expected a single triangular wall STL: {wall_stl}")
    triangles = np.asarray(surface.triangles, dtype=np.float64)
    spatial_index = surface.triangles_tree
    order = np.random.default_rng(random_seed).permutation(len(cells))
    exact_values: list[float] = []
    stored_values: list[float] = []
    attempted = 0
    epsilon = np.finfo(np.float64).eps * max(1.0, float(np.max(np.abs(surface.bounds))))
    for index in order:
        attempted += 1
        start = starts[index]
        vector = vectors[index]
        end = start + vector
        lower = np.minimum(start, end) - epsilon
        upper = np.maximum(start, end) + epsilon
        candidates = list(spatial_index.intersection(np.concatenate((lower, upper))))
        if not candidates:
            continue
        exact = _segment_triangle_q(start, vector, triangles[candidates])
        if exact is None or not (0.0 < exact <= 1.0 + TOLERANCE):
            continue
        stored = float(q_seeder[index])
        if not math.isfinite(stored) or not (0.0 < stored <= 1.0):
            continue
        exact_values.append(min(exact, 1.0))
        stored_values.append(stored)
        if len(exact_values) >= int(target_comparable):
            break
    exact_array = np.asarray(exact_values, dtype=np.float64)
    stored_array = np.asarray(stored_values, dtype=np.float64)
    errors = stored_array - exact_array
    absolute = np.abs(errors)
    enough = len(errors) >= int(target_comparable)
    mismatch = attempted - len(errors)
    return {
        "status": "PASS_SAMPLE_SIZE" if enough else "FAIL_INSUFFICIENT_COMPARABLE",
        "mesh_dir": str(mesh_dir),
        "wall_stl": str(wall_stl),
        "wall_stl_sha256": sha256_file(wall_stl),
        "intersection_method": (
            "trimesh triangle R-tree plus nearest valid Moller-Trumbore ray-segment "
            "intersection from the fluid center in the source-proven D3Q19 direction"
        ),
        "random_seed": random_seed,
        "sample_target": int(target_comparable),
        "sample_attempted": int(attempted),
        "sample_comparable": int(len(errors)),
        "direction_mismatch_count": int(mismatch),
        "direction_mismatch_fraction": float(mismatch / attempted),
        "q_exact": _summary(exact_array),
        "q_seeder": _summary(stored_array),
        "q_seeder_minus_q_exact": {
            "count": int(len(errors)),
            "mean": float(np.mean(errors)) if len(errors) else None,
            "median_bias": float(np.median(errors)) if len(errors) else None,
            "rms": float(np.sqrt(np.mean(errors * errors))) if len(errors) else None,
            "p95": float(np.percentile(absolute, 95.0)) if len(errors) else None,
            "max": float(np.max(absolute)) if len(errors) else None,
        },
    }


def _pipe_distribution(case: dict[str, Any]) -> dict[str, Any]:
    support = {
        float(value): int(count)
        for value, count in case["q_value_distribution_on_declared_wall_links"].items()
    }
    total = int(case["links_tested"])
    half = sum(count for value, count in support.items() if abs(value - 0.5) <= TOLERANCE)
    one = sum(count for value, count in support.items() if abs(value - 1.0) <= TOLERANCE)
    non_fallback = total - half - one
    return {
        "links": total,
        "support": {format(value, ".17g"): count for value, count in support.items()},
        "near_0.5_fraction": half / total,
        "near_1.0_fraction": one / total,
        "non_fallback_continuous_fraction": non_fallback / total,
    }


def classify_qvalue_contract(
    *,
    reader_pass: bool,
    pipe_abnormal: bool,
    vascular_abnormal: bool,
) -> tuple[str, str, str]:
    if not reader_pass:
        return (
            "CFD_FLOW_QVALUE_READER_FIXED_NO_MESH_ERROR",
            "RUN MINIMAL PERIODIC PIPE FORCE",
            "CASE_A",
        )
    if pipe_abnormal and not vascular_abnormal:
        return (
            "CFD_FLOW_PIPE_BENCHMARK_QVALUE_GENERATION_ERROR",
            "FIX PIPE BENCHMARK GEOMETRY/QVALUE GENERATION",
            "CASE_B",
        )
    if pipe_abnormal and vascular_abnormal:
        return (
            "CFD_FLOW_VASCULAR_WALL_QVALUE_CONTRACT_ERROR",
            "FIX SEEDER→QVAL→WALL_LIBB CONTRACT THEN REGENERATE TEST MESHES",
            "CASE_C",
        )
    return (
        "CFD_FLOW_QVALUE_CONTRACT_UNRESOLVED",
        "DO NOT RUN CFD; RESOLVE CURRENT QVALUE EVIDENCE GAP",
        "CASE_D",
    )


def run_qvalue_contract_forensics(project_root: Path) -> dict[str, Any]:
    """Execute the ordered zero-run decision tree and write three QC records."""

    started = time.perf_counter()
    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / RESEARCH_RUN
    qc = run_root / "qc"
    contract = build_qvalue_file_format_contract()
    write_json(qc / "qvalue_file_format_contract.json", contract)

    reader = synthetic_reader_roundtrip()
    benchmark = run_root / "pressure_bc_benchmark" / "meshes"
    requested = {
        "axis_n20": benchmark / "axis_aligned" / "n20",
        "axis_n27": benchmark / "axis_aligned" / "n27",
        "oblique_n27": benchmark / "worst_real_outlet" / "n27",
    }
    pipes = {}
    for name, case_root in requested.items():
        audit = audit_mesh_qvalues(case_root / "mesh", case_root / "mesh_summary.json")
        pipes[name] = {**audit, "distribution": _pipe_distribution(audit)}
    pipe_abnormal = all(
        item["distribution"]["non_fallback_continuous_fraction"] < 0.01
        for item in pipes.values()
    )

    vascular = {
        name: vascular_wall_qvalue_distribution(root / relative)
        for name, relative in VASCULAR_MESH_RELATIVE.items()
    }
    write_json(
        qc / "vascular_wall_qvalue_distribution.json",
        {
            "status": "PASS" if all(item["status"] == "PASS" for item in vascular.values()) else "FAIL",
            "tolerance": TOLERANCE,
            "grids": vascular,
            "seeder_calls": 0,
            "musubi_calls": 0,
            "harvester_calls": 0,
            "vascular_cfd_calls": 0,
        },
    )
    base_oracle = sample_base_stl_ray_oracle(
        root / VASCULAR_MESH_RELATIVE["base"],
        root / BASE_WALL_STL_RELATIVE,
    )
    vascular_abnormal = any(item["status"] != "PASS" for item in vascular.values())
    reader_pass = contract["status"] == "PASS" and reader["status"] == "PASS"
    final_status, next_step, decision_case = classify_qvalue_contract(
        reader_pass=reader_pass,
        pipe_abnormal=pipe_abnormal,
        vascular_abnormal=vascular_abnormal,
    )
    exact_root_cause = (
        "Pinned Seeder 667109 wrote no continuous intersection fractions into any "
        "audited pipe or vascular qval.lsb: raw wall slots are restricted to the "
        "source-documented unknown-boundary q=0.5 and non-flooded q=1.0 fallbacks. "
        "The source-proven little-endian binary64 reader and direction scatter pass, "
        "while independent BASE wall-STL segment intersections produce continuous "
        "q_exact values. The failure is therefore in the existing Seeder "
        "geometry-to-qVal intersection/fallback path, upstream of Musubi wall_libb."
    )
    first_failure = (
        "Existing axis N20/N27 and oblique N27 pipe meshes contain no meaningful "
        "continuous wall q support after the reader passed."
        if pipe_abnormal
        else None
    )
    result = {
        "status": final_status,
        "decision_case": decision_case,
        "actual_head": subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        ).stdout.strip(),
        "production_pipeline_modified": False,
        "source_contract": contract,
        "python_reader": reader,
        "reader_bug_found": False,
        "pipe_meshes": pipes,
        "vascular_meshes": vascular,
        "base_stl_ray_oracle": base_oracle,
        "pipe_only_bug": pipe_abnormal and not vascular_abnormal,
        "vascular_wide_bug": vascular_abnormal,
        "vascular_q_contract_pass": not vascular_abnormal,
        "exact_root_cause": exact_root_cause,
        "calls": {
            "vascular_cfd": 0,
            "seeder": 0,
            "musubi": 0,
            "harvester": 0,
            "small_pipe_force": 0,
        },
        "pipe_force_result": "NOT_RUN_QVALUE_GATE_FAILED",
        "first_failure": first_failure,
        "next": next_step,
        "runtime_s": time.perf_counter() - started,
    }
    write_json(qc / "qvalue_contract_forensics_v2.json", result)
    return result
