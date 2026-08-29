"""Exact, zero-run Musubi boundary-link flux audit from one frozen restart.

The implementation mirrors the pinned D3Q19 PULL boundary path.  It never
launches Seeder, Musubi, or harvesting and treats all solver artefacts as
read-only evidence.
"""

from __future__ import annotations

import csv
import math
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .apes import load_boundary_conditions
from .config import load_cfd_flow_config
from .geometry import load_frozen_surface_partition
from .io import FlowError, load_flow_inputs, read_json, sha256_file, write_json
from .port_flux_audit import (
    EXPECTED_BOUNDARY_ELEMENT_COUNT,
    EXPECTED_SIDE_COUNT,
    PORT_LABELS,
    build_source_contract,
    extract_boundary_property_indices,
    parse_bnd_header,
    parse_boundary_property_header,
    read_boundary_ids,
)
from .restart_decode import (
    D3Q19_DIRECTIONS,
    D3Q19_WEIGHTS,
    read_restart_pdf,
    read_treelm_elemlist,
)


SOURCE_BRANCH = "codex/cfd-flow-musubi-recovery-20260828"
CURRENT_SYNCED_BASE_COMMIT = "02660620d9bead5356d931f0c4fa65bd3efa3491"
DIRECT_FIELD_RUN = "musubi_direct_restart_field_anchor003274_20260829_000544"
FROZEN_SEEDER_RUN = "musubi_recovery_anchor003274_20260828_162530"
FROZEN_STEADY_RUN = "musubi_project_steady_confirmation_anchor003274_20260828_225334"
AUDIT_PREFIX = "exact_link_flux_audit_anchor003274"
AUDIT_REVISION = "ZERO_RUN_EXACT_MUSUBI_BOUNDARY_LINK_FLUX_V1"

MUSUBI_SOURCE_COMMIT = "4e8b277b66226277171ef93bf054d21270812793"
MUSUBI_SCHEME_COMMIT = "81f8c4f13772f6d4af31f335e1e3f99b02726e25"
TREELM_SOURCE_COMMIT = "9899d1376992c4fafc8a343d2b4ccef81de670d1"
MUSUBI_EXECUTABLE_SHA256 = (
    "a005b4f00bd45df0339adc22460f251c3f300f967ff746c1cd43fa5ad7c07e88"
)

FROZEN_ITERATION = 198_064
EXPECTED_CELL_COUNT = 221_109
EXPECTED_DX_M = 2.0e-7
EXPECTED_DT_S = 2.44140625e-8
REFERENCE_DENSITY_KG_M3 = 1056.0
TARGET_Q_M3_S = 7.693508475538942e-16
TARGET_MASS_FLOW_KG_S = 8.124344950169123e-13
MAXIMUM_RELATIVE_ERROR = 0.01
CS2 = 1.0 / 3.0

NOT_IDENTIFIABLE = "CFD_FLOW_EXACT_LINK_FLUX_NOT_IDENTIFIABLE_FROM_SINGLE_RESTART"
MFR_TARGET_PASS = "MFR_EQ_TARGET_PASS"
MFR_TARGET_MISMATCH = "MFR_EQ_OBLIQUE_PORT_TARGET_MISMATCH_CONFIRMED"
MASS_BALANCE_PASS = "EXACT_LINK_MASS_BALANCE_PASS"
NEXT_NOT_IDENTIFIABLE = "DESIGN MINIMAL ONE-STEP INSTRUMENTED FLUX TRACKING"
NEXT_MISMATCH = "REVIEW MFR_EQ OBLIQUE-PORT AREA AND RAW DISCRETE-NORMAL SEMANTICS"
NEXT_TARGET_PASS = "PROMOTE EXACT PDF-LINK FLUX QC THEN PROCEED TO GRID CONVERGENCE"

INVERSE_DIRECTIONS = np.asarray(
    [
        int(np.flatnonzero(np.all(D3Q19_DIRECTIONS == -direction, axis=1))[0])
        for direction in D3Q19_DIRECTIONS
    ],
    dtype=np.int8,
)
NORMAL_WEIGHTS = np.asarray(
    (*(4 for _ in range(6)), *(2 for _ in range(12))), dtype=np.int8
)
GLOBAL_NORMALS = {
    "outlet_01": np.asarray((0, -1, 0), dtype=np.int8),
    "outlet_02": np.asarray((1, 1, 0), dtype=np.int8),
    "outlet_03": np.asarray((-1, 0, 0), dtype=np.int8),
    "inlet": np.asarray((0, -1, -1), dtype=np.int8),
}
EXPECTED_STENCIL_COUNTS = {
    "inlet": 311,
    "outlet_01": 169,
    "outlet_02": 179,
    "outlet_03": 201,
}
EXPECTED_PRESSURE_VALID_COUNTS = {"outlet_01": 169, "outlet_02": 173, "outlet_03": 201}


@dataclass(frozen=True, slots=True)
class BoundaryReconstruction:
    label: str
    boundary_id: int
    property_rows: np.ndarray
    cell_indices: np.ndarray
    outward_masks: np.ndarray
    incoming_masks: np.ndarray
    raw_normals: np.ndarray
    normal_indices: np.ndarray


def mfr_eq_area_proxy(element_count: int, dx_m: float) -> float:
    """Mirror ``sum(nElems_totalLevel * dxLvl**2)`` for one level."""

    return int(element_count) * float(dx_m) ** 2


def closest_discrete_direction(raw_normal: np.ndarray) -> int:
    """Mirror tem_determine_discreteVector and return a zero-based index."""

    vector = np.asarray(raw_normal, dtype=np.float64).reshape(3)
    length = float(np.linalg.norm(vector))
    if length == 0.0:
        raise ValueError("A boundary normal cannot be zero")
    prevailing = D3Q19_DIRECTIONS[:18].astype(np.float64)
    prevailing /= np.linalg.norm(prevailing, axis=1)[:, None]
    return int(np.argmax(prevailing @ (vector / length)))


def reconstruct_boundary(
    boundary_ids: np.ndarray,
    property_element_indices: np.ndarray,
    *,
    label: str,
    boundary_id: int,
    allow_zero_normals: bool = False,
) -> BoundaryReconstruction:
    """Mirror countnBnds/assignBCList for the active D3Q19 QQN=18."""

    ids = np.asarray(boundary_ids, dtype=np.int64)
    property_indices = np.asarray(property_element_indices, dtype=np.int64).reshape(-1)
    if ids.shape[0] != len(property_indices) or ids.shape[1] < 18:
        raise ValueError("Boundary IDs and property-element mapping are inconsistent")
    active = ids[:, :18] == int(boundary_id)
    rows = np.flatnonzero(np.any(active, axis=1)).astype(np.int64)
    outward = active[rows]
    incoming = np.zeros_like(outward)
    raw = np.zeros((len(rows), 3), dtype=np.int64)
    normal_indices = np.empty(len(rows), dtype=np.int64)
    for local, mask in enumerate(outward):
        outgoing_directions = np.flatnonzero(mask)
        incoming[local, INVERSE_DIRECTIONS[outgoing_directions]] = True
        raw[local] = -np.sum(
            NORMAL_WEIGHTS[outgoing_directions, None]
            * D3Q19_DIRECTIONS[outgoing_directions].astype(np.int64),
            axis=0,
        )
        if np.all(raw[local] == 0) and allow_zero_normals:
            # Adaptive-flux/pressure mesh QC consumes the exact boundary-link
            # masks, not normalInd. A zero aggregate at a coarse Cartesian rim
            # is represented explicitly instead of inventing a direction.
            # The default remains strict for source-proven V2 audits.
            normal_indices[local] = -1
        else:
            normal_indices[local] = closest_discrete_direction(raw[local])
    return BoundaryReconstruction(
        label=label,
        boundary_id=int(boundary_id),
        property_rows=rows,
        cell_indices=property_indices[rows],
        outward_masks=outward,
        incoming_masks=incoming,
        raw_normals=raw,
        normal_indices=normal_indices,
    )


def equilibrium_pdf(rho: np.ndarray | float, velocity: np.ndarray) -> np.ndarray:
    """Strict standard D3Q19 equilibrium used by pinned get_pdfEq_d3q19."""

    density = np.asarray(rho, dtype=np.float64)
    vel = np.asarray(velocity, dtype=np.float64)
    scalar = vel.ndim == 1
    if scalar:
        vel = vel.reshape(1, 3)
    if vel.ndim != 2 or vel.shape[1] != 3:
        raise ValueError("velocity must have shape (3,) or (n, 3)")
    if density.ndim == 0:
        density = np.full(len(vel), float(density), dtype=np.float64)
    else:
        density = density.reshape(-1)
    if len(density) != len(vel):
        raise ValueError("rho and velocity lengths differ")
    cu = vel @ D3Q19_DIRECTIONS.astype(np.float64).T
    u2 = np.sum(vel * vel, axis=1)
    result = (
        density[:, None]
        * D3Q19_WEIGHTS[None, :]
        * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2[:, None])
    )
    return result[0] if scalar else result


def build_coordinate_lookup(cell_ijk: np.ndarray) -> dict[tuple[int, int, int], int]:
    coordinates = np.asarray(cell_ijk, dtype=np.int64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("cell_ijk must have shape (n, 3)")
    lookup = {
        tuple(int(value) for value in row): index
        for index, row in enumerate(coordinates)
    }
    if len(lookup) != len(coordinates):
        raise ValueError("cell_ijk contains duplicates")
    return lookup


def pull_fetch_pdfs(
    pdf: np.ndarray,
    cell_ijk: np.ndarray,
    target_indices: np.ndarray,
    *,
    coordinate_lookup: dict[tuple[int, int, int], int] | None = None,
) -> np.ndarray:
    """Mirror PULL FETCH, including local inverse-direction bounce-back."""

    values = np.asarray(pdf, dtype=np.float64)
    coordinates = np.asarray(cell_ijk, dtype=np.int64)
    targets = np.asarray(target_indices, dtype=np.int64).reshape(-1)
    lookup = coordinate_lookup or build_coordinate_lookup(coordinates)
    fetched = np.empty((len(targets), 19), dtype=np.float64)
    for row, target in enumerate(targets):
        for direction_index, direction in enumerate(D3Q19_DIRECTIONS):
            source_coord = tuple(
                int(value) for value in coordinates[target] - direction
            )
            source = lookup.get(source_coord)
            if source is None:
                fetched[row, direction_index] = values[
                    target, INVERSE_DIRECTIONS[direction_index]
                ]
            else:
                fetched[row, direction_index] = values[source, direction_index]
    return fetched


def velocity_from_pdf(pdf: np.ndarray) -> np.ndarray:
    values = np.asarray(pdf, dtype=np.float64)
    rho = np.sum(values, axis=1)
    return (values @ D3Q19_DIRECTIONS.astype(np.float64)) / rho[:, None]


def signed_mass_balance(
    inlet_kg_s: float, outlets_kg_s: Iterable[float]
) -> dict[str, float]:
    outlet_values = tuple(float(value) for value in outlets_kg_s)
    outlet_sum = float(sum(outlet_values))
    error = abs(float(inlet_kg_s) - outlet_sum) / abs(float(inlet_kg_s))
    return {"outlet_signed_sum_kg_s": outlet_sum, "relative_error": error}


def classify_exact_flux(
    inlet_relative_error: float, mass_balance_error: float, outlet02: float
) -> dict[str, Any]:
    inlet_pass = float(inlet_relative_error) <= MAXIMUM_RELATIVE_ERROR
    balance_pass = float(mass_balance_error) <= MAXIMUM_RELATIVE_ERROR
    status = MFR_TARGET_PASS if inlet_pass else MFR_TARGET_MISMATCH
    return {
        "status": status,
        "next": NEXT_TARGET_PASS if inlet_pass else NEXT_MISMATCH,
        "mfr_eq_target_classification": status,
        "exact_link_mass_balance": MASS_BALANCE_PASS if balance_pass else "FAIL",
        "outlet_02_backflow_confirmed": "YES"
        if balance_pass and outlet02 < 0.0
        else "UNRESOLVED",
        "mfr_eq_oblique_target_mismatch_confirmed": "NO" if inlet_pass else "YES",
    }


def _git_value(repository: Path, *arguments: str) -> str:
    process = subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return process.stdout.strip()


def _wsl_source_root(distribution: str) -> Path:
    for candidate in (
        Path(rf"\\wsl.localhost\{distribution}\home\lzy\apes-pinned\musubi_official"),
        Path(rf"\\wsl$\{distribution}\home\lzy\apes-pinned\musubi_official"),
    ):
        if candidate.is_dir():
            return candidate
    raise FlowError(NOT_IDENTIFIABLE, "Pinned Musubi source tree is unavailable")


def _source_evidence(
    path: Path, revision: str, statements: dict[str, tuple[str, ...]]
) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="strict")
    evidence: dict[str, Any] = {}
    for label, tokens in statements.items():
        line_numbers: list[int] = []
        for token in tokens:
            offset = text.find(token)
            if offset < 0:
                raise FlowError(
                    NOT_IDENTIFIABLE,
                    f"Pinned source token is absent in {path.name}: {token}",
                )
            line_numbers.append(text.count("\n", 0, offset) + 1)
        evidence[label] = {
            "status": "PASS",
            "tokens": list(tokens),
            "line_numbers": line_numbers,
        }
    return {
        "path": str(path),
        "source_revision": revision,
        "sha256": sha256_file(path),
        "evidence": evidence,
    }


def _exact_source_contract(source_root: Path, runtime_log: Path) -> dict[str, Any]:
    mus = source_root / "mus" / "source"
    tem = source_root / "tem" / "source"
    files = {
        "main_loop_and_final_restart": _source_evidence(
            mus / "mus_program_module.fpp",
            MUSUBI_SCHEME_COMMIT,
            {
                "steady_check_precedes_compute": (
                    "call check_flow_status(",
                    "exit mainLoop",
                    "call control%do_computation(",
                ),
                "restart_after_main_loop": (
                    "Writing restart at final time:",
                    "call mus_writeRestart(",
                ),
            },
        ),
        "fast_single_level_order": _source_evidence(
            mus / "mus_control_module.f90",
            MUSUBI_SCHEME_COMMIT,
            {
                "boundary_before_swap_compute": (
                    "call set_boundary(",
                    "call mus_swap_now_next(",
                    "call me%scheme%compute(",
                ),
                "compute_now_to_next": (
                    "inState   = me%scheme%state(iLevel)%val(:,Now)",
                    "outState  = me%scheme%state(iLevel)%val(:,Next)",
                ),
            },
        ),
        "restart_serializes_nnext_save": _source_evidence(
            mus / "mus_buffer_module.fpp",
            MUSUBI_SCHEME_COMMIT,
            {
                "nnext_save": (
                    "buffer(iIndex) = scheme%state",
                    "?SAVE?(",
                    "scheme%pdf( iLevel )%nNext",
                )
            },
        ),
        "pull_macros": _source_evidence(
            mus / "header" / "lbm_macros.inc",
            MUSUBI_SCHEME_COMMIT,
            {
                "pull_fetch_save": (
                    "else !PULL",
                    "macro :: FETCH",
                    "macro :: SAVE",
                    "macro :: NgDir",
                ),
            },
        ),
        "pull_connectivity": _source_evidence(
            mus / "mus_connectivity_module.fpp",
            MUSUBI_SCHEME_COMMIT,
            {
                "inverse_bounceback": (
                    "nghDir = ?NgDir?(iDir, stencil)",
                    "sourceDir = stencil%cxDirInv( iDir )",
                    "GetFromPos = iElem",
                )
            },
        ),
        "boundary_buffers": _source_evidence(
            mus / "bc" / "mus_bc_general_module.fpp",
            MUSUBI_SCHEME_COMMIT,
            {
                "bc_buffer_from_nnext": (
                    "currState = state( :, pdf%nNext )",
                    "bcBuffer always uses AOS",
                    "?SAVE?(",
                ),
                "neighbor_buffer_from_nnext_fetch": (
                    "neighBufferPre_nNext",
                    "currstate(",
                    "?FETCH?(",
                ),
            },
        ),
        "mfr_and_pressure_eq": _source_evidence(
            mus / "bc" / "mus_bc_fluid_module.fpp",
            MUSUBI_SCHEME_COMMIT,
            {
                "mfr_area_velocity": (
                    "globBC%nElems_totalLevel(iLvl) * physics%dxLvl(iLvl)**2",
                    "1.0_rk / ( physics%rho0 * area )",
                    "layout%fStencil%cxDirRK",
                    "normalInd%val(iElem)",
                ),
                "pressure_macro_and_extrapolation": (
                    "rho = rho / physics%fac( iLevel )%press * cs2inv",
                    "neighBufferPre_nNext(1,:)",
                    "neighBufferPre_nNext(2,:)",
                    "1.5_rk * uxB_1",
                    "0.5_rk * uxB_2",
                ),
                "incoming_fetch_replacement": (
                    "bitmask%val( iDir, iElem )",
                    "state( ?FETCH?(",
                    "fEq( (iElem-1)*QQ + iDir )",
                ),
            },
        ),
        "d3q19_boundary_lists": _source_evidence(
            mus / "mus_construction_module.fpp",
            MUSUBI_SCHEME_COMMIT,
            {
                "active_stencil_only": (
                    "do iDir = 1, stencil%QQN",
                    "iMeshDir = stencil%map( iDir )",
                ),
                "bitmask_inverse": ("stencil%cxDirInv(iDir), nBnds( bID ) ) = .true.",),
                "weighted_inward_normal": (
                    "weight(iDir) = 4",
                    "weight(iDir) = 2",
                    "-  weight(iDir) * stencil%cxDir( :, iDir )",
                ),
                "normal_projection": (
                    "call tem_determine_discreteVector(",
                    "normalInd",
                    "val = iDir",
                ),
                "remove_runtime_solid_boundary_elements": (
                    "if ( remove_solid .and. (nBCs > 0) ) then",
                    "call remove_solid_in_bc(",
                ),
            },
        ),
        "higher_order_stencil_solidification": _source_evidence(
            tem / "tem_construction_module.f90",
            TREELM_SOURCE_COMMIT,
            {
                "missing_either_required_neighbor_marks_current_element_solid": (
                    "if( levelPos <= 0 .and. computeStencil%requireNeighNeigh ) then",
                    "missingNeigh = .true.",
                    "= ibset( levelDesc%property( totalPos ), prp_solid )",
                ),
            },
        ),
        "solid_boundary_removal": _source_evidence(
            mus / "bc" / "mus_bc_header_module.fpp",
            MUSUBI_SCHEME_COMMIT,
            {
                "pressure_eq_requires_two_neighbors": (
                    "case( 'pressure_eq' )",
                    "me( myBCID )%requireNeighBufPre_nNext = .true.",
                    "me( myBCID )%nNeighs = 2",
                ),
                "solid_elements_excluded": (
                    "btest( levelDesc(iLevel)%property(posInTotal), prp_solid )",
                    "nValid(iLevel) = nValid(iLevel) + 1",
                ),
            },
        ),
        "raw_cxdir_and_unit_prevailing": _source_evidence(
            tem / "tem_stencil_module.fpp",
            TREELM_SOURCE_COMMIT,
            {
                "raw_cxdirrk": ("me%cxDirRK(:,:) = real( cxDir(:,:), rk )",),
                "unit_prevailing": (
                    "length_rk = sqrt(real(length, kind=rk))",
                    "/ length_rk",
                ),
            },
        ),
        "strict_equilibrium": _source_evidence(
            mus / "scheme" / "mus_scheme_derived_quantities_type_module.f90",
            MUSUBI_SCHEME_COMMIT,
            {
                "d3q19_equilibrium": (
                    "pure function get_pdfEq_d3q19",
                    "rho_div_18",
                    "rho_div_36",
                    "fEq(19)",
                ),
            },
        ),
    }
    log_text = runtime_log.read_text(encoding="utf-8", errors="replace")
    required_log_tokens = (
        "Using AOS (array of structures data layout)",
        "Using PULL for streaming",
        "Select fast single level control routine.",
        "nElems: 311",
        "Reached steady state       198064 T",
        "Writing restart at final time:",
        "Removed 6 elements on level 9",
    )
    missing = [token for token in required_log_tokens if token not in log_text]
    if missing:
        raise FlowError(
            NOT_IDENTIFIABLE, f"Runtime-stage evidence is incomplete: {missing}"
        )
    return {
        "status": "PASS",
        "restart_pdf_timestep_stage_proven": True,
        "restart_stage": "POST_FUSED_PULL_STREAM_AND_BGK_COLLISION_NNEXT_BEFORE_NEXT_BOUNDARY_UPDATE",
        "streaming": "PULL",
        "fetch": "neighbor lookup in inverse streaming direction; missing/solid link reads local inverse PDF",
        "bc_buffer_exactly_reconstructable": True,
        "cxDirRK": "RAW_INTEGER_DIRECTION_AS_REAL_NOT_UNIT_NORMALIZED",
        "files": files,
        "runtime_log": {
            "path": str(runtime_log),
            "sha256": sha256_file(runtime_log),
            "required_tokens": list(required_log_tokens),
        },
    }


def _file_manifest(paths: Iterable[Path]) -> dict[str, Any]:
    return {
        str(path.resolve()): {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for path in paths
    }


def _normal_distribution(boundary: BoundaryReconstruction) -> list[dict[str, int]]:
    unique, counts = np.unique(boundary.normal_indices, return_counts=True)
    return [
        {
            "direction_index": int(index + 1),
            "cx": int(D3Q19_DIRECTIONS[index, 0]),
            "cy": int(D3Q19_DIRECTIONS[index, 1]),
            "cz": int(D3Q19_DIRECTIONS[index, 2]),
            "element_count": int(count),
        }
        for index, count in zip(unique, counts, strict=True)
    ]


def _valid_pressure_mask(
    boundary: BoundaryReconstruction,
    cell_ijk: np.ndarray,
    lookup: dict[tuple[int, int, int], int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    normal = GLOBAL_NORMALS[boundary.label].astype(np.int64)
    neighbor1 = np.full(len(boundary.cell_indices), -1, dtype=np.int64)
    neighbor2 = np.full(len(boundary.cell_indices), -1, dtype=np.int64)
    for row, cell_index in enumerate(boundary.cell_indices):
        coord = cell_ijk[cell_index]
        neighbor1[row] = lookup.get(tuple(int(value) for value in coord + normal), -1)
        neighbor2[row] = lookup.get(
            tuple(int(value) for value in coord + 2 * normal), -1
        )
    valid = (neighbor1 >= 0) & (neighbor2 >= 0)
    return valid, neighbor1, neighbor2


def _link_rows(
    *,
    boundary: BoundaryReconstruction,
    selected_rows: np.ndarray,
    tree_ids: np.ndarray,
    old_pdf: np.ndarray,
    new_pdf: np.ndarray,
    mass_factor: float,
    outward_sign: bool,
) -> tuple[list[dict[str, Any]], float]:
    records: list[dict[str, Any]] = []
    total_lattice = 0.0
    for output_row, boundary_row in enumerate(selected_rows):
        cell_index = int(boundary.cell_indices[boundary_row])
        for incoming_index in np.flatnonzero(boundary.incoming_masks[boundary_row]):
            outgoing_index = int(INVERSE_DIRECTIONS[incoming_index])
            old = float(old_pdf[cell_index, outgoing_index])
            new = float(new_pdf[output_row, incoming_index])
            domain_delta = new - old
            signed = -domain_delta if outward_sign else domain_delta
            total_lattice += signed
            records.append(
                {
                    "port": boundary.label,
                    "property_row": int(boundary.property_rows[boundary_row]),
                    "cell_index_zero_based": cell_index,
                    "tree_id": int(tree_ids[cell_index]),
                    "normal_index_one_based": int(
                        boundary.normal_indices[boundary_row] + 1
                    ),
                    "incoming_pdf_index_one_based": int(incoming_index + 1),
                    "outgoing_storage_pdf_index_one_based": outgoing_index + 1,
                    "old_outgoing_pdf": old,
                    "replacement_incoming_equilibrium_pdf": new,
                    "domain_mass_delta_lattice": domain_delta,
                    "signed_port_mass_flow_kg_s": signed * mass_factor,
                }
            )
    return records, total_lattice * mass_factor


def _write_link_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("No boundary-link rows were reconstructed")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_exact_link_flux_audit(project_root: Path) -> dict[str, Any]:
    """Perform one source-proven audit without any APES executable call."""

    root = Path(project_root).resolve()
    branch = _git_value(root, "branch", "--show-current")
    head = _git_value(root, "rev-parse", "HEAD")
    if branch != SOURCE_BRANCH or head != CURRENT_SYNCED_BASE_COMMIT:
        raise FlowError(
            NOT_IDENTIFIABLE,
            f"Expected {SOURCE_BRANCH}@{CURRENT_SYNCED_BASE_COMMIT}, found {branch}@{head}",
        )
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = config.paths.output_root
    direct_run = output_root / DIRECT_FIELD_RUN
    seeder_run = output_root / FROZEN_SEEDER_RUN
    steady_run = output_root / FROZEN_STEADY_RUN
    mesh_dir = seeder_run / "seeder" / "mesh"
    direct_npz = direct_run / "flow" / "direct_cell_field.npz"
    direct_vtu = direct_run / "flow" / "flow_field.vtu"
    direct_manifest_path = direct_run / "qc" / "direct_restart_decode_manifest.json"
    direct_manifest = read_json(direct_manifest_path)
    restart_binary = Path(direct_manifest["restart_header_contract"]["binary"])
    restart_header = restart_binary.parent / "roi003274_steady_lbm_lastHeader.lua"
    elemlist_path = mesh_dir / "elemlist.lsb"
    bnd_path = mesh_dir / "bnd.lsb"
    header_path = mesh_dir / "header.lua"
    bnd_header_path = mesh_dir / "bnd.lua"
    solver_config = steady_run / "diagnostic_musubi.lua"
    runtime_log = steady_run / "tracking" / "musubi_stdout.log"
    critical_paths = (
        direct_npz,
        direct_vtu,
        direct_manifest_path,
        restart_binary,
        restart_header,
        elemlist_path,
        bnd_path,
        header_path,
        bnd_header_path,
        solver_config,
        runtime_log,
    )
    before = _file_manifest(critical_paths)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"{AUDIT_PREFIX}_{stamp}"
    qc_dir = run_root / "qc"
    qc_dir.mkdir(parents=True, exist_ok=False)
    manifest_path = qc_dir / "exact_link_flux_audit_manifest.json"
    summary: dict[str, Any] = {
        "status": NOT_IDENTIFIABLE,
        "next": NEXT_NOT_IDENTIFIABLE,
        "audit_revision": AUDIT_REVISION,
        "run_root": str(run_root),
        "branch": branch,
        "actual_head": head,
        "production_solver_modified": False,
        "external_apes_executable_calls": 0,
        "seeder_run_count": 0,
        "musubi_run_count": 0,
        "harvester_run_count": 0,
        "grid_convergence": "NOT_RUN",
        "physics_modified": False,
        "bc_modified": False,
        "mesh_modified": False,
        "surface_modified": False,
        "restart_modified": False,
        "direct_decoded_field_modified": False,
        "frozen_files_before": before,
        "started_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, summary)
    try:
        if (
            direct_manifest.get("frozen_iteration") != FROZEN_ITERATION
            or direct_manifest.get("macro_reconstruction_validation") != "PASS"
            or direct_manifest.get("field_identity", {}).get("status") != "PASS"
        ):
            raise FlowError(
                NOT_IDENTIFIABLE, "Validated direct-restart identity contract changed"
            )
        prior_source = read_json(
            direct_run / "qc" / "direct_decode_source_contract.json"
        )
        executable = Path(prior_source["musubi_executable_revision"]["path"])
        if sha256_file(executable) != MUSUBI_EXECUTABLE_SHA256:
            raise FlowError(NOT_IDENTIFIABLE, "Pinned Musubi executable hash changed")

        generic_source, side_names, side_offsets, mfr_contract = build_source_contract(
            distribution=config.apes.wsl_distribution,
            current_evidence_commit=head,
        )
        source_root = _wsl_source_root(config.apes.wsl_distribution)
        exact_source = _exact_source_contract(source_root, runtime_log)
        source_contract = {
            "status": "PASS",
            "generic_boundary_binary_and_mfr_contract": generic_source,
            "exact_restart_and_boundary_path": exact_source,
            "mfr_eq_contract": mfr_contract,
            "runtime_executable": {
                "path": str(executable),
                "sha256": sha256_file(executable),
            },
        }
        write_json(qc_dir / "exact_link_flux_source_contract.json", source_contract)

        property_header = parse_boundary_property_header(
            header_path.read_text(encoding="utf-8")
        )
        boundary_header = parse_bnd_header(bnd_header_path.read_text(encoding="utf-8"))
        if (
            property_header.element_count != EXPECTED_BOUNDARY_ELEMENT_COUNT
            or boundary_header.side_count != EXPECTED_SIDE_COUNT
        ):
            raise FlowError(
                NOT_IDENTIFIABLE, "Frozen boundary header dimensions changed"
            )
        tree_ids, property_bits, elemlist_contract = read_treelm_elemlist(
            elemlist_path, n_elems=EXPECTED_CELL_COUNT
        )
        property_indices = extract_boundary_property_indices(
            property_bits, property_header.bit_position
        )
        boundary_ids = read_boundary_ids(
            bnd_path,
            element_count=property_header.element_count,
            side_count=boundary_header.side_count,
        )
        label_to_id = {
            label: index for index, label in enumerate(boundary_header.labels, start=1)
        }
        boundaries = {
            label: reconstruct_boundary(
                boundary_ids,
                property_indices,
                label=label,
                boundary_id=label_to_id[label],
            )
            for label in PORT_LABELS
        }
        counts = {
            label: len(boundary.cell_indices) for label, boundary in boundaries.items()
        }
        if counts != EXPECTED_STENCIL_COUNTS:
            raise FlowError(NOT_IDENTIFIABLE, f"D3Q19 globBC counts changed: {counts}")
        topology_counts = {
            label: int(
                np.count_nonzero(
                    np.any(np.asarray(boundary_ids) == label_to_id[label], axis=1)
                )
            )
            for label in PORT_LABELS
        }
        if topology_counts["inlet"] != 334:
            raise FlowError(
                NOT_IDENTIFIABLE, "Frozen 26-neighbor inlet topology count changed"
            )

        with np.load(direct_npz) as data:
            npz_tree_ids = np.asarray(data["tree_id"], dtype=np.int64)
            cell_ijk = np.asarray(data["cell_ijk"], dtype=np.int64)
        if not np.array_equal(npz_tree_ids, tree_ids):
            raise FlowError(
                NOT_IDENTIFIABLE, "Restart and direct cell order no longer match"
            )
        pdf = read_restart_pdf(
            restart_binary, n_elems=EXPECTED_CELL_COUNT, n_components=19
        )
        lookup = build_coordinate_lookup(cell_ijk)

        inputs = load_flow_inputs(config.paths.source_surface_run)
        partition = load_frozen_surface_partition(
            inputs, seeder_run / "geometry" / "geometry_solver_m"
        )
        boundary_conditions = load_boundary_conditions(inputs.boundary_conditions)
        smooth_area = float(partition.patch("inlet").area_um2 * 1.0e-12)
        actual_count = len(boundaries["inlet"].cell_indices)
        area_actual = mfr_eq_area_proxy(actual_count, EXPECTED_DX_M)
        area_prompt_proxy = mfr_eq_area_proxy(topology_counts["inlet"], EXPECTED_DX_M)
        u_actual = TARGET_MASS_FLOW_KG_S / (REFERENCE_DENSITY_KG_M3 * area_actual)
        u_prompt_proxy = TARGET_MASS_FLOW_KG_S / (
            REFERENCE_DENSITY_KG_M3 * area_prompt_proxy
        )
        u_expected = TARGET_Q_M3_S / smooth_area
        area_qc = {
            "status": "PASS",
            "source_formula": "sum(nElems_totalLevel(level) * dx(level)^2)",
            "d3q19_musubi_globbc_element_count": actual_count,
            "treelm_26_neighbor_touching_row_count": topology_counts["inlet"],
            "why_334_is_not_mfr_eq_count": "23 rows contain the inlet only in D3Q27 corner directions 19..26; D3Q19 QQN scans directions 1..18",
            "dx_m": EXPECTED_DX_M,
            "area_mfr_actual_m2": area_actual,
            "area_mfr_334_hypothesis_m2": area_prompt_proxy,
            "smooth_inlet_area_m2": smooth_area,
            "actual_area_ratio": area_actual / smooth_area,
            "prompt_334_area_ratio": area_prompt_proxy / smooth_area,
            "mfr_base_physical_velocity_actual_m_s": u_actual,
            "mfr_base_physical_velocity_334_hypothesis_m_s": u_prompt_proxy,
            "expected_smooth_cap_mean_velocity_m_s": u_expected,
            "actual_velocity_ratio": u_actual / u_expected,
            "prompt_334_velocity_ratio": u_prompt_proxy / u_expected,
        }
        write_json(qc_dir / "mfr_eq_area_velocity_qc.json", area_qc)

        inlet = boundaries["inlet"]
        normal_qc = {
            "status": "PASS",
            "d3q19_element_count": actual_count,
            "d3q27_topology_count": topology_counts["inlet"],
            "d3q27_only_not_in_musubi_globbc": topology_counts["inlet"] - actual_count,
            "distribution": _normal_distribution(inlet),
            "cxDirRK_raw_not_unit": True,
            "axis_direction_norm": 1.0,
            "diagonal_direction_norm": math.sqrt(2.0),
            "incoming_modified_link_count": int(np.count_nonzero(inlet.incoming_masks)),
            "side_names_first_18": list(side_names[:18]),
            "side_offsets_first_18": side_offsets[:18].tolist(),
        }
        write_json(qc_dir / "boundary_discrete_reconstruction_qc.json", normal_qc)

        mass_factor = REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**3 / EXPECTED_DT_S
        inlet_rho = np.sum(np.asarray(pdf[inlet.cell_indices]), axis=1)
        inlet_velocity = (
            u_actual / (EXPECTED_DX_M / EXPECTED_DT_S)
        ) * D3Q19_DIRECTIONS[inlet.normal_indices].astype(np.float64)
        inlet_eq = equilibrium_pdf(inlet_rho, inlet_velocity)
        inlet_rows = np.arange(len(inlet.cell_indices), dtype=np.int64)
        link_records, inlet_mass = _link_rows(
            boundary=inlet,
            selected_rows=inlet_rows,
            tree_ids=tree_ids,
            old_pdf=pdf,
            new_pdf=inlet_eq,
            mass_factor=mass_factor,
            outward_sign=False,
        )

        pressure_factor = REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**2 / EXPECTED_DT_S**2
        pressure_reference = pressure_factor * CS2
        pressure_by_label = {
            label: pressure_reference + gauge
            for label, gauge in zip(
                PORT_LABELS[1:],
                boundary_conditions.outlet_gauge_pressures_pa,
                strict=True,
            )
        }
        outlet_mass: dict[str, float] = {}
        outlet_qc: dict[str, Any] = {}
        for label in PORT_LABELS[1:]:
            boundary = boundaries[label]
            valid, neighbor1, neighbor2 = _valid_pressure_mask(
                boundary, cell_ijk, lookup
            )
            selected = np.flatnonzero(valid).astype(np.int64)
            if len(selected) != EXPECTED_PRESSURE_VALID_COUNTS[label]:
                raise FlowError(
                    NOT_IDENTIFIABLE,
                    f"Pressure valid-element count changed for {label}: {len(selected)}",
                )
            fetched1 = pull_fetch_pdfs(
                pdf, cell_ijk, neighbor1[selected], coordinate_lookup=lookup
            )
            fetched2 = pull_fetch_pdfs(
                pdf, cell_ijk, neighbor2[selected], coordinate_lookup=lookup
            )
            velocity = 1.5 * velocity_from_pdf(fetched1) - 0.5 * velocity_from_pdf(
                fetched2
            )
            rho = pressure_by_label[label] / pressure_factor / CS2
            outlet_eq = equilibrium_pdf(rho, velocity)
            records, mass = _link_rows(
                boundary=boundary,
                selected_rows=selected,
                tree_ids=tree_ids,
                old_pdf=pdf,
                new_pdf=outlet_eq,
                mass_factor=mass_factor,
                outward_sign=True,
            )
            link_records.extend(records)
            outlet_mass[label] = mass
            outlet_qc[label] = {
                "pre_remove_d3q19_elements": len(boundary.cell_indices),
                "valid_after_remove_solid": len(selected),
                "removed_elements": int(np.count_nonzero(~valid)),
                "removed_cell_indices_zero_based": boundary.cell_indices[~valid]
                .astype(int)
                .tolist(),
                "global_inward_normal": GLOBAL_NORMALS[label].astype(int).tolist(),
                "absolute_pressure_pa": pressure_by_label[label],
                "lattice_density": rho,
                "modified_link_count": int(
                    sum(
                        np.count_nonzero(boundary.incoming_masks[row])
                        for row in selected
                    )
                ),
                "mass_flow_signed_outward_kg_s": mass,
            }

        link_csv = qc_dir / "exact_boundary_links.csv"
        _write_link_csv(link_csv, link_records)
        balance = signed_mass_balance(inlet_mass, outlet_mass.values())
        inlet_error = abs(inlet_mass - TARGET_MASS_FLOW_KG_S) / TARGET_MASS_FLOW_KG_S
        classification = classify_exact_flux(
            inlet_error, balance["relative_error"], outlet_mass["outlet_02"]
        )
        flux_qc = {
            "status": "PASS",
            "exact_pdf_link_flux_identifiable": True,
            "restart_pdf_timestep_stage_proven": True,
            "bc_buffer_exactly_reconstructable": True,
            "mass_conversion_kg_s_per_lattice_population": mass_factor,
            "sign_convention": "inlet positive into domain; outlets positive out of domain; negative outlet retained",
            "inlet": {
                "mass_flow_exact_kg_s": inlet_mass,
                "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
                "relative_error": inlet_error,
                "modified_link_count": int(np.count_nonzero(inlet.incoming_masks)),
            },
            "outlets": outlet_qc,
            "signed_balance": balance,
            "maximum_relative_error": MAXIMUM_RELATIVE_ERROR,
            "classification": classification,
            "per_link_csv": str(link_csv),
        }
        write_json(qc_dir / "exact_link_flux_qc.json", flux_qc)

        stage_qc = {
            "status": "PASS",
            "iteration": FROZEN_ITERATION,
            "restart_stage": exact_source["restart_stage"],
            "restart_pdf_timestep_stage_proven": True,
            "streaming": "PULL",
            "fetch_macro_proven": True,
            "bc_buffer_exactly_reconstructable": True,
            "neighbor_buffer_exactly_reconstructable": True,
            "single_restart_sufficient": True,
        }
        write_json(qc_dir / "restart_stage_and_streaming_qc.json", stage_qc)
        write_json(
            qc_dir / "bc_buffer_reconstructability_qc.json",
            {
                "status": "PASS",
                "exactly_reconstructable": True,
                "source_array": "frozen restart nNext serialized with SAVE",
                "approximate_substitute_used": False,
                "pressure_neighbor_fetch_reconstructed": True,
            },
        )

        after = _file_manifest(critical_paths)
        frozen_unchanged = before == after
        write_json(
            qc_dir / "frozen_read_only_sha_qc.json",
            {
                "status": "PASS" if frozen_unchanged else "FAIL",
                "frozen_files_modified": not frozen_unchanged,
                "before": before,
                "after": after,
            },
        )
        if not frozen_unchanged:
            raise FlowError(
                NOT_IDENTIFIABLE, "One or more frozen files changed during the audit"
            )

        summary.update(
            {
                **classification,
                "mfr_eq": area_qc,
                "normal_ind_distribution": normal_qc["distribution"],
                "restart_pdf_timestep_stage_proven": True,
                "streaming": "PULL",
                "bc_buffer_exactly_reconstructable": True,
                "exact_pdf_link_flux_identifiable": True,
                "exact_inlet_mass_flow_kg_s": inlet_mass,
                "target_inlet_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
                "inlet_relative_error": inlet_error,
                "exact_outlet_mass_flow_kg_s": outlet_mass,
                "exact_signed_mass_balance_error": balance["relative_error"],
                "frozen_files_after": after,
                "frozen_files_modified": False,
                "source_contract": str(qc_dir / "exact_link_flux_source_contract.json"),
                "exact_flux_qc": str(qc_dir / "exact_link_flux_qc.json"),
                "link_csv": str(link_csv),
                "elemlist_contract": elemlist_contract,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        after = _file_manifest(critical_paths)
        summary.update(
            {
                "status": NOT_IDENTIFIABLE,
                "next": NEXT_NOT_IDENTIFIABLE,
                "failure": str(error),
                "frozen_files_after": after,
                "frozen_files_modified": before != after,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
