"""Python-only audit of frozen Musubi port fluxes on the TreElm lattice."""

from __future__ import annotations

import csv
import math
import re
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .apes import load_boundary_conditions
from .config import load_cfd_flow_config
from .geometry import SurfacePartition, load_frozen_surface_partition
from .io import FlowError, load_flow_inputs, read_json, sha256_file, write_json
from .qc import evaluate_mass_conservation
from .restart_decode import read_treelm_elemlist


SOURCE_BRANCH = "codex/cfd-flow-musubi-recovery-20260828"
CURRENT_SYNCED_BASE_COMMIT = "06ba08ac63f4153146caabb55d1d2c5739cf793e"
DIRECT_DECODE_EXECUTION_BASE_COMMIT = "459a3b4c15f2d333eccbd2ce928391bd50573c93"
TREELM_SOURCE_COMMIT = "9899d1376992c4fafc8a343d2b4ccef81de670d1"
MUSUBI_SOURCE_COMMIT = "4e8b277b66226277171ef93bf054d21270812793"
MUSUBI_SCHEME_COMMIT = "81f8c4f13772f6d4af31f335e1e3f99b02726e25"

DIRECT_FIELD_RUN = "musubi_direct_restart_field_anchor003274_20260829_000544"
FROZEN_SEEDER_RUN = "musubi_recovery_anchor003274_20260828_162530"
AUDIT_PREFIX = "port_flux_audit_anchor003274"
AUDIT_REVISION = "ZERO_RUN_LATTICE_PORT_FLUX_AUDIT_V1"

PORT_LABELS = ("inlet", "outlet_01", "outlet_02", "outlet_03")
DIRECT_SIDE_NAMES = ("W", "S", "B", "E", "N", "T")
EXPECTED_SIDE_NAMES = (
    "W",
    "S",
    "B",
    "E",
    "N",
    "T",
    "BS",
    "TS",
    "BN",
    "TN",
    "BW",
    "BE",
    "TW",
    "TE",
    "SW",
    "NW",
    "SE",
    "NE",
    "BSW",
    "BSE",
    "BNW",
    "BNE",
    "TSW",
    "TSE",
    "TNW",
    "TNE",
)

EXPECTED_CELL_COUNT = 221_109
EXPECTED_BOUNDARY_ELEMENT_COUNT = 75_358
EXPECTED_SIDE_COUNT = 26
EXPECTED_BND_BYTES = EXPECTED_BOUNDARY_ELEMENT_COUNT * EXPECTED_SIDE_COUNT * 8
EXPECTED_DX_M = 2.0e-7
REFERENCE_DENSITY_KG_M3 = 1056.0
TARGET_Q_M3_S = 7.693508475538942e-16
TARGET_MASS_FLOW_KG_S = 8.124344950169123e-13

MAXIMUM_CUTSET_DEPTH = 12
MINIMUM_STABLE_DEPTH_COUNT = 3
MAXIMUM_INLET_MASS_RELATIVE_ERROR = 0.01
MAXIMUM_SIGNED_MASS_BALANCE_ERROR = 0.01
MAXIMUM_INLET_MASS_SPREAD = 0.02
MAXIMUM_OUTLET_MASS_SPREAD = 0.05

TOPOLOGY_FAILED = "CFD_FLOW_PORT_TOPOLOGY_CONTRACT_FAILED"
AUDIT_UNRESOLVED = "CFD_FLOW_PORT_FLUX_AUDIT_UNRESOLVED"
INTEGRATION_ARTIFACT = "CFD_FLOW_PORT_INTEGRATION_ARTIFACT_CONFIRMED"
BACKFLOW_CONFIRMED = "CFD_FLOW_PORT_MASS_BALANCE_PASS_OUTLET02_BACKFLOW_CONFIRMED"

NEXT_BY_STATUS = {
    TOPOLOGY_FAILED: "REVIEW PINNED TREELM BOUNDARY TOPOLOGY CONTRACT",
    AUDIT_UNRESOLVED: "SOURCE-PROVEN EXACT MFR_EQ BOUNDARY-LINK FLUX AUDIT",
    INTEGRATION_ARTIFACT: "PROMOTE LATTICE_CUTSET PORT QC THEN PROCEED TO GRID CONVERGENCE",
    BACKFLOW_CONFIRMED: "REVIEW 1D_TO_3D OUTLET PRESSURE BC MAPPING",
}


@dataclass(frozen=True, slots=True)
class BoundaryPropertyHeader:
    label: str
    bit_position: int
    element_count: int


@dataclass(frozen=True, slots=True)
class BoundaryHeader:
    side_count: int
    boundary_type_count: int
    labels: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PortFaceTopology:
    label: str
    boundary_id: int
    touching_property_rows: np.ndarray
    seed_cell_indices: np.ndarray
    face_cell_indices: np.ndarray
    face_normals: np.ndarray
    qc: dict[str, Any]


def _git_value(repository: Path, *arguments: str) -> str:
    process = subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return process.stdout.strip()


def _wsl_git_value(distribution: str, repository: str, *arguments: str) -> str:
    process = subprocess.run(
        ["wsl.exe", "-d", distribution, "--", "git", "-C", repository, *arguments],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return process.stdout.strip()


def _wsl_source_root(distribution: str) -> Path:
    candidates = (
        Path(rf"\\wsl.localhost\{distribution}\home\lzy\apes-pinned\musubi_official"),
        Path(rf"\\wsl$\{distribution}\home\lzy\apes-pinned\musubi_official"),
    )
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    raise FlowError(TOPOLOGY_FAILED, "Pinned Musubi/TreElm source tree is unavailable")


def _line_number(text: str, token: str) -> int:
    offset = text.find(token)
    if offset < 0:
        raise FlowError(TOPOLOGY_FAILED, f"Pinned source token is missing: {token}")
    return text.count("\n", 0, offset) + 1


def _source_evidence(
    path: Path,
    revision: str,
    statements: dict[str, tuple[str, ...]],
) -> dict[str, Any]:
    text = Path(path).read_text(encoding="utf-8", errors="strict")
    evidence: dict[str, Any] = {}
    for label, tokens in statements.items():
        evidence[label] = {
            "status": "PASS",
            "line_numbers": [_line_number(text, token) for token in tokens],
            "tokens": list(tokens),
        }
    return {
        "path": str(path),
        "source_revision": revision,
        "sha256": sha256_file(path),
        "evidence": evidence,
    }


def parse_boundary_property_header(
    header_text: str,
    *,
    property_label: str = "has boundaries",
) -> BoundaryPropertyHeader:
    property_match = re.search(r"property\s*=\s*\{(.*)\}\s*$", header_text, re.DOTALL)
    if not property_match:
        raise ValueError("Mesh header does not contain a property table")
    blocks = re.findall(r"\{(.*?)\}", property_match.group(1), re.DOTALL)
    for block in blocks:
        label_match = re.search(r"label\s*=\s*['\"]([^'\"]+)['\"]", block)
        if not label_match or label_match.group(1).strip() != property_label:
            continue
        bit_match = re.search(r"bitpos\s*=\s*(\d+)", block)
        count_match = re.search(r"nElems\s*=\s*(\d+)", block)
        if not bit_match or not count_match:
            raise ValueError(f"Incomplete property header for {property_label}")
        return BoundaryPropertyHeader(
            label=property_label,
            bit_position=int(bit_match.group(1)),
            element_count=int(count_match.group(1)),
        )
    raise ValueError(f"Property label is missing: {property_label}")


def extract_boundary_property_indices(property_bits: np.ndarray, bit_position: int) -> np.ndarray:
    """Mirror Fortran BTEST and gather_property's ascending element-order scan."""

    bits = np.asarray(property_bits, dtype=np.int64).reshape(-1)
    mask = (bits & (np.int64(1) << np.int64(bit_position))) != 0
    return np.flatnonzero(mask).astype(np.int64)


def parse_bnd_header(text: str) -> BoundaryHeader:
    side_match = re.search(r"\bnSides\s*=\s*(\d+)", text)
    type_match = re.search(r"\bnBCtypes\s*=\s*(\d+)", text)
    labels_match = re.search(r"\bbclabel\s*=\s*\{(.*?)\}", text, re.DOTALL)
    if not side_match or not type_match or not labels_match:
        raise ValueError("Incomplete bnd.lua header")
    labels = tuple(re.findall(r"['\"]([^'\"]+)['\"]", labels_match.group(1)))
    header = BoundaryHeader(int(side_match.group(1)), int(type_match.group(1)), labels)
    if len(header.labels) != header.boundary_type_count:
        raise ValueError("bnd.lua boundary label count does not match nBCtypes")
    return header


def parse_treelm_side_contract(source_text: str) -> tuple[tuple[str, ...], np.ndarray]:
    uncommented = "\n".join(line.split("!", 1)[0] for line in source_text.splitlines())
    normalized = uncommented.replace("&", " ")
    names_match = re.search(
        r"qDirName\s*\(\s*qQQQ\s*\)\s*=\s*\[(.*?)\]",
        normalized,
        re.DOTALL,
    )
    offsets_match = re.search(
        r"qOffset\s*=\s*reshape\s*\(\s*\(/(.*?)/\)\s*,\s*\(/\s*qQQQ\s*,\s*3\s*/\)\s*\)",
        normalized,
        re.DOTALL,
    )
    if not names_match or not offsets_match:
        raise ValueError("Could not parse pinned qDirName/qOffset definitions")
    names = tuple(value.strip() for value in re.findall(r"['\"]([^'\"]+)['\"]", names_match.group(1)))
    values = [int(value) for value in re.findall(r"(?<![A-Za-z0-9_])[-+]?\d+", offsets_match.group(1))]
    if len(names) != EXPECTED_SIDE_COUNT or len(values) != EXPECTED_SIDE_COUNT * 3:
        raise ValueError("Pinned side contract is not 26 directions by 3 components")
    offsets = np.asarray(values, dtype=np.int8).reshape((EXPECTED_SIDE_COUNT, 3), order="F")
    return names, offsets


def boundary_binary_contract(
    path: Path,
    *,
    element_count: int,
    side_count: int,
) -> dict[str, Any]:
    expected = int(element_count) * int(side_count) * np.dtype("<i8").itemsize
    actual = Path(path).stat().st_size
    return {
        "status": "PASS" if actual == expected else "FAIL",
        "path": str(path),
        "actual_bytes": actual,
        "expected_bytes": expected,
        "dtype": "<i8",
        "layout": "element-major rows; 26 contiguous boundary_ID values per property element",
        "element_count": int(element_count),
        "side_count": int(side_count),
    }


def read_boundary_ids(path: Path, *, element_count: int, side_count: int) -> np.memmap:
    contract = boundary_binary_contract(path, element_count=element_count, side_count=side_count)
    if contract["status"] != "PASS":
        raise FlowError(TOPOLOGY_FAILED, f"bnd.lsb byte contract failed: {contract}")
    return np.memmap(path, dtype="<i8", mode="r", shape=(element_count, side_count), order="C")


def build_source_contract(
    *,
    distribution: str,
    current_evidence_commit: str,
) -> tuple[dict[str, Any], tuple[str, ...], np.ndarray, dict[str, Any]]:
    root = _wsl_source_root(distribution)
    posix_root = "/home/lzy/apes-pinned/musubi_official"
    repositories = {
        "musubi": (posix_root, MUSUBI_SOURCE_COMMIT),
        "musubi_scheme": (f"{posix_root}/mus", MUSUBI_SCHEME_COMMIT),
        "treelm": (f"{posix_root}/tem", TREELM_SOURCE_COMMIT),
    }
    revision_qc: dict[str, Any] = {}
    for label, (repository, expected) in repositories.items():
        actual = _wsl_git_value(distribution, repository, "rev-parse", "HEAD")
        dirty = _wsl_git_value(
            distribution,
            repository,
            "status",
            "--short",
            "--untracked-files=no",
        )
        if actual != expected or dirty:
            raise FlowError(
                TOPOLOGY_FAILED,
                f"Pinned {label} revision/cleanliness mismatch: {actual}, status={dirty!r}",
            )
        revision_qc[label] = {
            "path": repository,
            "commit": actual,
            "expected_commit": expected,
            "tracked_worktree_clean": True,
        }

    tem = root / "tem" / "source"
    mus = root / "mus" / "source"
    property_source = tem / "tem_property_module.f90"
    boundary_source = tem / "tem_bc_prop_module.f90"
    parameter_source = tem / "tem_param_module.f90"
    environment_source = tem / "env_module.f90"
    mfr_source = mus / "bc" / "mus_bc_fluid_module.fpp"

    side_names, side_offsets = parse_treelm_side_contract(
        parameter_source.read_text(encoding="utf-8", errors="strict")
    )
    if side_names != EXPECTED_SIDE_NAMES:
        raise FlowError(TOPOLOGY_FAILED, f"Pinned 26-side ordering changed: {side_names}")
    expected_direct = np.asarray(
        ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
        dtype=np.int8,
    )
    if not np.array_equal(side_offsets[:6], expected_direct):
        raise FlowError(TOPOLOGY_FAILED, "Pinned direct-face qOffset ordering changed")

    source_files = {
        "property_element_mapping": _source_evidence(
            property_source,
            TREELM_SOURCE_COMMIT,
            {
                "btest_count": ("Property%nElems = count(btest(BitField, Header%BitPos))",),
                "ascending_tree_element_scan": (
                    "do iElem=1,nElems",
                    "if (btest(BitField(iElem), Header%BitPos)) then",
                    "Property%ElemID(PropElem) = iElem",
                ),
            },
        ),
        "boundary_binary_layout": _source_evidence(
            boundary_source,
            TREELM_SOURCE_COMMIT,
            {
                "contiguous_sides_per_element": (
                    "call MPI_TYPE_CONTIGUOUS( me%nSides, long_k_mpi, etype, iError )",
                    "call MPI_FILE_READ_ALL( fh, buffer, nElems, etype, iostatus, iError )",
                ),
                "row_mapping": (
                    "me%boundary_ID(:,i) = buffer( ((i-1)*me%nSides+1) : (i*me%nSides) )",
                ),
                "dump_matches_load": (
                    "call MPI_FILE_WRITE_ALL( fh, me%boundary_ID, nElems, etype, iostatus, iError )",
                ),
            },
        ),
        "side_ordering": _source_evidence(
            parameter_source,
            TREELM_SOURCE_COMMIT,
            {
                "direct_face_indices": (
                    "q__W     = 1",
                    "q__S     = 2",
                    "q__B     = 3",
                    "q__E     = 4",
                    "q__N     = 5",
                    "q__T     = 6",
                ),
                "offset_definition": ("qOffset =", "reshape((/-1, 0, 0, 1, 0, 0"),
            },
        ),
        "integer_width_and_endian": _source_evidence(
            environment_source,
            TREELM_SOURCE_COMMIT,
            {
                "integer8": ("long_k_mpi = MPI_INTEGER8",),
                "little_endian_suffix": (
                    "isLittleEndian = Sys_is_Little_Endian()",
                    "suffix = '.lsb'",
                ),
            },
        ),
    }

    mfr_evidence = _source_evidence(
        mfr_source,
        MUSUBI_SCHEME_COMMIT,
        {
            "mfr_eq_routine": ("subroutine mfr_eq(", "res     = massFlowRate"),
            "total_boundary_area": (
                "area = area + globBC%nElems_totalLevel(iLvl) * physics%dxLvl(iLvl)**2",
            ),
            "mass_to_physical_velocity": (
                "massFlowRateToVel = 1.0_rk / ( physics%rho0 * area )",
            ),
            "physical_to_lattice_velocity": ("/ physics%fac( iLevel )%vel",),
            "normal_velocity_per_boundary_element": (
                "velocity(:,iElem) = massFlowRate(iElem) * massFlowRateToVel",
                "normalInd%val(iElem)",
            ),
        },
    )
    mfr_contract = {
        "status": "PASS",
        "source": mfr_evidence,
        "mass_flowrate_input_unit": "kg/s",
        "unit_basis": "massFlowRate/(rho0[kg/m3] * area[m2]) produces physical velocity[m/s]",
        "physical_to_lattice_conversion": "divide physical velocity by physics%fac(iLevel)%vel",
        "boundary_velocity": "normal-aligned velocity using cxDirRK(:, normalInd)",
        "boundary_area": "sum nElems_totalLevel(level) * dx(level)^2 over levels",
        "target_distribution": (
            "the same configured mass-flow value produces one uniform normal speed from the total "
            "boundary-element area; equilibrium PDFs are then assigned per boundary element"
        ),
        "implementation_modified": False,
    }
    contract = {
        "status": "PASS",
        "audit_revision": AUDIT_REVISION,
        "direct_decode_execution_base_commit": DIRECT_DECODE_EXECUTION_BASE_COMMIT,
        "current_evidence_commit": current_evidence_commit,
        "pinned_source_root": str(root),
        "git_pull_performed": False,
        "submodule_update_performed": False,
        "metadata_git_calls_only": True,
        "repositories": revision_qc,
        "source_files": source_files,
        "property_element_mapping": {
            "status": "PASS",
            "fortran_bit_test": "BTEST(ElemPropertyBits, bitpos)",
            "python_mirror": "(bits & (1 << bitpos)) != 0",
            "ordering": "ascending tree-element order",
        },
        "boundary_binary": {
            "status": "PASS",
            "dtype": "little-endian signed int64 (<i8)",
            "ordering": "26 contiguous side values per boundary-property element",
        },
        "side_ordering": {
            "status": "PASS",
            "names": list(side_names),
            "q_offset": side_offsets.tolist(),
            "direct_six": list(DIRECT_SIDE_NAMES),
        },
    }
    return contract, side_names, side_offsets, mfr_contract


def legacy_flux_semantics(legacy_qc: dict[str, Any]) -> dict[str, Any]:
    q_in = float(legacy_qc["q_in_m3_s"])
    q_out = tuple(float(legacy_qc["q_out_m3_s"][label]) for label in PORT_LABELS[1:])
    signed = evaluate_mass_conservation(q_in, q_out)
    absolute_sum = float(sum(abs(value) for value in q_out))
    magnitude_mismatch = abs(abs(q_in) - absolute_sum) / abs(q_in)
    return {
        "status": "PASS",
        "estimator": "LEGACY_SMOOTH_CAP_INTERPOLATED_ESTIMATOR",
        "algorithm": [
            "CellData to point data",
            "translate smooth VMTK cap inward by 2dx",
            "three barycentric samples per triangle",
            "VTK interpolation",
            "integrate u dot outward-normal",
        ],
        "conservation_property": "NOT_CONSERVATIVE_BY_CONSTRUCTION",
        "conservation_property_interpretation": (
            "no discrete conservation guarantee; this statement alone does not prove a code bug"
        ),
        "q_in_m3_s": q_in,
        "q_out_signed_m3_s": {
            label: value for label, value in zip(PORT_LABELS[1:], q_out, strict=True)
        },
        "signed_q_out_sum_m3_s": signed["q_out_sum_m3_s"],
        "signed_volumetric_balance_error": signed["relative_error"],
        "sum_absolute_outlet_flow_m3_s": absolute_sum,
        "absolute_flow_magnitude_mismatch": magnitude_mismatch,
        "absolute_magnitude_role": "ABSOLUTE_MAGNITUDE_DIAGNOSTIC",
        "old_0p56_value_is_mass_conservation_error": False,
    }


def _patch_reference(partition: SurfacePartition, label: str) -> tuple[float, np.ndarray]:
    patch = partition.patch(label)
    return float(patch.area_um2 * 1.0e-12), np.asarray(patch.outward_normal, dtype=np.float64)


def build_port_topology(
    *,
    boundary_ids: np.ndarray,
    property_element_indices: np.ndarray,
    boundary_header: BoundaryHeader,
    side_names: tuple[str, ...],
    side_offsets: np.ndarray,
    partition: SurfacePartition,
    dx_m: float,
) -> tuple[dict[str, Any], dict[str, PortFaceTopology]]:
    ids = np.asarray(boundary_ids, dtype=np.int64)
    prop_indices = np.asarray(property_element_indices, dtype=np.int64).reshape(-1)
    if ids.shape != (len(prop_indices), len(side_names)):
        raise FlowError(TOPOLOGY_FAILED, "bnd rows do not match boundary-property element mapping")
    label_to_id = {label: index for index, label in enumerate(boundary_header.labels, start=1)}
    if not set(PORT_LABELS).issubset(label_to_id):
        raise FlowError(TOPOLOGY_FAILED, f"Frozen bnd labels do not contain four ports: {label_to_id}")

    port_results: dict[str, Any] = {}
    topology: dict[str, PortFaceTopology] = {}
    all_pass = True
    for label in PORT_LABELS:
        boundary_id = label_to_id[label]
        touching_rows = np.flatnonzero(np.any(ids == boundary_id, axis=1)).astype(np.int64)
        direct_rows, direct_columns = np.nonzero(ids[:, :6] == boundary_id)
        diagonal_count = int(np.count_nonzero(ids[:, 6:] == boundary_id))
        face_cell_indices = prop_indices[direct_rows]
        face_normals = np.asarray(side_offsets[direct_columns], dtype=np.int8)
        seed_indices = np.unique(face_cell_indices)
        area_vector = np.sum(face_normals.astype(np.float64), axis=0) * float(dx_m) ** 2
        area_magnitude = float(np.linalg.norm(area_vector))
        smooth_area, smooth_normal = _patch_reference(partition, label)
        dot_product = float(np.dot(area_vector, smooth_normal))
        cosine = dot_product / area_magnitude if area_magnitude > 0.0 else math.nan
        angle = float(math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))) if np.isfinite(cosine) else None
        projected_area = dot_product
        projected_difference = abs(projected_area - smooth_area) / smooth_area
        orientation_counts = {
            name: int(np.count_nonzero(direct_columns == index))
            for index, name in enumerate(DIRECT_SIDE_NAMES)
        }
        port_pass = len(direct_rows) > 0 and area_magnitude > 0.0 and dot_product > 0.0
        all_pass = all_pass and port_pass
        port_qc = {
            "status": "PASS" if port_pass else "FAIL",
            "boundary_id": boundary_id,
            "boundary_elements_touching_port": int(len(touching_rows)),
            "direct_face_count": int(len(direct_rows)),
            "diagonal_link_count": diagonal_count,
            "direct_face_orientation_counts": orientation_counts,
            "lattice_area_vector_m2": area_vector.tolist(),
            "lattice_area_vector_magnitude_m2": area_magnitude,
            "smooth_cap_area_m2": smooth_area,
            "smooth_cap_outward_normal": smooth_normal.tolist(),
            "dot_product": dot_product,
            "cosine_alignment": float(cosine),
            "angle_deg": angle,
            "projected_area_on_smooth_normal_m2": projected_area,
            "relative_projected_area_difference": projected_difference,
            "area_difference_role": "TOPOLOGY_DIAGNOSTIC_ONLY",
            "seed_cell_count": int(len(seed_indices)),
        }
        port_results[label] = port_qc
        topology[label] = PortFaceTopology(
            label=label,
            boundary_id=boundary_id,
            touching_property_rows=touching_rows,
            seed_cell_indices=seed_indices,
            face_cell_indices=np.asarray(face_cell_indices, dtype=np.int64),
            face_normals=face_normals,
            qc=port_qc,
        )
    return {
        "status": "PASS" if all_pass else "FAIL",
        "hard_gates": [
            "direct face count > 0",
            "lattice area-vector magnitude > 0",
            "dot(lattice area vector, smooth outward normal) > 0",
        ],
        "ports": port_results,
    }, topology


def boundary_cell_fluxes(
    topology: dict[str, PortFaceTopology],
    velocity_m_s: np.ndarray,
    density_lattice: np.ndarray,
    *,
    dx_m: float,
    rho0_kg_m3: float,
) -> dict[str, Any]:
    velocity = np.asarray(velocity_m_s, dtype=np.float64)
    density = np.asarray(density_lattice, dtype=np.float64).reshape(-1) * float(rho0_kg_m3)
    ports: dict[str, Any] = {}
    for label in PORT_LABELS:
        item = topology[label]
        indices = item.face_cell_indices
        normals = item.face_normals.astype(np.float64)
        axial_velocity = np.einsum("ij,ij->i", velocity[indices], normals)
        q_outward = float(np.sum(axial_velocity) * dx_m**2)
        mass_outward = float(np.sum(density[indices] * axial_velocity) * dx_m**2)
        sign = -1.0 if label == "inlet" else 1.0
        ports[label] = {
            "q_outward_m3_s": q_outward,
            "mass_outward_kg_s": mass_outward,
            "q_m3_s": sign * q_outward,
            "mass_flow_kg_s": sign * mass_outward,
            "face_count": int(len(indices)),
        }
    return {
        "status": "DIAGNOSTIC",
        "estimator": "LATTICE_BOUNDARY_CELL_CENTER_ESTIMATOR",
        "final_gate_role": "DIAGNOSTIC_ONLY_NOT_EXACT_BOUNDARY_FACE_VELOCITY",
        "ports": ports,
    }


def build_face_neighbor_graph(cell_ijk: np.ndarray) -> np.ndarray:
    coordinates = np.asarray(cell_ijk, dtype=np.int64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("cell_ijk must have shape (n, 3)")
    lookup = {tuple(row): index for index, row in enumerate(coordinates)}
    if len(lookup) != len(coordinates):
        raise ValueError("cell_ijk contains duplicate cells")
    direct = np.asarray(
        ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
        dtype=np.int64,
    )
    graph = np.full((len(coordinates), 6), -1, dtype=np.int64)
    for index, coordinate in enumerate(coordinates):
        for direction, offset in enumerate(direct):
            graph[index, direction] = lookup.get(tuple(coordinate + offset), -1)
    return graph


def _region_overlap(regions: dict[str, set[int]]) -> tuple[bool, list[dict[str, Any]]]:
    details: list[dict[str, Any]] = []
    labels = tuple(regions)
    for left_index, left in enumerate(labels):
        for right in labels[left_index + 1 :]:
            overlap = regions[left] & regions[right]
            if overlap:
                details.append({"ports": [left, right], "cell_count": len(overlap)})
    return bool(details), details


def _cutset_flux(
    *,
    port_indices: np.ndarray,
    union_mask: np.ndarray,
    graph: np.ndarray,
    cell_ijk: np.ndarray,
    velocity_m_s: np.ndarray,
    mass_flux_vector: np.ndarray,
    dx_m: float,
) -> dict[str, Any]:
    neighbor_matrix = graph[port_indices]
    valid = neighbor_matrix >= 0
    safe_neighbors = np.where(valid, neighbor_matrix, 0)
    cut_mask = valid & ~union_mask[safe_neighbors]
    repeated_port = np.repeat(port_indices, graph.shape[1])
    flat_mask = cut_mask.reshape(-1)
    core_indices = neighbor_matrix.reshape(-1)[flat_mask]
    cut_port_indices = repeated_port[flat_mask]
    normals = (
        np.asarray(cell_ijk[cut_port_indices], dtype=np.int64)
        - np.asarray(cell_ijk[core_indices], dtype=np.int64)
    ).astype(np.float64)
    if len(normals) and not np.all(np.sum(np.abs(normals), axis=1) == 1):
        raise FlowError(TOPOLOGY_FAILED, "Cutset contains a non-face neighbor")
    u_face = 0.5 * (velocity_m_s[core_indices] + velocity_m_s[cut_port_indices])
    j_face = 0.5 * (mass_flux_vector[core_indices] + mass_flux_vector[cut_port_indices])
    q_outward = float(np.sum(np.einsum("ij,ij->i", u_face, normals)) * dx_m**2)
    mass_outward = float(np.sum(np.einsum("ij,ij->i", j_face, normals)) * dx_m**2)
    return {
        "cut_face_count": int(len(core_indices)),
        "q_outward_m3_s": q_outward,
        "mass_outward_kg_s": mass_outward,
    }


def lattice_internal_cutset_sweep(
    *,
    graph: np.ndarray,
    cell_ijk: np.ndarray,
    velocity_m_s: np.ndarray,
    density_lattice: np.ndarray,
    seeds_by_port: dict[str, np.ndarray],
    dx_m: float,
    rho0_kg_m3: float,
    maximum_depth: int,
) -> tuple[list[dict[str, Any]], int | None]:
    """Evaluate fixed BFS depths; depth one contains boundary seed cells only."""

    velocity = np.asarray(velocity_m_s, dtype=np.float64)
    physical_density = np.asarray(density_lattice, dtype=np.float64).reshape(-1) * rho0_kg_m3
    mass_flux = velocity * physical_density[:, None]
    visited = {label: set(int(value) for value in seeds) for label, seeds in seeds_by_port.items()}
    frontier = {label: set(values) for label, values in visited.items()}
    rows: list[dict[str, Any]] = []
    first_overlap: int | None = None
    for depth in range(1, maximum_depth + 1):
        overlap, overlap_details = _region_overlap(visited)
        if overlap and first_overlap is None:
            first_overlap = depth
        union_mask = np.zeros(len(graph), dtype=bool)
        for values in visited.values():
            if values:
                union_mask[np.fromiter(values, dtype=np.int64)] = True
        for label, region in visited.items():
            region_indices = np.fromiter(region, dtype=np.int64)
            flux = _cutset_flux(
                port_indices=region_indices,
                union_mask=union_mask,
                graph=graph,
                cell_ijk=cell_ijk,
                velocity_m_s=velocity,
                mass_flux_vector=mass_flux,
                dx_m=dx_m,
            )
            sign = -1.0 if label == "inlet" else 1.0
            rows.append(
                {
                    "depth_cells": depth,
                    "depth_um": depth * dx_m * 1.0e6,
                    "port": label,
                    "visited_cell_count": len(region),
                    "cut_face_count": flux["cut_face_count"],
                    "q_outward_m3_s": flux["q_outward_m3_s"],
                    "mass_outward_kg_s": flux["mass_outward_kg_s"],
                    "Q_m3_s": sign * flux["q_outward_m3_s"],
                    "Mdot_kg_s": sign * flux["mass_outward_kg_s"],
                    "cross_port_overlap": overlap,
                    "overlap_details": overlap_details,
                    "global_balance_valid": first_overlap is None,
                }
            )
        if depth == maximum_depth:
            break
        next_frontier: dict[str, set[int]] = {}
        for label, current_frontier in frontier.items():
            candidates: set[int] = set()
            for index in current_frontier:
                candidates.update(int(value) for value in graph[index] if value >= 0)
            candidates.difference_update(visited[label])
            visited[label].update(candidates)
            next_frontier[label] = candidates
        frontier = next_frontier
    return rows, first_overlap


def signed_balance(q_in: float, q_out: Iterable[float]) -> dict[str, Any]:
    """Use the production helper's signed-outlet conservation semantics."""

    return evaluate_mass_conservation(float(q_in), tuple(float(value) for value in q_out))


def absolute_magnitude_diagnostic(q_in: float, q_out: Iterable[float]) -> float:
    values = tuple(float(value) for value in q_out)
    return abs(abs(float(q_in)) - sum(abs(value) for value in values)) / abs(float(q_in))


def summarize_cutset_depths(
    rows: list[dict[str, Any]],
    *,
    q_target_m3_s: float,
    mass_target_kg_s: float,
) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    depths = sorted({int(row["depth_cells"]) for row in rows})
    for depth in depths:
        per_port = {row["port"]: row for row in rows if int(row["depth_cells"]) == depth}
        missing = set(PORT_LABELS) - set(per_port)
        if missing:
            raise ValueError(f"Depth {depth} lacks ports: {sorted(missing)}")
        q_in = float(per_port["inlet"]["Q_m3_s"])
        mass_in = float(per_port["inlet"]["Mdot_kg_s"])
        q_out = tuple(float(per_port[label]["Q_m3_s"]) for label in PORT_LABELS[1:])
        mass_out = tuple(float(per_port[label]["Mdot_kg_s"]) for label in PORT_LABELS[1:])
        q_balance = signed_balance(q_in, q_out)
        mass_balance = signed_balance(mass_in, mass_out)
        overlap = any(bool(item["cross_port_overlap"]) for item in per_port.values())
        valid = all(bool(item["global_balance_valid"]) for item in per_port.values()) and not overlap
        outlet_two = mass_out[1]
        summaries.append(
            {
                "depth_cells": depth,
                "depth_um": float(per_port["inlet"]["depth_um"]),
                "global_balance_valid": valid,
                "cross_port_overlap": overlap,
                "q_in_m3_s": q_in,
                "q_out_m3_s": dict(zip(PORT_LABELS[1:], q_out, strict=True)),
                "mass_in_kg_s": mass_in,
                "mass_out_kg_s": dict(zip(PORT_LABELS[1:], mass_out, strict=True)),
                "inlet_q_relative_error": abs(q_in - q_target_m3_s) / abs(q_target_m3_s),
                "inlet_mass_flow_relative_error": abs(mass_in - mass_target_kg_s)
                / abs(mass_target_kg_s),
                "signed_volumetric_balance_error": q_balance["relative_error"],
                "signed_mass_balance_error": mass_balance["relative_error"],
                "sum_absolute_outlet_flow_m3_s": sum(abs(value) for value in q_out),
                "absolute_magnitude_role": "ABSOLUTE_MAGNITUDE_DIAGNOSTIC",
                "outlet_02_sign": (
                    "NEGATIVE" if outlet_two < 0.0 else "POSITIVE" if outlet_two > 0.0 else "ZERO"
                ),
            }
        )
    return summaries


def _relative_spread(values: Iterable[float]) -> float:
    array = np.asarray(tuple(values), dtype=np.float64)
    mean = float(np.mean(array))
    if mean == 0.0:
        return math.inf
    return float((np.max(array) - np.min(array)) / abs(mean))


def find_stable_window(
    depth_summaries: list[dict[str, Any]],
    *,
    minimum_depth_count: int = MINIMUM_STABLE_DEPTH_COUNT,
) -> dict[str, Any]:
    """Select the first shallow-to-deep qualifying triple, then extend it."""

    ordered = sorted(depth_summaries, key=lambda item: int(item["depth_cells"]))
    by_depth = {int(item["depth_cells"]): item for item in ordered}

    def depth_gate(item: dict[str, Any]) -> bool:
        return bool(
            item["global_balance_valid"]
            and not item["cross_port_overlap"]
            and item["inlet_mass_flow_relative_error"] <= MAXIMUM_INLET_MASS_RELATIVE_ERROR
            and item["signed_mass_balance_error"] <= MAXIMUM_SIGNED_MASS_BALANCE_ERROR
        )

    def stability_gate(window: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
        inlet_spread = _relative_spread(item["mass_in_kg_s"] for item in window)
        outlet_spreads = {
            label: _relative_spread(item["mass_out_kg_s"][label] for item in window)
            for label in PORT_LABELS[1:]
        }
        passed = inlet_spread <= MAXIMUM_INLET_MASS_SPREAD and all(
            value <= MAXIMUM_OUTLET_MASS_SPREAD for value in outlet_spreads.values()
        )
        return passed, {
            "inlet_mass_flow_relative_spread": inlet_spread,
            "outlet_mass_flow_relative_spread": outlet_spreads,
        }

    depths = sorted(by_depth)
    selected: list[dict[str, Any]] | None = None
    stability: dict[str, Any] | None = None
    for start in depths:
        candidate_depths = list(range(start, start + minimum_depth_count))
        if any(depth not in by_depth for depth in candidate_depths):
            continue
        candidate = [by_depth[depth] for depth in candidate_depths]
        stable, candidate_stability = stability_gate(candidate)
        if not all(depth_gate(item) for item in candidate) or not stable:
            continue
        selected = candidate
        stability = candidate_stability
        next_depth = candidate_depths[-1] + 1
        while next_depth in by_depth:
            extended = [*selected, by_depth[next_depth]]
            extended_stable, extended_metrics = stability_gate(extended)
            if not depth_gate(by_depth[next_depth]) or not extended_stable:
                break
            selected = extended
            stability = extended_metrics
            next_depth += 1
        break

    policy = {
        "selection": "first shallow-to-deep qualifying 3-depth window, extended while all gates hold",
        "minimum_consecutive_depths": minimum_depth_count,
        "cross_port_overlap_allowed": False,
        "maximum_inlet_mass_flow_relative_error": MAXIMUM_INLET_MASS_RELATIVE_ERROR,
        "maximum_signed_mass_balance_error": MAXIMUM_SIGNED_MASS_BALANCE_ERROR,
        "maximum_inlet_mass_flow_relative_spread": MAXIMUM_INLET_MASS_SPREAD,
        "maximum_each_outlet_mass_flow_relative_spread": MAXIMUM_OUTLET_MASS_SPREAD,
    }
    if selected is None:
        return {"found": False, "status": "NO_STABLE_WINDOW", "policy": policy}

    mean_q = {
        label: float(np.mean([item["q_in_m3_s"] if label == "inlet" else item["q_out_m3_s"][label] for item in selected]))
        for label in PORT_LABELS
    }
    mean_mass = {
        label: float(np.mean([item["mass_in_kg_s"] if label == "inlet" else item["mass_out_kg_s"][label] for item in selected]))
        for label in PORT_LABELS
    }
    q_balance = signed_balance(mean_q["inlet"], (mean_q[label] for label in PORT_LABELS[1:]))
    mass_balance = signed_balance(
        mean_mass["inlet"],
        (mean_mass[label] for label in PORT_LABELS[1:]),
    )
    outlet_two_values = [item["mass_out_kg_s"]["outlet_02"] for item in selected]
    return {
        "found": True,
        "status": "PORT_CUTSET_STABLE_WINDOW_FOUND",
        "policy": policy,
        "depth_start": int(selected[0]["depth_cells"]),
        "depth_end": int(selected[-1]["depth_cells"]),
        "depth_count": len(selected),
        "depths": [int(item["depth_cells"]) for item in selected],
        "stability": stability,
        "mean_q_m3_s": mean_q,
        "mean_mass_flow_kg_s": mean_mass,
        "inlet_q_relative_error": abs(mean_q["inlet"] - TARGET_Q_M3_S) / abs(TARGET_Q_M3_S),
        "inlet_mass_flow_relative_error": abs(mean_mass["inlet"] - TARGET_MASS_FLOW_KG_S)
        / abs(TARGET_MASS_FLOW_KG_S),
        "signed_volumetric_balance_error": q_balance["relative_error"],
        "signed_mass_balance_error": mass_balance["relative_error"],
        "outlet_02_mass_flow_values_kg_s": outlet_two_values,
        "outlet_02_all_negative": all(value < 0.0 for value in outlet_two_values),
        "outlet_02_all_positive": all(value > 0.0 for value in outlet_two_values),
    }


def classify_audit(
    *,
    topology_pass: bool,
    stable_window: dict[str, Any],
    legacy_outlet_02_m3_s: float,
) -> dict[str, Any]:
    if not topology_pass:
        status = TOPOLOGY_FAILED
        backflow = "UNRESOLVED"
        legacy_artifact = "UNRESOLVED"
    elif not stable_window.get("found", False):
        status = AUDIT_UNRESOLVED
        backflow = "UNRESOLVED"
        legacy_artifact = "UNRESOLVED"
    elif stable_window["outlet_02_all_negative"]:
        status = BACKFLOW_CONFIRMED
        backflow = "YES"
        legacy_artifact = "NO"
    else:
        status = INTEGRATION_ARTIFACT
        backflow = "NO" if stable_window["outlet_02_all_positive"] else "UNRESOLVED"
        legacy_artifact = (
            "YES"
            if stable_window["outlet_02_all_positive"] and legacy_outlet_02_m3_s < 0.0
            else "UNRESOLVED"
        )
    return {
        "status": status,
        "next": NEXT_BY_STATUS[status],
        "outlet_02_backflow_confirmed": backflow,
        "legacy_geometric_estimator_classified_as_artifact": legacy_artifact,
    }


def _write_depth_csv(
    path: Path,
    rows: list[dict[str, Any]],
    summaries: list[dict[str, Any]],
) -> None:
    summary_by_depth = {int(item["depth_cells"]): item for item in summaries}
    fields = (
        "depth_cells",
        "depth_um",
        "port",
        "visited_cell_count",
        "cut_face_count",
        "Q_m3_s",
        "Mdot_kg_s",
        "q_outward_m3_s",
        "mass_outward_kg_s",
        "inlet_q_relative_error",
        "inlet_mass_flow_relative_error",
        "signed_volumetric_balance_error",
        "signed_mass_balance_error",
        "outlet_02_sign",
        "cross_port_overlap",
        "global_balance_valid",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            summary = summary_by_depth[int(row["depth_cells"])]
            writer.writerow(
                {
                    **{field: row[field] for field in fields[:9]},
                    "inlet_q_relative_error": summary["inlet_q_relative_error"],
                    "inlet_mass_flow_relative_error": summary[
                        "inlet_mass_flow_relative_error"
                    ],
                    "signed_volumetric_balance_error": summary[
                        "signed_volumetric_balance_error"
                    ],
                    "signed_mass_balance_error": summary["signed_mass_balance_error"],
                    "outlet_02_sign": summary["outlet_02_sign"],
                    "cross_port_overlap": row["cross_port_overlap"],
                    "global_balance_valid": row["global_balance_valid"],
                }
            )


def _file_manifest(paths: Iterable[Path]) -> dict[str, Any]:
    return {
        str(path.resolve()): {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for path in paths
    }


def run_port_flux_audit(project_root: Path) -> dict[str, Any]:
    """Run the fixed-depth audit without launching Seeder, Musubi, or Harvester."""

    root = Path(project_root).resolve()
    branch = _git_value(root, "branch", "--show-current")
    head = _git_value(root, "rev-parse", "HEAD")
    if branch != SOURCE_BRANCH or head != CURRENT_SYNCED_BASE_COMMIT:
        raise FlowError(
            AUDIT_UNRESOLVED,
            f"Expected {SOURCE_BRANCH}@{CURRENT_SYNCED_BASE_COMMIT}, found {branch}@{head}",
        )
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = config.paths.output_root
    direct_run = output_root / DIRECT_FIELD_RUN
    seeder_run = output_root / FROZEN_SEEDER_RUN
    mesh_dir = seeder_run / "seeder" / "mesh"
    direct_npz = direct_run / "flow" / "direct_cell_field.npz"
    direct_vtu = direct_run / "flow" / "flow_field.vtu"
    direct_manifest = read_json(direct_run / "qc" / "direct_restart_decode_manifest.json")
    restart_binary = Path(direct_manifest["restart_header_contract"]["binary"])
    elemlist_path = mesh_dir / "elemlist.lsb"
    bnd_path = mesh_dir / "bnd.lsb"
    critical_paths = [direct_npz, direct_vtu, restart_binary, elemlist_path, bnd_path]
    before = _file_manifest(critical_paths)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"{AUDIT_PREFIX}_{stamp}"
    qc_dir = run_root / "qc"
    qc_dir.mkdir(parents=True, exist_ok=False)
    manifest_path = qc_dir / "port_flux_audit_manifest.json"
    summary: dict[str, Any] = {
        "status": AUDIT_UNRESOLVED,
        "audit_revision": AUDIT_REVISION,
        "run_root": str(run_root),
        "branch": branch,
        "direct_decode_execution_base_commit": DIRECT_DECODE_EXECUTION_BASE_COMMIT,
        "current_evidence_commit": head,
        "external_apes_executable_calls": 0,
        "seeder_run_count": 0,
        "musubi_run_count": 0,
        "harvester_run_count": 0,
        "grid_convergence": "NOT_RUN",
        "production_solver_code_modified": False,
        "production_port_method_modified": False,
        "direct_restart_decoder_modified": False,
        "bc_modified": False,
        "surface_modified": False,
        "proteus_runtime_executed": False,
        "figures": [],
        "stable_window_policy_fixed_before_data_analysis": {
            "minimum_consecutive_depths": MINIMUM_STABLE_DEPTH_COUNT,
            "maximum_depth_cells": MAXIMUM_CUTSET_DEPTH,
            "maximum_inlet_mass_flow_relative_error": MAXIMUM_INLET_MASS_RELATIVE_ERROR,
            "maximum_signed_mass_balance_error": MAXIMUM_SIGNED_MASS_BALANCE_ERROR,
            "maximum_inlet_mass_flow_relative_spread": MAXIMUM_INLET_MASS_SPREAD,
            "maximum_each_outlet_mass_flow_relative_spread": MAXIMUM_OUTLET_MASS_SPREAD,
        },
        "frozen_files_before": before,
        "started_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, summary)
    try:
        if (
            direct_manifest.get("macro_reconstruction_validation") != "PASS"
            or direct_manifest.get("field_identity", {}).get("status") != "PASS"
        ):
            raise FlowError(AUDIT_UNRESOLVED, "Frozen direct field is not validated")

        source_contract, side_names, side_offsets, mfr_contract = build_source_contract(
            distribution=config.apes.wsl_distribution,
            current_evidence_commit=head,
        )
        write_json(qc_dir / "port_flux_source_contract.json", source_contract)
        write_json(qc_dir / "mfr_eq_boundary_contract.json", mfr_contract)

        mesh_header_text = (mesh_dir / "header.lua").read_text(encoding="utf-8")
        property_header = parse_boundary_property_header(mesh_header_text)
        bnd_header = parse_bnd_header((mesh_dir / "bnd.lua").read_text(encoding="utf-8"))
        if (
            property_header.bit_position != 3
            or property_header.element_count != EXPECTED_BOUNDARY_ELEMENT_COUNT
            or bnd_header.side_count != EXPECTED_SIDE_COUNT
            or bnd_header.boundary_type_count != 5
        ):
            raise FlowError(TOPOLOGY_FAILED, "Frozen boundary header contract changed")
        bnd_contract = boundary_binary_contract(
            bnd_path,
            element_count=property_header.element_count,
            side_count=bnd_header.side_count,
        )
        if bnd_contract["status"] != "PASS" or bnd_contract["expected_bytes"] != EXPECTED_BND_BYTES:
            raise FlowError(TOPOLOGY_FAILED, f"Frozen bnd.lsb contract failed: {bnd_contract}")

        tree_ids, property_bits, elemlist_contract = read_treelm_elemlist(
            elemlist_path,
            n_elems=EXPECTED_CELL_COUNT,
        )
        property_indices = extract_boundary_property_indices(
            property_bits,
            property_header.bit_position,
        )
        if len(property_indices) != property_header.element_count:
            raise FlowError(
                TOPOLOGY_FAILED,
                f"Boundary-property bit count is {len(property_indices)}, expected {property_header.element_count}",
            )
        boundary_ids = read_boundary_ids(
            bnd_path,
            element_count=property_header.element_count,
            side_count=bnd_header.side_count,
        )

        with np.load(direct_npz) as data:
            npz_tree_ids = np.asarray(data["tree_id"], dtype=np.int64)
            cell_ijk = np.asarray(data["cell_ijk"], dtype=np.int64)
            density_lattice = np.asarray(data["density_lattice"], dtype=np.float64)
            velocity_m_s = np.asarray(data["velocity_phy"], dtype=np.float64)
        if (
            len(npz_tree_ids) != EXPECTED_CELL_COUNT
            or not np.array_equal(npz_tree_ids, tree_ids)
            or cell_ijk.shape != (EXPECTED_CELL_COUNT, 3)
            or velocity_m_s.shape != (EXPECTED_CELL_COUNT, 3)
            or density_lattice.shape != (EXPECTED_CELL_COUNT,)
        ):
            raise FlowError(TOPOLOGY_FAILED, "Direct field and frozen TreElm element order differ")

        inputs = load_flow_inputs(config.paths.source_surface_run)
        partition = load_frozen_surface_partition(
            inputs,
            seeder_run / "geometry" / "geometry_solver_m",
        )
        boundary_conditions = load_boundary_conditions(inputs.boundary_conditions)
        if (
            not math.isclose(boundary_conditions.inlet_flow_m3_s, TARGET_Q_M3_S, rel_tol=0.0, abs_tol=0.0)
            or not math.isclose(
                boundary_conditions.density_kg_m3 * boundary_conditions.inlet_flow_m3_s,
                TARGET_MASS_FLOW_KG_S,
                rel_tol=1.0e-15,
                abs_tol=0.0,
            )
        ):
            raise FlowError(AUDIT_UNRESOLVED, "Frozen inlet target contract changed")

        topology_qc, topology = build_port_topology(
            boundary_ids=boundary_ids,
            property_element_indices=property_indices,
            boundary_header=bnd_header,
            side_names=side_names,
            side_offsets=side_offsets,
            partition=partition,
            dx_m=EXPECTED_DX_M,
        )
        topology_qc.update(
            {
                "boundary_property": {
                    "label": property_header.label,
                    "bit_position": property_header.bit_position,
                    "element_count": property_header.element_count,
                    "python_btest_count": int(len(property_indices)),
                    "mapping_status": "PASS",
                },
                "bnd_binary_contract": bnd_contract,
                "bnd_header": {
                    "side_count": bnd_header.side_count,
                    "boundary_type_count": bnd_header.boundary_type_count,
                    "labels": list(bnd_header.labels),
                },
                "elemlist_contract": elemlist_contract,
                "tree_element_order_matches_direct_field": True,
            }
        )
        write_json(qc_dir / "lattice_boundary_topology_qc.json", topology_qc)
        if topology_qc["status"] != "PASS":
            raise FlowError(TOPOLOGY_FAILED, "One or more lattice port topology gates failed")

        legacy_qc = read_json(direct_run / "qc" / "port_flux_qc.json")
        legacy = legacy_flux_semantics(legacy_qc)
        write_json(qc_dir / "legacy_geometric_flux_review.json", legacy)

        boundary_cell = boundary_cell_fluxes(
            topology,
            velocity_m_s,
            density_lattice,
            dx_m=EXPECTED_DX_M,
            rho0_kg_m3=REFERENCE_DENSITY_KG_M3,
        )
        write_json(qc_dir / "boundary_cell_flux_diagnostic.json", boundary_cell)

        graph = build_face_neighbor_graph(cell_ijk)
        seeds = {label: topology[label].seed_cell_indices for label in PORT_LABELS}
        sweep_rows, first_overlap = lattice_internal_cutset_sweep(
            graph=graph,
            cell_ijk=cell_ijk,
            velocity_m_s=velocity_m_s,
            density_lattice=density_lattice,
            seeds_by_port=seeds,
            dx_m=EXPECTED_DX_M,
            rho0_kg_m3=REFERENCE_DENSITY_KG_M3,
            maximum_depth=MAXIMUM_CUTSET_DEPTH,
        )
        depth_summaries = summarize_cutset_depths(
            sweep_rows,
            q_target_m3_s=TARGET_Q_M3_S,
            mass_target_kg_s=TARGET_MASS_FLOW_KG_S,
        )
        csv_path = qc_dir / "port_cutset_depth_sweep.csv"
        _write_depth_csv(csv_path, sweep_rows, depth_summaries)
        stable = find_stable_window(depth_summaries)
        cutset_qc = {
            "status": "PASS" if stable.get("found") else "UNRESOLVED",
            "estimator": "LATTICE_INTERNAL_CUTSET_ESTIMATOR",
            "depth_definition": (
                "depth 1 contains boundary seed cells; each additional depth adds one 6-neighbor layer"
            ),
            "depths_evaluated": list(range(1, MAXIMUM_CUTSET_DEPTH + 1)),
            "maximum_depth_cells": MAXIMUM_CUTSET_DEPTH,
            "maximum_depth_um": MAXIMUM_CUTSET_DEPTH * EXPECTED_DX_M * 1.0e6,
            "first_cross_port_overlap_depth": first_overlap,
            "depth_summaries": depth_summaries,
            "stable_window": stable,
            "csv": str(csv_path),
        }
        write_json(qc_dir / "lattice_internal_cutset_qc.json", cutset_qc)

        classification = classify_audit(
            topology_pass=True,
            stable_window=stable,
            legacy_outlet_02_m3_s=float(legacy_qc["q_out_m3_s"]["outlet_02"]),
        )
        after = _file_manifest(critical_paths)
        frozen_unchanged = before == after
        sha_qc = {
            "status": "PASS" if frozen_unchanged else "FAIL",
            "frozen_files_unchanged": frozen_unchanged,
            "before": before,
            "after": after,
        }
        write_json(qc_dir / "frozen_read_only_sha_qc.json", sha_qc)
        if not frozen_unchanged:
            classification = {
                "status": AUDIT_UNRESOLVED,
                "next": NEXT_BY_STATUS[AUDIT_UNRESOLVED],
                "outlet_02_backflow_confirmed": "UNRESOLVED",
                "legacy_geometric_estimator_classified_as_artifact": "UNRESOLVED",
            }

        summary.update(
            {
                **classification,
                "source_contract": str(qc_dir / "port_flux_source_contract.json"),
                "mfr_eq_contract": str(qc_dir / "mfr_eq_boundary_contract.json"),
                "legacy_flux": legacy,
                "boundary_topology": topology_qc,
                "boundary_cell_estimator": boundary_cell,
                "cutset_audit": cutset_qc,
                "frozen_files_after": after,
                "frozen_field_modified": False,
                "frozen_restart_modified": False,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        after = _file_manifest(critical_paths)
        status = error.status if isinstance(error, FlowError) else AUDIT_UNRESOLVED
        if status not in NEXT_BY_STATUS:
            status = AUDIT_UNRESOLVED
        summary.update(
            {
                "status": status,
                "next": NEXT_BY_STATUS[status],
                "failure": str(error),
                "frozen_files_after": after,
                "frozen_field_modified": before[str(direct_npz.resolve())] != after[str(direct_npz.resolve())],
                "frozen_restart_modified": before[str(restart_binary.resolve())]
                != after[str(restart_binary.resolve())],
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
