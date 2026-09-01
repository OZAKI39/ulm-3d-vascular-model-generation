"""Python-only audit of frozen Musubi port fluxes on the TreElm lattice."""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .io import FlowError, sha256_file
from .validated_contract import RHO0_KG_M3


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

EXPECTED_BOUNDARY_ELEMENT_COUNT = 75_358
EXPECTED_SIDE_COUNT = 26
EXPECTED_BND_BYTES = EXPECTED_BOUNDARY_ELEMENT_COUNT * EXPECTED_SIDE_COUNT * 8
REFERENCE_DENSITY_KG_M3 = RHO0_KG_M3

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
