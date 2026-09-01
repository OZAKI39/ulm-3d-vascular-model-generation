"""Decode one frozen Musubi D3Q19 restart without running APES executables.

The decoder is deliberately tied to the locally pinned Musubi/TreElm revisions.
It first proves the binary, stencil, macroscopic-quantity, and mesh-order
contracts from those sources.  Only after that proof passes does it read the
restart PDFs and reconstruct the Cartesian hexahedral field.
"""

from __future__ import annotations

import math
import re
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np
from .io import FlowError, sha256_file
from .validated_contract import (
    BASE_DX_M,
    LATTICE_CS2,
    MAXIMUM_LATTICE_SPEED,
    RHO0_KG_M3,
    ValidatedTau1Contract,
    pressure_reference_pa,
)


SOURCE_COMMIT = "459a3b4c15f2d333eccbd2ce928391bd50573c93"
SOURCE_BRANCH = "codex/cfd-flow-musubi-recovery-20260828"
SOURCE_RUN_NAME = "musubi_project_steady_confirmation_anchor003274_20260828_225334"
FROZEN_SEEDER_NAME = "musubi_recovery_anchor003274_20260828_162530"
HISTORICAL_ASCII_EXPORT_NAME = "musubi_steady_field_export_anchor003274_20260828_233011"
HISTORICAL_VTK_DIAGNOSTIC_NAME = "musubi_diagnostic_anchor003274_20260828_182034"
EXPORT_PREFIX = "musubi_direct_restart_field_anchor003274"
EXPORT_REVISION = "FROZEN_MUSUBI_D3Q19_DIRECT_RESTART_DECODE_V1"

SOURCE_STATE = "FROZEN_PROJECT_STEADY_RESTART"
FROZEN_ITERATION = 198_064
EXPECTED_CELL_COUNT = 182_320
EXPECTED_PDF_COMPONENTS = 19
EXPECTED_NDOFS = 1
EXPECTED_RESTART_BYTES = 27_712_640
EXPECTED_ELEMLIST_BYTES = 2_917_120
EXPECTED_LEVEL = 9
EXPECTED_DX_M = BASE_DX_M
EXPECTED_DT_S = ValidatedTau1Contract().dt_s
REFERENCE_DENSITY_KG_M3 = RHO0_KG_M3
PRESSURE_REFERENCE_PA = ValidatedTau1Contract().pressure_reference_pa
PRESSURE_ABSOLUTE_TOLERANCE_PA = 1.0e-6
SPEED_ABSOLUTE_TOLERANCE_M_S = 1.0e-12
REFERENCE_RELATIVE_TOLERANCE = 1.0e-8
TOTAL_DENSITY_PRINT_TOLERANCE = 5.0e-8
MAXIMUM_LATTICE_MACH = MAXIMUM_LATTICE_SPEED
CS2 = LATTICE_CS2

SUCCESS_STATUS = "CFD_FLOW_DIRECT_RESTART_FIELD_EXPORT_PASS_PENDING_GRID_CONVERGENCE"
SUCCESS_NEXT = "RUN CFD GRID-SPACING CONVERGENCE STUDY"
PORT_REVIEW_STATUS = "CFD_FLOW_DIRECT_FIELD_PASS_PORT_QC_REVIEW_NEEDED"
PORT_REVIEW_NEXT = "REVIEW PORT FLUX INTEGRATION ONLY"
CONTRACT_UNPROVEN = "CFD_FLOW_DIRECT_RESTART_DECODE_CONTRACT_UNPROVEN"
BINARY_SIZE_INVALID = "CFD_FLOW_DIRECT_RESTART_BINARY_SIZE_INVALID"
MACRO_VALIDATION_FAILED = "CFD_FLOW_DIRECT_RESTART_MACRO_VALIDATION_FAILED"
MESH_LAYOUT_UNPROVEN = "CFD_FLOW_DIRECT_TREELM_MESH_LAYOUT_UNPROVEN"
DIRECT_DECODE_INVALID = "CFD_FLOW_DIRECT_RESTART_DECODE_INVALID"

MUSUBI_SOURCE_COMMIT = "4e8b277b66226277171ef93bf054d21270812793"
MUSUBI_SCHEME_COMMIT = "81f8c4f13772f6d4af31f335e1e3f99b02726e25"
TREELM_SOURCE_COMMIT = "9899d1376992c4fafc8a343d2b4ccef81de670d1"
MUSUBI_EXECUTABLE_SHA256 = "a005b4f00bd45df0339adc22460f251c3f300f967ff746c1cd43fa5ad7c07e88"
MUSUBI_SOLVER_TAG = "Musubi_v2.0.0-4-g4e8b27"

D3Q19_DIRECTIONS = np.asarray(
    (
        (-1, 0, 0),
        (0, -1, 0),
        (0, 0, -1),
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (0, -1, -1),
        (0, -1, 1),
        (0, 1, -1),
        (0, 1, 1),
        (-1, 0, -1),
        (1, 0, -1),
        (-1, 0, 1),
        (1, 0, 1),
        (-1, -1, 0),
        (-1, 1, 0),
        (1, -1, 0),
        (1, 1, 0),
        (0, 0, 0),
    ),
    dtype=np.int8,
)
D3Q19_WEIGHTS = np.asarray((*(1.0 / 18.0 for _ in range(6)), *(1.0 / 36.0 for _ in range(12)), 1.0 / 3.0))


@dataclass(frozen=True, slots=True)
class DirectDecodeLayout:
    root: Path
    input: Path
    flow: Path
    qc: Path
    proteus: Path
    figures: Path


@dataclass(frozen=True, slots=True)
class RestartHeader:
    binary_path: Path
    binary_name_wsl: str
    solver_config: Path
    solver_config_wsl: str
    iteration: int
    n_elems: int
    n_dofs: int
    variable_name: str
    n_components: int
    n_scalars: int
    solver_tag: str


@dataclass(frozen=True, slots=True)
class DirectMacroscopicField:
    density_lattice: np.ndarray
    velocity_lattice: np.ndarray
    velocity_phy: np.ndarray
    pressure_phy: np.ndarray


def _create_layout(output_root: Path) -> DirectDecodeLayout:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    root = Path(output_root) / f"{EXPORT_PREFIX}_{stamp}"
    if root.exists():
        raise FlowError(DIRECT_DECODE_INVALID, f"Output already exists: {root}")
    directories = {name: root / name for name in ("input", "flow", "qc", "proteus", "figures")}
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=False)
    return DirectDecodeLayout(root=root, **directories)


def _git_value(root: Path, *arguments: str) -> str:
    process = subprocess.run(
        ["git", *arguments],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return process.stdout.strip()


def _wsl_git_value(distribution: str, repository: str, *arguments: str) -> str:
    process = subprocess.run(
        ["wsl.exe", "-d", distribution, "--", "git", "-C", repository, *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
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
    raise FlowError(CONTRACT_UNPROVEN, "Pinned Musubi source directory is unavailable")


def _wsl_to_windows_path(value: str, *, relative_root: Path | None = None) -> Path:
    match = re.fullmatch(r"/mnt/([A-Za-z])/(.+)", value)
    if match:
        return Path(f"{match.group(1).upper()}:/{match.group(2)}")
    if value.startswith("/"):
        return Path(r"\\wsl.localhost\Ubuntu" + value.replace("/", "\\"))
    if relative_root is not None:
        return Path(relative_root).joinpath(*PurePosixPath(value).parts).resolve()
    raise FlowError(CONTRACT_UNPROVEN, f"Unsupported frozen WSL path: {value}")


def _extract_string(text: str, key: str) -> str:
    match = re.search(rf"\b{re.escape(key)}\s*=\s*['\"]([^'\"]+)['\"]", text)
    if not match:
        raise FlowError(CONTRACT_UNPROVEN, f"Missing Lua string: {key}")
    return match.group(1)


def _extract_integer(text: str, key: str) -> int:
    match = re.search(rf"\b{re.escape(key)}\s*=\s*(\d+)", text)
    if not match:
        raise FlowError(CONTRACT_UNPROVEN, f"Missing Lua integer: {key}")
    return int(match.group(1))


def parse_restart_header(path: Path) -> RestartHeader:
    """Parse only source-proven fields from the frozen restart header."""

    text = Path(path).read_text(encoding="utf-8", errors="strict")
    binary_match = re.search(r"binary_name\s*=\s*\{\s*['\"]([^'\"]+)['\"]", text, re.DOTALL)
    variable_match = re.search(
        r"variable\s*=\s*\{.*?name\s*=\s*['\"]([^'\"]+)['\"].*?ncomponents\s*=\s*(\d+)",
        text,
        re.DOTALL,
    )
    if not binary_match or not variable_match:
        raise FlowError(CONTRACT_UNPROVEN, "Restart binary or variable metadata is missing")
    binary_wsl = binary_match.group(1)
    solver_config_wsl = _extract_string(text, "solver_configFile")
    # Musubi normally writes paths relative to the configuration directory.
    # A header is one directory below that root (for example restart/header.lua).
    relative_root = Path(path).parent.parent
    return RestartHeader(
        binary_path=_wsl_to_windows_path(binary_wsl, relative_root=relative_root),
        binary_name_wsl=binary_wsl,
        solver_config=_wsl_to_windows_path(
            solver_config_wsl, relative_root=relative_root
        ),
        solver_config_wsl=solver_config_wsl,
        iteration=_extract_integer(text, "iter"),
        n_elems=_extract_integer(text, "nElems"),
        n_dofs=_extract_integer(text, "nDofs"),
        variable_name=variable_match.group(1),
        n_components=int(variable_match.group(2)),
        n_scalars=_extract_integer(text, "nScalars"),
        solver_tag=_extract_string(text, "solver"),
    )


def parse_d3q19_layout(source_text: str) -> np.ndarray:
    """Parse the 19 ordered direction rows from TreElm's pinned source block."""

    uncommented = "\n".join(line.split("!", 1)[0] for line in source_text.splitlines())
    normalized = uncommented.replace("&", " ")
    block_match = re.search(
        r"integer\s*,\s*parameter\s*::\s*d3q19_cxDir\s*\(\s*3\s*,\s*19\s*\).*?"
        r"reshape\s*\(\s*\[(.*?)\]\s*,\s*\[\s*3\s*,\s*19\s*\]",
        normalized,
        re.DOTALL | re.IGNORECASE,
    )
    if not block_match:
        raise FlowError(CONTRACT_UNPROVEN, "Could not locate pinned d3q19_cxDir source block")
    values = [int(value) for value in re.findall(r"(?<![A-Za-z0-9_])[-+]?\d+", block_match.group(1))]
    if len(values) != 3 * EXPECTED_PDF_COMPONENTS:
        raise FlowError(CONTRACT_UNPROVEN, f"Expected 57 D3Q19 values, found {len(values)}")
    return np.asarray(values, dtype=np.int8).reshape(EXPECTED_PDF_COMPONENTS, 3)


def _line_number(text: str, token: str) -> int:
    index = text.find(token)
    if index < 0:
        raise FlowError(CONTRACT_UNPROVEN, f"Pinned source token is absent: {token}")
    return text.count("\n", 0, index) + 1


def _source_file_evidence(
    path: Path,
    revision: str,
    statements: dict[str, tuple[str, ...]],
) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="strict")
    evidence: dict[str, Any] = {}
    for label, tokens in statements.items():
        lines = [_line_number(text, token) for token in tokens]
        evidence[label] = {"status": "PASS", "line_numbers": lines, "tokens": list(tokens)}
    return {
        "path": str(path),
        "source_revision": revision,
        "sha256": sha256_file(path),
        "evidence": evidence,
    }


def build_direct_decode_source_contract(
    *,
    distribution: str,
    solver_config: Path,
    restart_header: RestartHeader,
) -> tuple[dict[str, Any], np.ndarray]:
    """Prove every binary/macro/mesh assumption from the pinned source tree."""

    root = _wsl_source_root(distribution)
    mus_root_posix = "/home/lzy/apes-pinned/musubi_official"
    repositories = {
        "musubi": (mus_root_posix, MUSUBI_SOURCE_COMMIT),
        "musubi_scheme": (f"{mus_root_posix}/mus", MUSUBI_SCHEME_COMMIT),
        "treelm": (f"{mus_root_posix}/tem", TREELM_SOURCE_COMMIT),
    }
    revisions: dict[str, str] = {}
    tracked_status: dict[str, str] = {}
    for label, (repository, expected) in repositories.items():
        actual = _wsl_git_value(distribution, repository, "rev-parse", "HEAD")
        status = _wsl_git_value(distribution, repository, "status", "--short", "--untracked-files=no")
        revisions[label] = actual
        tracked_status[label] = status
        if actual != expected or status:
            raise FlowError(
                CONTRACT_UNPROVEN,
                f"Pinned {label} revision/cleanliness mismatch: {actual}, status={status!r}",
            )

    executable = Path(rf"\\wsl.localhost\{distribution}\home\lzy\.local\bin\musubi")
    if not executable.is_file():
        executable = Path(rf"\\wsl$\{distribution}\home\lzy\.local\bin\musubi")
    executable_sha = sha256_file(executable)
    if executable_sha != MUSUBI_EXECUTABLE_SHA256:
        raise FlowError(CONTRACT_UNPROVEN, "Pinned Musubi executable SHA-256 changed")

    tem = root / "tem" / "source"
    mus = root / "mus" / "source"
    source_paths = {
        "restart_serialization": tem / "tem_restart_module.f90",
        "real_kind_and_endian": tem / "env_module.f90",
        "treelm_mesh_io": tem / "treelmesh_module.f90",
        "treelm_topology": tem / "tem_topology_module.f90",
        "d3q19_layout": tem / "tem_stencil_module.fpp",
        "musubi_restart_order": mus / "mus_restart_module.f90",
        "musubi_pdf_serialization": mus / "mus_buffer_module.fpp",
        "density_velocity_derivation": mus / "derived" / "mus_derQuan_module.fpp",
        "d3q19_velocity_derivation": mus / "scheme" / "mus_scheme_derived_quantities_type_module.f90",
        "pressure_derivation": mus / "derived" / "mus_derQuan_module.fpp",
        "physical_derivation": mus / "derived" / "mus_derQuanPhysics_module.fpp",
        "physical_conversion": mus / "mus_physics_module.f90",
    }
    file_evidence = {
        "restart_serialization": _source_file_evidence(
            source_paths["restart_serialization"],
            TREELM_SOURCE_COMMIT,
            {
                "contiguous_nscalars_ndofs_rk": (
                    "call MPI_Type_contiguous( me%varMap%nScalars*me%write_file%nDofs",
                    "rk_mpi",
                ),
                "eight_byte_element_size": ("* 8_MPI_OFFSET_KIND",),
                "elementwise_variable_system_order": (
                    "Within each variable system the data is",
                    "organized elementwise.",
                ),
                "native_mpi_io": ('me%read_file%ftype, "native"',),
            },
        ),
        "real_kind_and_endian": _source_file_evidence(
            source_paths["real_kind_and_endian"],
            TREELM_SOURCE_COMMIT,
            {
                "rk_is_double": (
                    "rk_prec = double_prec",
                    "rk_mpi = mpi_double_precision",
                ),
                "long_integer_is_eight_bytes": ("long_k_mpi = MPI_INTEGER8",),
                "lsb_means_little_endian": (
                    "isLittleEndian = Sys_is_Little_Endian()",
                    "suffix = '.lsb'",
                ),
            },
        ),
        "treelm_mesh_io": _source_file_evidence(
            source_paths["treelm_mesh_io"],
            TREELM_SOURCE_COMMIT,
            {
                "two_interleaved_integer8_values": (
                    "call MPI_Type_contiguous( 2, long_k_mpi, etype, iError )",
                    "me%treeID(iElem)           = buffer( (iElem-1)*2 + 1 )",
                    "me%ElemPropertyBits(iElem) = buffer( (iElem-1)*2 + 2 )",
                ),
                "dump_order_matches_load_order": (
                    "buffer((iElem-1)*2+1) = me%treeID(iElem)",
                    "buffer((iElem-1)*2+2) = me%ElemPropertyBits(iElem)",
                ),
                "elemlist_endian_suffix": (
                    "Dump treeID and propertyBits to elemlist.lsb",
                    "tem_create_EndianSuffix()",
                ),
            },
        ),
        "treelm_topology": _source_file_evidence(
            source_paths["treelm_topology"],
            TREELM_SOURCE_COMMIT,
            {
                "level_formula": ("tElem = ishft( (tElem-1), -3 )",),
                "first_id_formula": ("res = ( 8_long_k**level - 1_long_k ) / 7_long_k",),
                "morton_decode": (
                    "tElem = TreeID - FirstIdAtLevel(coord(4))",
                    "coord(i) = coord(i) + bitlevel * int(mod(tElem / fak(i), 2_long_k))",
                ),
            },
        ),
        "musubi_restart_order": _source_file_evidence(
            source_paths["musubi_restart_order"],
            MUSUBI_SCHEME_COMMIT,
            {
                "restart_uses_treeid_chunks": (
                    "tree%treeID(elemOff+1:elemOff+restart%nChunkElems)",
                    "call tem_restart_writeData( restart, buffer )",
                ),
                "read_uses_same_treeid_chunks": (
                    "call tem_restart_readData( restart, buffer )",
                    "call mus_pdf_unserialize",
                ),
            },
        ),
        "musubi_pdf_serialization": _source_file_evidence(
            source_paths["musubi_pdf_serialization"],
            MUSUBI_SCHEME_COMMIT,
            {
                "original_treeid_order": (
                    "into chunks for writing it in original treeIDlist order to disk",
                    "do iElem = 1, nElems",
                    "do iComp = 1, nComp",
                    "buffer(iIndex) = scheme%state",
                ),
                "inverse_mapping_is_identical": (
                    "do iElem =  1, nElems",
                    "scheme%state(iLevel)%val",
                    ") = buffer(iIndex)",
                ),
            },
        ),
        "density_velocity_derivation": _source_file_evidence(
            source_paths["density_velocity_derivation"],
            MUSUBI_SCHEME_COMMIT,
            {
                "density_from_ordered_state": (
                    "subroutine deriveRho_FromState",
                    "state( iDir+(iElem-1)*varSys%nScalars )",
                ),
                "velocity_from_ordered_state": (
                    "subroutine deriveVel_FromState",
                    "rho = sum( pdf )",
                    "vel = layout%quantities%vel_from_pdf_ptr(pdf = pdf, dens = rho)",
                ),
            },
        ),
        "d3q19_velocity_derivation": _source_file_evidence(
            source_paths["d3q19_velocity_derivation"],
            MUSUBI_SCHEME_COMMIT,
            {
                "fluid_d3q19_function_selection": (
                    "getQuantities%vel_from_pdf_ptr => get_vel_from_pdf_d3q19",
                ),
                "exact_d3q19_velocity": (
                    "pure function get_vel_from_pdf_d3q19",
                    "vel(1) = pdf(4) - pdf(1)",
                    "vel = vel / dens",
                ),
            },
        ),
        "pressure_and_physical_conversion": {
            "pressure_derivation": _source_file_evidence(
                source_paths["pressure_derivation"],
                MUSUBI_SCHEME_COMMIT,
                {"pressure_is_density_cs2": ("res(1:nElems) = res(1:nElems) * cs2",)},
            ),
            "physical_derivation": _source_file_evidence(
                source_paths["physical_derivation"],
                MUSUBI_SCHEME_COMMIT,
                {
                    "pressure_phy_uses_press_factor": (
                        "recursive subroutine derivePressurePhy",
                        "%fac( level( iElem ) )%press",
                    ),
                    "velocity_phy_uses_velocity_factor": (
                        "recursive subroutine deriveVelocityPhy",
                        "%fac( level( iElem ) )%vel",
                    ),
                },
            ),
            "physical_conversion": _source_file_evidence(
                source_paths["physical_conversion"],
                MUSUBI_SCHEME_COMMIT,
                {
                    "velocity_factor": (
                        "me%fac( iLevel )%vel = me%dxLvl( iLevel )/me%dtLvl( iLevel )",
                    ),
                    "pressure_factor": (
                        "me%fac( iLevel )%press = me%rho0 * me%dxLvl( iLevel )**2",
                        "/ me%dtLvl( iLevel )**2",
                    ),
                },
            ),
        },
    }

    stencil_text = source_paths["d3q19_layout"].read_text(encoding="utf-8", errors="strict")
    source_directions = parse_d3q19_layout(stencil_text)
    if not np.array_equal(source_directions, D3Q19_DIRECTIONS):
        raise FlowError(CONTRACT_UNPROVEN, "Pinned D3Q19 ordering differs from the decoder constant")
    solver_text = Path(solver_config).read_text(encoding="utf-8", errors="strict")
    force_assignments = re.findall(r"(?mi)^\s*(?:source|body_force|force)\s*=", solver_text)
    if force_assignments:
        raise FlowError(CONTRACT_UNPROVEN, "The frozen solver config contains a force/source assignment")
    solver_commit_match = re.search(r"-g([0-9a-f]+)$", restart_header.solver_tag)
    solver_commit_prefix = solver_commit_match.group(1) if solver_commit_match else ""
    if (
        restart_header.solver_tag != MUSUBI_SOLVER_TAG
        or not solver_commit_prefix
        or not MUSUBI_SOURCE_COMMIT.startswith(solver_commit_prefix)
    ):
        raise FlowError(CONTRACT_UNPROVEN, "Restart solver tag does not match pinned Musubi source")

    d3q19_path = source_paths["d3q19_layout"]
    contract = {
        "status": "PASS",
        "contract_revision": EXPORT_REVISION,
        "pinned_source_root": str(root),
        "git_pull_performed": False,
        "submodule_update_performed": False,
        "metadata_git_calls_only": True,
        "repositories": {
            label: {
                "path": repository,
                "commit": revisions[label],
                "expected_commit": expected,
                "tracked_worktree_clean": tracked_status[label] == "",
            }
            for label, (repository, expected) in repositories.items()
        },
        "musubi_executable_revision": {
            "path": str(executable),
            "solver_tag_from_restart": restart_header.solver_tag,
            "sha256": executable_sha,
            "solver_tag_commit_prefix": solver_commit_prefix,
            "source_commit_prefix_match": MUSUBI_SOURCE_COMMIT.startswith(solver_commit_prefix),
            "executed": False,
        },
        "restart_serialization_source_file": str(source_paths["restart_serialization"]),
        "d3q19_layout_source_file": str(d3q19_path),
        "density_derivation_source": str(source_paths["density_velocity_derivation"]),
        "velocity_derivation_source": str(source_paths["d3q19_velocity_derivation"]),
        "pressure_derivation_source": str(source_paths["pressure_derivation"]),
        "physical_conversion_source": str(source_paths["physical_conversion"]),
        "treelm_mesh_io_source": str(source_paths["treelm_mesh_io"]),
        "treelm_topology_source": str(source_paths["treelm_topology"]),
        "source_files": file_evidence,
        "restart_binary_contract": {
            "dtype": "little-endian IEEE-754 float64",
            "numpy_dtype": "<f8",
            "rk_bytes": 8,
            "values_per_element": EXPECTED_PDF_COMPONENTS * EXPECTED_NDOFS,
            "organization": "one variable system; elementwise; 19 contiguous rk values per element",
            "status": "PASS",
        },
        "treelm_elemlist_contract": {
            "dtype": "little-endian signed INTEGER8",
            "numpy_dtype": "<i8",
            "bytes_per_element": 16,
            "organization": "interleaved [treeID, ElemPropertyBits] per element",
            "status": "PASS",
        },
        "restart_mesh_index_correspondence": {
            "status": "PASS",
            "statement": "restart row k and elemlist record k both follow tree%treeID element order",
        },
        "body_force_contract": {
            "status": "PASS",
            "source_or_force_assignment_in_solver_config": False,
            "velocity_force_correction": "NONE_REQUIRED_FOR_THIS_FROZEN_FLUID_CONFIG",
        },
        "d3q19": {
            "status": "PASS",
            "qq": EXPECTED_PDF_COMPONENTS,
            "source_file": str(d3q19_path),
            "source_revision": TREELM_SOURCE_COMMIT,
            "source_sha256": sha256_file(d3q19_path),
            "rows": [
                {"pdf_index": index, "cx": int(row[0]), "cy": int(row[1]), "cz": int(row[2])}
                for index, row in enumerate(source_directions, start=1)
            ],
        },
    }
    return contract, source_directions


def restart_binary_size_contract(path: Path, *, n_elems: int, n_components: int, n_dofs: int) -> dict[str, Any]:
    expected = int(n_elems) * int(n_components) * int(n_dofs) * np.dtype("<f8").itemsize
    actual = Path(path).stat().st_size
    return {
        "status": "PASS" if actual == expected else "FAIL",
        "path": str(path),
        "actual_bytes": actual,
        "expected_bytes": expected,
        "n_elems": int(n_elems),
        "n_components": int(n_components),
        "n_dofs": int(n_dofs),
        "rk_bytes": 8,
        "dtype": "<f8",
    }


def read_restart_pdf(path: Path, *, n_elems: int, n_components: int) -> np.memmap:
    contract = restart_binary_size_contract(
        path,
        n_elems=n_elems,
        n_components=n_components,
        n_dofs=EXPECTED_NDOFS,
    )
    if contract["status"] != "PASS":
        raise FlowError(BINARY_SIZE_INVALID, f"Restart size contract failed: {contract}")
    return np.memmap(path, mode="r", dtype="<f8", shape=(n_elems, n_components), order="C")


def reconstruct_macroscopic_field(
    pdf: np.ndarray,
    *,
    directions: np.ndarray = D3Q19_DIRECTIONS,
    dx_m: float = EXPECTED_DX_M,
    dt_s: float = EXPECTED_DT_S,
    rho0_kg_m3: float = REFERENCE_DENSITY_KG_M3,
) -> DirectMacroscopicField:
    """Mirror pinned Musubi density, D3Q19 velocity, and pressure derivations."""

    values = np.asarray(pdf)
    stencil = np.asarray(directions, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] != EXPECTED_PDF_COMPONENTS:
        raise ValueError(f"Expected (n, 19) PDF array, found {values.shape}")
    if stencil.shape != (EXPECTED_PDF_COMPONENTS, 3):
        raise ValueError(f"Expected (19, 3) directions, found {stencil.shape}")
    density = np.sum(values, axis=1, dtype=np.float64)
    momentum = values @ stencil
    velocity_lattice = momentum / density[:, None]
    velocity_phy = velocity_lattice * (float(dx_m) / float(dt_s))
    pressure_phy = density * pressure_reference_pa(
        float(dx_m), float(dt_s), rho0_kg_m3=float(rho0_kg_m3), cs2=CS2
    )
    return DirectMacroscopicField(density, velocity_lattice, velocity_phy, pressure_phy)


def reference_reproduction(
    actual: float,
    reference: float,
    *,
    absolute_tolerance: float,
    relative_tolerance: float = REFERENCE_RELATIVE_TOLERANCE,
) -> dict[str, Any]:
    absolute_error = abs(float(actual) - float(reference))
    relative_error = absolute_error / abs(float(reference)) if reference else math.inf
    passed = absolute_error <= absolute_tolerance or relative_error <= relative_tolerance
    return {
        "status": "PASS" if passed else "FAIL",
        "actual": float(actual),
        "reference": float(reference),
        "absolute_error": float(absolute_error),
        "relative_error": float(relative_error),
        "absolute_tolerance": float(absolute_tolerance),
        "relative_tolerance": float(relative_tolerance),
        "policy": "absolute OR relative tolerance",
    }


def first_id_at_level(level: int) -> int:
    return (8**int(level) - 1) // 7


def tree_levels(tree_ids: np.ndarray) -> np.ndarray:
    ids = np.asarray(tree_ids, dtype=np.int64).reshape(-1)
    if np.any(ids < 0):
        raise ValueError("treeID must be non-negative")
    work = ids.copy()
    levels = np.zeros(len(ids), dtype=np.int16)
    while np.any(work):
        active = work != 0
        work[active] = (work[active] - 1) >> 3
        levels[active] += 1
    return levels


def tree_ids_to_ijk(tree_ids: np.ndarray, levels: np.ndarray | None = None) -> np.ndarray:
    """Mirror tem_CoordOfId bit de-interleaving for an array of treeIDs."""

    ids = np.asarray(tree_ids, dtype=np.int64).reshape(-1)
    level_values = tree_levels(ids) if levels is None else np.asarray(levels, dtype=np.int16)
    if level_values.shape != ids.shape:
        raise ValueError("levels must match tree_ids")
    result = np.empty((len(ids), 3), dtype=np.int32)
    for level in np.unique(level_values):
        mask = level_values == level
        morton = ids[mask] - first_id_at_level(int(level))
        coordinates = np.zeros((int(np.count_nonzero(mask)), 3), dtype=np.int32)
        for bit in range(int(level)):
            coordinates[:, 0] |= (((morton >> (3 * bit)) & 1) << bit).astype(np.int32)
            coordinates[:, 1] |= (((morton >> (3 * bit + 1)) & 1) << bit).astype(np.int32)
            coordinates[:, 2] |= (((morton >> (3 * bit + 2)) & 1) << bit).astype(np.int32)
        result[mask] = coordinates
    return result


def read_treelm_elemlist(path: Path, *, n_elems: int) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    expected_bytes = int(n_elems) * 2 * np.dtype("<i8").itemsize
    actual_bytes = Path(path).stat().st_size
    contract = {
        "status": "PASS" if actual_bytes == expected_bytes else "FAIL",
        "path": str(path),
        "actual_bytes": actual_bytes,
        "expected_bytes": expected_bytes,
        "bytes_per_element": 16,
        "dtype": "<i8",
        "layout": "interleaved [treeID, ElemPropertyBits]",
    }
    if contract["status"] != "PASS":
        raise FlowError(MESH_LAYOUT_UNPROVEN, f"elemlist size contract failed: {contract}")
    records = np.memmap(path, mode="r", dtype="<i8", shape=(n_elems, 2), order="C")
    return np.asarray(records[:, 0]).copy(), np.asarray(records[:, 1]).copy(), contract


def unique_cell_gate(indices: np.ndarray, *, expected_count: int = EXPECTED_CELL_COUNT) -> dict[str, Any]:
    cells = np.asarray(indices, dtype=np.int64)
    unique_count = int(len(np.unique(cells, axis=0)))
    duplicate_count = int(len(cells) - unique_count)
    return {
        "status": "PASS"
        if len(cells) == expected_count and unique_count == expected_count and duplicate_count == 0
        else "FAIL",
        "cell_count": int(len(cells)),
        "unique_cell_count": unique_count,
        "duplicate_cell_count": duplicate_count,
        "expected_cell_count": int(expected_count),
    }
