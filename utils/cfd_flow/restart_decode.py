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
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv

from .apes import load_boundary_conditions
from .config import load_cfd_flow_config
from .geometry import SurfacePartition, load_frozen_surface_partition
from .io import FlowError, load_flow_inputs, read_json, sha256_file, write_json
from .qc import numerical_port_fluxes, percentile_summary, reynolds_diagnostics
from .steady_export import (
    LatticeMapping,
    _scatter_figure,
    build_proteus_metadata,
    parse_uniform_mesh_lattice,
    reconstruct_hexahedral_field,
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
EXPECTED_CELL_COUNT = 221_109
EXPECTED_PDF_COMPONENTS = 19
EXPECTED_NDOFS = 1
EXPECTED_RESTART_BYTES = 33_608_568
EXPECTED_ELEMLIST_BYTES = 3_537_744
EXPECTED_LEVEL = 9
EXPECTED_DX_M = 2.0e-7
EXPECTED_DT_S = 2.44140625e-8
REFERENCE_DENSITY_KG_M3 = 1056.0
PRESSURE_REFERENCE_PA = 23_622.32012800001
MUSUBI_MEAN_PRESSURE_PA = 23_682.89481126405
MUSUBI_MEAN_SPEED_M_S = 6.119872696805561e-5
MUSUBI_TOTAL_DENSITY = 221_675.9895069
PRESSURE_ABSOLUTE_TOLERANCE_PA = 1.0e-6
SPEED_ABSOLUTE_TOLERANCE_M_S = 1.0e-12
REFERENCE_RELATIVE_TOLERANCE = 1.0e-8
TOTAL_DENSITY_PRINT_TOLERANCE = 5.0e-8
MAXIMUM_LATTICE_MACH = 0.05
CS2 = 1.0 / 3.0

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


def _wsl_to_windows_path(value: str) -> Path:
    match = re.fullmatch(r"/mnt/([A-Za-z])/(.+)", value)
    if not match:
        raise FlowError(CONTRACT_UNPROVEN, f"Unsupported frozen WSL path: {value}")
    return Path(f"{match.group(1).upper()}:/{match.group(2)}")


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
    return RestartHeader(
        binary_path=_wsl_to_windows_path(binary_wsl),
        binary_name_wsl=binary_wsl,
        solver_config=_wsl_to_windows_path(solver_config_wsl),
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
    pressure_factor = float(rho0_kg_m3) * float(dx_m) ** 2 / float(dt_s) ** 2
    pressure_phy = density * CS2 * pressure_factor
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


def _file_manifest(paths: list[Path]) -> dict[str, Any]:
    return {
        str(path.resolve()): {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for path in paths
    }


def _historical_harvester_review(output_root: Path) -> dict[str, Any]:
    ascii_qc_path = output_root / HISTORICAL_ASCII_EXPORT_NAME / "qc" / "harvester_run_qc.json"
    vtk_qc_path = output_root / HISTORICAL_VTK_DIAGNOSTIC_NAME / "qc" / "harvester_failure_qc.json"
    ascii_qc = read_json(ascii_qc_path)
    vtk_qc = read_json(vtk_qc_path)
    companion_written = bool(ascii_qc.get("preserved_partial_files"))
    return {
        "status": "MUS_HARVESTING_POST_RESTART_PRE_OUTPUT_SIGSEGV",
        "interpretation": "FORMAT-INDEPENDENT POST-RESTART/PRE-DATA FAILURE; NOT CFD GEOMETRY/PHYSICS FAILURE",
        "vtk": {
            "source_qc": str(vtk_qc_path),
            "restart_read_completed": bool(vtk_qc["restart_read_completed"]),
            "complete_vtk_dataset_written": bool(vtk_qc["complete_vtk_dataset_written"]),
            "signal": vtk_qc["signal"],
        },
        "ascii_spatial": {
            "source_qc": str(ascii_qc_path),
            "restart_read_completed": bool(ascii_qc["restart_read_completed"]),
            "companion_metadata_written": companion_written,
            "ascii_data_file_started_writing": False,
            "ascii_data_file_generated": False,
            "signal": ascii_qc["signal"],
        },
        "historical_downstream_semantics": {
            "flow_direction": "NOT_EVALUATED_DUE_TO_NO_FIELD",
            "proteus_compatibility": "NOT_EVALUATED_DUE_TO_NO_FIELD",
        },
    }


def _port_diagnostics(
    fluxes: dict[str, float],
    partition: SurfacePartition,
    kinematic_viscosity_m2_s: float,
) -> dict[str, Any]:
    reynolds = reynolds_diagnostics(fluxes, partition, kinematic_viscosity_m2_s)
    return {
        patch.label: {
            "q_m3_s": float(fluxes[patch.label]),
            "area_m2": float(patch.area_um2 * 1.0e-12),
            "mean_velocity_m_s": float(
                abs(fluxes[patch.label]) / (patch.area_um2 * 1.0e-12)
            ),
            "equivalent_diameter_m": float(2.0 * patch.equivalent_radius_um * 1.0e-6),
            "reynolds": float(reynolds[patch.label]),
        }
        for patch in partition.patches
        if patch.label != "wall"
    }


def _existing_revision(output_root: Path) -> dict[str, Any] | None:
    for run_root in sorted(Path(output_root).glob(f"{EXPORT_PREFIX}_*"), reverse=True):
        manifest = run_root / "qc" / "direct_restart_decode_manifest.json"
        if not manifest.is_file():
            continue
        value = read_json(manifest)
        if value.get("export_revision") == EXPORT_REVISION:
            return value
    return None


def _source_read_only_audit(
    *,
    paths: list[Path],
    before: dict[str, Any],
    destination: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    after = _file_manifest(paths)
    unchanged = before == after
    result = {
        "status": "PASS" if unchanged else "FAIL",
        "frozen_source_files_unchanged": unchanged,
        "before": before,
        "after": after,
    }
    write_json(destination, result)
    return result, after


def run_direct_restart_decode(project_root: Path) -> dict[str, Any]:
    """Run one Python-only direct restart decode with zero APES executable calls."""

    root = Path(project_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = config.paths.output_root
    previous = _existing_revision(output_root)
    if previous is not None:
        return previous

    branch = _git_value(root, "branch", "--show-current")
    head = _git_value(root, "rev-parse", "HEAD")
    if branch != SOURCE_BRANCH or head != SOURCE_COMMIT:
        raise FlowError(
            CONTRACT_UNPROVEN,
            f"Expected {SOURCE_BRANCH}@{SOURCE_COMMIT}, found {branch}@{head}",
        )

    source_run = output_root / SOURCE_RUN_NAME
    frozen_seeder = output_root / FROZEN_SEEDER_NAME
    restart_header_path = source_run / "restart" / "roi003274_steady_lbm_lastHeader.lua"
    mesh_dir = frozen_seeder / "seeder" / "mesh"
    mesh_header_path = mesh_dir / "header.lua"
    elemlist_path = mesh_dir / "elemlist.lsb"
    header = parse_restart_header(restart_header_path)
    layout = _create_layout(output_root)
    manifest_path = layout.qc / "direct_restart_decode_manifest.json"
    summary: dict[str, Any] = {
        "status": "CFD_FLOW_DIRECT_RESTART_DECODE_PREFLIGHT",
        "export_revision": EXPORT_REVISION,
        "run_root": str(layout.root),
        "branch": branch,
        "source_commit": head,
        "source_state": SOURCE_STATE,
        "frozen_iteration": FROZEN_ITERATION,
        "external_apes_executable_calls": 0,
        "seeder_run_count": 0,
        "musubi_run_count": 0,
        "harvester_run_count": 0,
        "grid_sweep_performed": False,
        "physics_or_bc_modified": False,
        "source_steady_restart_modified": False,
        "microbubble_simulation_run": False,
        "proteus_runtime_executed": False,
        "started_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, summary)

    critical_paths: list[Path] = []
    critical_before: dict[str, Any] = {}
    try:
        steady = read_json(source_run / "qc" / "project_steady_confirmation.json")
        steady_restart = read_json(source_run / "qc" / "steady_restart_manifest.json")
        if (
            steady.get("status") != "CFD_FLOW_PROJECT_STEADY_0P1PCT_CONFIRMED"
            or not steady["official_steady_termination"]["official_steady_termination"]
            or int(steady["official_steady_termination"]["confirmation_iteration"])
            != FROZEN_ITERATION
            or steady_restart.get("status") != SOURCE_STATE
        ):
            raise FlowError(CONTRACT_UNPROVEN, "Frozen steady evidence does not match")
        header_checks = {
            "iteration": header.iteration == FROZEN_ITERATION,
            "n_elems": header.n_elems == EXPECTED_CELL_COUNT,
            "n_dofs": header.n_dofs == EXPECTED_NDOFS,
            "variable_name": header.variable_name == "pdf",
            "n_components": header.n_components == EXPECTED_PDF_COMPONENTS,
            "n_scalars": header.n_scalars == EXPECTED_PDF_COMPONENTS,
            "solver_tag": header.solver_tag == MUSUBI_SOLVER_TAG,
            "binary_name": header.binary_path.name == "roi003274_steady_lbm_4.836E-03.lsb",
            "binary_exists": header.binary_path.is_file(),
            "solver_config_exists": header.solver_config.is_file(),
        }
        header_contract = {
            "status": "PASS" if all(header_checks.values()) else "FAIL",
            "checks": header_checks,
            "iteration": header.iteration,
            "n_elems": header.n_elems,
            "n_dofs": header.n_dofs,
            "variable_name": header.variable_name,
            "n_components": header.n_components,
            "n_scalars": header.n_scalars,
            "solver_tag": header.solver_tag,
            "binary": str(header.binary_path),
            "solver_config": str(header.solver_config),
        }
        write_json(layout.qc / "restart_header_contract.json", header_contract)
        if header_contract["status"] != "PASS":
            raise FlowError(CONTRACT_UNPROVEN, f"Restart header contract failed: {header_checks}")

        source_contract, source_directions = build_direct_decode_source_contract(
            distribution=config.apes.wsl_distribution,
            solver_config=header.solver_config,
            restart_header=header,
        )
        write_json(layout.qc / "direct_decode_source_contract.json", source_contract)
        write_json(layout.qc / "d3q19_layout_contract.json", source_contract["d3q19"])

        lattice = parse_uniform_mesh_lattice(mesh_header_path)
        lattice_checks = {
            "source_status": lattice["status"] == "PASS",
            "n_elems": int(lattice["n_elems"]) == EXPECTED_CELL_COUNT,
            "minimum_level": int(lattice["minimum_level"]) == EXPECTED_LEVEL,
            "maximum_level": int(lattice["maximum_level"]) == EXPECTED_LEVEL,
            "dx": math.isclose(float(lattice["dx_m"]), EXPECTED_DX_M, rel_tol=0.0, abs_tol=1.0e-21),
        }
        if not all(lattice_checks.values()):
            raise FlowError(MESH_LAYOUT_UNPROVEN, f"Frozen mesh header mismatch: {lattice_checks}")
        write_json(
            layout.qc / "frozen_mesh_header_contract.json",
            {"status": "PASS", "checks": lattice_checks, **lattice},
        )

        binary_size = restart_binary_size_contract(
            header.binary_path,
            n_elems=header.n_elems,
            n_components=header.n_components,
            n_dofs=header.n_dofs,
        )
        write_json(layout.qc / "restart_binary_size_contract.json", binary_size)
        if binary_size["status"] != "PASS" or binary_size["expected_bytes"] != EXPECTED_RESTART_BYTES:
            raise FlowError(BINARY_SIZE_INVALID, f"Restart size contract failed: {binary_size}")
        if elemlist_path.stat().st_size != EXPECTED_ELEMLIST_BYTES:
            raise FlowError(MESH_LAYOUT_UNPROVEN, "Frozen elemlist byte size changed")

        critical_paths = [restart_header_path, header.binary_path, elemlist_path, mesh_header_path]
        critical_before = _file_manifest(critical_paths)
        summary.update(
            {
                "status": "CFD_FLOW_DIRECT_RESTART_SOURCE_CONTRACT_PASS",
                "steady_evidence": {
                    "status": steady["status"],
                    "official_steady_termination": True,
                    "confirmation_iteration": FROZEN_ITERATION,
                },
                "restart_header_contract": header_contract,
                "direct_decode_source_contract": str(
                    layout.qc / "direct_decode_source_contract.json"
                ),
                "d3q19_layout_contract": str(layout.qc / "d3q19_layout_contract.json"),
                "restart_binary_size_contract": binary_size,
                "critical_source_manifest_before": critical_before,
            }
        )
        write_json(manifest_path, summary)

        # The binary is intentionally not opened until every static source contract above passes.
        pdf = read_restart_pdf(
            header.binary_path,
            n_elems=header.n_elems,
            n_components=header.n_components,
        )
        nan_count = int(np.isnan(pdf).sum())
        inf_count = int(np.isinf(pdf).sum())
        field = reconstruct_macroscopic_field(
            pdf,
            directions=source_directions,
            dx_m=EXPECTED_DX_M,
            dt_s=EXPECTED_DT_S,
            rho0_kg_m3=REFERENCE_DENSITY_KG_M3,
        )
        speed = np.linalg.norm(field.velocity_phy, axis=1)
        mean_pressure = float(np.mean(field.pressure_phy, dtype=np.float64))
        mean_speed = float(np.mean(speed, dtype=np.float64))
        total_density = float(np.sum(field.density_lattice, dtype=np.float64))
        pressure_validation = reference_reproduction(
            mean_pressure,
            MUSUBI_MEAN_PRESSURE_PA,
            absolute_tolerance=PRESSURE_ABSOLUTE_TOLERANCE_PA,
        )
        speed_validation = reference_reproduction(
            mean_speed,
            MUSUBI_MEAN_SPEED_M_S,
            absolute_tolerance=SPEED_ABSOLUTE_TOLERANCE_M_S,
        )
        density_validation = reference_reproduction(
            total_density,
            MUSUBI_TOTAL_DENSITY,
            absolute_tolerance=TOTAL_DENSITY_PRINT_TOLERANCE,
        )
        finite_pass = nan_count == 0 and inf_count == 0
        macro_pass = (
            finite_pass
            and pressure_validation["status"] == "PASS"
            and speed_validation["status"] == "PASS"
            and density_validation["status"] == "PASS"
        )
        macro_qc = {
            "status": "PASS" if macro_pass else "FAIL",
            "binary_sha256": sha256_file(header.binary_path),
            "binary_bytes": header.binary_path.stat().st_size,
            "pdf_dtype": str(pdf.dtype),
            "pdf_shape": list(pdf.shape),
            "pdf_elementwise_mapping": "f[k, 0:19] is restart element k",
            "nan_count": nan_count,
            "inf_count": inf_count,
            "finite": finite_pass,
            "mean_pressure_validation": pressure_validation,
            "mean_speed_validation": speed_validation,
            "total_density_validation": density_validation,
            "physical_conversion": {
                "dx_m": EXPECTED_DX_M,
                "dt_s": EXPECTED_DT_S,
                "rho0_kg_m3": REFERENCE_DENSITY_KG_M3,
                "velocity_factor_m_s": EXPECTED_DX_M / EXPECTED_DT_S,
                "pressure_factor_pa": REFERENCE_DENSITY_KG_M3
                * EXPECTED_DX_M**2
                / EXPECTED_DT_S**2,
                "cs2": CS2,
            },
        }
        write_json(layout.qc / "direct_macro_validation.json", macro_qc)
        if not macro_pass:
            raise FlowError(MACRO_VALIDATION_FAILED, f"Direct macro validation failed: {macro_qc}")

        tree_ids, property_bits, elemlist_contract = read_treelm_elemlist(
            elemlist_path,
            n_elems=EXPECTED_CELL_COUNT,
        )
        levels = tree_levels(tree_ids)
        cell_ijk = tree_ids_to_ijk(tree_ids, levels)
        cell_gate = unique_cell_gate(cell_ijk)
        unique_tree_ids = int(len(np.unique(tree_ids)))
        index_min = np.min(cell_ijk, axis=0)
        index_max = np.max(cell_ijk, axis=0)
        level_pass = bool(np.all(levels == EXPECTED_LEVEL))
        bounds_pass = bool(np.all(cell_ijk >= 0) and np.all(cell_ijk < 2**EXPECTED_LEVEL))
        tree_pass = (
            len(tree_ids) == EXPECTED_CELL_COUNT
            and unique_tree_ids == EXPECTED_CELL_COUNT
            and level_pass
            and bounds_pass
            and cell_gate["status"] == "PASS"
        )
        origin = np.asarray(lattice["origin_m"], dtype=np.float64)
        centers_m = origin[None, :] + (cell_ijk.astype(np.float64) + 0.5) * EXPECTED_DX_M
        center_unique = int(len(np.unique(centers_m, axis=0)))
        coordinate_bounds_m = {
            "minimum": np.min(centers_m, axis=0).tolist(),
            "maximum": np.max(centers_m, axis=0).tolist(),
        }
        coordinate_pass = bool(
            np.all(np.isfinite(centers_m))
            and center_unique == EXPECTED_CELL_COUNT
            and np.all(centers_m >= origin[None, :])
            and np.all(centers_m <= origin[None, :] + float(lattice["side_m"]))
        )
        mesh_qc = {
            "status": "PASS" if tree_pass and coordinate_pass else "FAIL",
            "elemlist": elemlist_contract,
            "elemlist_sha256": sha256_file(elemlist_path),
            "tree_id_count": int(len(tree_ids)),
            "unique_tree_id_count": unique_tree_ids,
            "property_bits_count": int(len(property_bits)),
            "all_levels": {str(int(value)): int(count) for value, count in zip(*np.unique(levels, return_counts=True), strict=True)},
            "all_level_9": level_pass,
            "cell_index_minimum": index_min.tolist(),
            "cell_index_maximum": index_max.tolist(),
            "indices_inside_level_bounds": bounds_pass,
            "cell_identity": cell_gate,
            "unique_center_count": center_unique,
            "coordinate_bounds_m": coordinate_bounds_m,
            "coordinates_finite_and_inside_root": coordinate_pass,
            "dx_m": EXPECTED_DX_M,
            "restart_mesh_index_correspondence": source_contract[
                "restart_mesh_index_correspondence"
            ],
        }
        write_json(layout.qc / "treelm_mesh_decode_qc.json", mesh_qc)
        if mesh_qc["status"] != "PASS":
            raise FlowError(MESH_LAYOUT_UNPROVEN, f"TreElm mesh decode failed: {mesh_qc}")

        npz_path = layout.flow / "direct_cell_field.npz"
        np.savez_compressed(
            npz_path,
            tree_id=tree_ids,
            cell_ijk=cell_ijk,
            center_m=centers_m,
            density_lattice=field.density_lattice,
            pressure_phy=field.pressure_phy,
            velocity_phy=field.velocity_phy,
        )

        mapping = LatticeMapping(
            cell_indices=cell_ijk.astype(np.int64),
            maximum_alignment_error_m=0.0,
            duplicate_cell_count=0,
            unique_cell_count=EXPECTED_CELL_COUNT,
        )
        grid = reconstruct_hexahedral_field(
            mapping=mapping,
            origin_m=origin,
            dx_m=EXPECTED_DX_M,
            pressure_pa=field.pressure_phy,
            velocity_m_s=field.velocity_phy,
            pressure_reference_pa=PRESSURE_REFERENCE_PA,
        )
        grid.cell_data["tree_id"] = tree_ids
        grid.cell_data["density_lattice"] = field.density_lattice
        flow_vtu = layout.flow / "flow_field.vtu"
        grid.save(flow_vtu, binary=True)
        reloaded = pv.read(flow_vtu).cast_to_unstructured_grid()
        cell_types, cell_type_counts = np.unique(reloaded.celltypes, return_counts=True)
        velocity_reload = np.asarray(reloaded.cell_data["velocity_phy"], dtype=np.float64)
        pressure_reload = np.asarray(reloaded.cell_data["pressure_phy"], dtype=np.float64)
        gauge_reload = np.asarray(reloaded.cell_data["pressure_gauge_pa"], dtype=np.float64)
        speed_reload = np.linalg.norm(velocity_reload, axis=1)
        vtu_binary = b'format="binary"' in flow_vtu.read_bytes()[:8192] or b'format="appended"' in flow_vtu.read_bytes()[:8192]
        field_finite = bool(
            np.all(np.isfinite(reloaded.points))
            and np.all(np.isfinite(velocity_reload))
            and np.all(np.isfinite(pressure_reload))
            and np.all(np.isfinite(gauge_reload))
        )
        field_identity_pass = bool(
            vtu_binary
            and reloaded.n_cells == EXPECTED_CELL_COUNT
            and len(cell_types) == 1
            and int(cell_types[0]) == int(pv.CellType.HEXAHEDRON)
            and velocity_reload.shape == (EXPECTED_CELL_COUNT, 3)
            and pressure_reload.shape == (EXPECTED_CELL_COUNT,)
            and gauge_reload.shape == (EXPECTED_CELL_COUNT,)
            and field_finite
        )
        field_qc = {
            "status": "PASS" if field_identity_pass else "FAIL",
            "npz_path": str(npz_path),
            "npz_sha256": sha256_file(npz_path),
            "vtu_path": str(flow_vtu),
            "vtu_sha256": sha256_file(flow_vtu),
            "vtu_binary": vtu_binary,
            "cell_count": int(reloaded.n_cells),
            "point_count": int(reloaded.n_points),
            "cell_type": "VTK_HEXAHEDRON",
            "cell_type_counts": {
                str(int(kind)): int(count)
                for kind, count in zip(cell_types, cell_type_counts, strict=True)
            },
            "velocity_phy_cell_data": "velocity_phy" in reloaded.cell_data,
            "velocity_phy_shape": list(velocity_reload.shape),
            "pressure_phy_cell_data": "pressure_phy" in reloaded.cell_data,
            "pressure_phy_shape": list(pressure_reload.shape),
            "pressure_gauge_pa_cell_data": "pressure_gauge_pa" in reloaded.cell_data,
            "coordinates_unit": "m",
            "coordinate_bounds_m": list(reloaded.bounds),
            "dx_m": EXPECTED_DX_M,
            "finite": field_finite,
            "velocity_m_s": percentile_summary(speed_reload, (50, 95, 99)),
            "pressure_gauge_pa": percentile_summary(gauge_reload, (1, 50, 99)),
        }
        write_json(layout.qc / "field_identity_qc.json", field_qc)
        if not field_identity_pass:
            raise FlowError(DIRECT_DECODE_INVALID, "VTU field identity contract failed")

        actual_mach = float(np.max(speed_reload) * EXPECTED_DT_S / EXPECTED_DX_M / math.sqrt(CS2))
        mach_qc = {
            "status": "PASS" if actual_mach < MAXIMUM_LATTICE_MACH else "FAIL",
            "maximum_velocity_m_s": float(np.max(speed_reload)),
            "actual_lattice_mach": actual_mach,
            "maximum_allowed": MAXIMUM_LATTICE_MACH,
            "dx_m": EXPECTED_DX_M,
            "dt_s": EXPECTED_DT_S,
            "d3q19_cs": math.sqrt(CS2),
        }
        write_json(layout.qc / "mach_qc.json", mach_qc)
        if mach_qc["status"] != "PASS":
            raise FlowError(DIRECT_DECODE_INVALID, "Actual lattice Mach is not below 0.05")

        inputs = load_flow_inputs(config.paths.source_surface_run)
        partition = load_frozen_surface_partition(
            inputs,
            frozen_seeder / "geometry" / "geometry_solver_m",
        )
        bc = load_boundary_conditions(inputs.boundary_conditions)
        metadata = build_proteus_metadata(
            flow_vtu=flow_vtu,
            inlet_area_m2=partition.patch("inlet").area_um2 * 1.0e-12,
            dx_m=EXPECTED_DX_M,
        )
        metadata_path = layout.proteus / "proteus_flow_metadata.json"
        write_json(metadata_path, metadata)
        proteus_qc = {
            "status": "PASS",
            "coordinates_meter": True,
            "cartesian_hexahedral_topology": True,
            "uniform_dx_m": EXPECTED_DX_M,
            "velocity_phy_cell_data": True,
            "velocity_components": 3,
            "velocity_unit": "m/s",
            "pressure_phy_cell_data": True,
            "finite_values": field_finite,
            "inlet_normal": None,
            "inlet_normal_policy": "AUTO_DETECT_BY_BACKPROPAGATION_LATER",
            "proteus_runtime_executed": False,
        }
        write_json(layout.qc / "proteus_field_contract_qc.json", proteus_qc)

        port_failure: FlowError | None = None
        try:
            fluxes, measured_pressures = numerical_port_fluxes(
                reloaded,
                partition,
                EXPECTED_DX_M,
            )
            q_in = float(fluxes["inlet"])
            q_out = [float(fluxes[f"outlet_{index:02d}"]) for index in range(1, 4)]
            inlet_error = abs(q_in - bc.inlet_flow_m3_s) / abs(bc.inlet_flow_m3_s)
            mass_error = abs(abs(q_in) - sum(abs(value) for value in q_out)) / abs(q_in)
            directions_pass = bool(q_in > 0.0 and all(value > 0.0 for value in q_out))
            flux_pass = inlet_error <= 0.01 and mass_error <= 0.01 and directions_pass
            flux_qc = {
                "status": "PASS" if flux_pass else "FAIL",
                "failure_kind": None if flux_pass else "PORT_INTEGRATION_OR_ACTUAL_FLOW_REVIEW",
                "method": "existing numerical_port_fluxes",
                "q_in_m3_s": q_in,
                "q_out_m3_s": {
                    f"outlet_{index:02d}": value
                    for index, value in enumerate(q_out, start=1)
                },
                "inlet_target_m3_s": bc.inlet_flow_m3_s,
                "inlet_relative_error": inlet_error,
                "mass_conservation_error": mass_error,
                "flow_directions_pass": directions_pass,
                "maximum_allowed_inlet_relative_error": 0.01,
                "maximum_allowed_mass_conservation_error": 0.01,
                "port_diagnostics": _port_diagnostics(
                    fluxes,
                    partition,
                    bc.kinematic_viscosity_m2_s,
                ),
                "expected_1d_vs_3d_outlet_flow_diagnostic_only": [
                    {
                        "label": f"outlet_{index:02d}",
                        "expected_1d_m3_s": expected,
                        "measured_3d_m3_s": q_out[index - 1],
                        "role": "DIAGNOSTIC_ONLY",
                    }
                    for index, expected in enumerate(
                        bc.outlet_expected_1d_flows_m3_s,
                        start=1,
                    )
                ],
            }
            pressure_qc = {
                "status": "DIAGNOSTIC",
                "method": "boundary-adjacent internal cap plane from numerical_port_fluxes",
                "outlets": [
                    {
                        "label": f"outlet_{index:02d}",
                        "target_gauge_pa": target,
                        "measured_boundary_adjacent_gauge_pa": measured_pressures[
                            f"outlet_{index:02d}"
                        ],
                        "difference_pa": measured_pressures[f"outlet_{index:02d}"] - target,
                        "role": "DIAGNOSTIC_ONLY",
                    }
                    for index, target in enumerate(bc.outlet_gauge_pressures_pa, start=1)
                ],
            }
            if not flux_pass:
                port_failure = FlowError(PORT_REVIEW_STATUS, "Port flux gates require review")
        except FlowError as error:
            flux_qc = {
                "status": "FAIL",
                "failure_kind": "PORT_INTEGRATION_SAMPLING_FAILURE",
                "macro_decode_remains_valid": True,
                "failure": str(error),
            }
            pressure_qc = {
                "status": "NOT_EVALUATED_DUE_TO_PORT_SAMPLING_FAILURE",
            }
            port_failure = FlowError(PORT_REVIEW_STATUS, str(error))
        write_json(layout.qc / "port_flux_qc.json", flux_qc)
        write_json(layout.qc / "outlet_pressure_diagnostic.json", pressure_qc)

        harvester_review = _historical_harvester_review(output_root)
        write_json(layout.qc / "harvester_historical_issue_review.json", harvester_review)

        figures: list[str] = []
        if port_failure is None:
            velocity_figure = layout.figures / "velocity_magnitude_review.png"
            pressure_figure = layout.figures / "gauge_pressure_review.png"
            _scatter_figure(
                centers_m,
                speed_reload,
                velocity_figure,
                title="Direct-decoded steady velocity magnitude",
                colorbar_label="|u| (m/s)",
            )
            _scatter_figure(
                centers_m,
                gauge_reload,
                pressure_figure,
                title="Direct-decoded steady gauge pressure",
                colorbar_label="gauge pressure (Pa)",
            )
            figures = [str(velocity_figure), str(pressure_figure)]

        source_audit, critical_after = _source_read_only_audit(
            paths=critical_paths,
            before=critical_before,
            destination=layout.qc / "frozen_source_read_only_qc.json",
        )
        if source_audit["status"] != "PASS":
            raise FlowError(DIRECT_DECODE_INVALID, "Frozen source files changed during decode")

        direct_qc = {
            "status": "PASS",
            "binary_sha256": sha256_file(header.binary_path),
            "binary_bytes": header.binary_path.stat().st_size,
            "pdf_shape": list(pdf.shape),
            "pdf_dtype": str(pdf.dtype),
            "pdf_finite": finite_pass,
            "d3q19_source_contract": source_contract["d3q19"],
            "mean_pressure_validation": pressure_validation,
            "mean_speed_validation": speed_validation,
            "total_density_validation": density_validation,
            "tree_id_validation": mesh_qc,
            "mesh_order_proof": source_contract["restart_mesh_index_correspondence"],
            "vtu_field_identity": field_qc,
            "mach": mach_qc,
        }
        write_json(layout.qc / "direct_restart_decode_qc.json", direct_qc)

        final_status = SUCCESS_STATUS if port_failure is None else PORT_REVIEW_STATUS
        final_next = SUCCESS_NEXT if port_failure is None else PORT_REVIEW_NEXT
        summary.update(
            {
                "status": final_status,
                "next": final_next,
                "restart_serialization_contract": "PASS",
                "treelm_elemlist_binary_contract": "PASS",
                "macro_reconstruction_validation": "PASS",
                "pdf_dtype": str(pdf.dtype),
                "pdf_shape": list(pdf.shape),
                "pdf_finite": finite_pass,
                "macro_validation": macro_qc,
                "mesh_decode": mesh_qc,
                "direct_cell_field_npz": str(npz_path),
                "flow_vtu": str(flow_vtu),
                "field_identity": field_qc,
                "mach": mach_qc,
                "port_flux": flux_qc,
                "outlet_pressure": pressure_qc,
                "proteus_field_contract": proteus_qc,
                "proteus_metadata": str(metadata_path),
                "harvester_historical_issue": harvester_review,
                "figures": figures,
                "critical_source_manifest_after": critical_after,
                "source_steady_restart_modified": False,
                "grid_convergence": "NOT_RUN",
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        if critical_paths and critical_before:
            try:
                audit, after = _source_read_only_audit(
                    paths=critical_paths,
                    before=critical_before,
                    destination=layout.qc / "frozen_source_read_only_qc.json",
                )
                summary["critical_source_manifest_after"] = after
                summary["source_steady_restart_modified"] = audit["status"] != "PASS"
            except Exception as audit_error:  # pragma: no cover - last-resort evidence path
                summary["source_read_only_audit_failure"] = str(audit_error)
        status = error.status if isinstance(error, FlowError) else DIRECT_DECODE_INVALID
        summary.update(
            {
                "status": status,
                "failure": str(error),
                "next": "REVIEW DIRECT RESTART DECODE EVIDENCE WITHOUT APES FALLBACK",
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
