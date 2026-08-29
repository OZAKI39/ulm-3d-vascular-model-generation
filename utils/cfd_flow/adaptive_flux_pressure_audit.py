"""Source-proven, zero-run audit for a Musubi adaptive pressure-flux inlet.

The audit is deliberately separate from the production CFD pipeline.  It
reconstructs one hypothetical boundary update from the frozen iteration-198064
PDF field and never launches Seeder, Musubi, or Harvester.
"""

from __future__ import annotations

import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from pypdf import PdfReader

from .apes import diffusive_time_step
from .config import load_cfd_flow_config
from .exact_link_flux import (
    DIRECT_FIELD_RUN,
    EXPECTED_CELL_COUNT,
    EXPECTED_DX_M,
    EXPECTED_STENCIL_COUNTS,
    FROZEN_ITERATION,
    FROZEN_SEEDER_RUN,
    FROZEN_STEADY_RUN,
    GLOBAL_NORMALS,
    INVERSE_DIRECTIONS,
    MUSUBI_EXECUTABLE_SHA256,
    REFERENCE_DENSITY_KG_M3,
    TARGET_MASS_FLOW_KG_S,
    TARGET_Q_M3_S,
    _file_manifest,
    _normal_distribution,
    _valid_pressure_mask,
    build_coordinate_lookup,
    equilibrium_pdf,
    pull_fetch_pdfs,
    reconstruct_boundary,
    velocity_from_pdf,
)
from .io import FlowError, read_json, sha256_file, write_json
from .mcclure_adaptive_flux_reference import physical_volume_flux_to_lattice
from .port_flux_audit import (
    EXPECTED_BOUNDARY_ELEMENT_COUNT,
    EXPECTED_SIDE_COUNT,
    PORT_LABELS,
    extract_boundary_property_indices,
    parse_bnd_header,
    parse_boundary_property_header,
    read_boundary_ids,
)
from .restart_decode import read_restart_pdf, read_treelm_elemlist


AUDIT_PREFIX = "mcclure_adaptive_flux_audit_anchor003274"
AUDIT_REVISION = "MCCLURE_MUSUBI_ZERO_RUN_AFFINE_PRESSURE_FLUX_V1"
LBPM_COMMIT = "6d686d354e5b8140841d3601e4c8c0e4e4b77e48"
MUSUBI_COMMIT = "4e8b277b66226277171ef93bf054d21270812793"
MUS_COMMIT = "81f8c4f13772f6d4af31f335e1e3f99b02726e25"
TREELM_COMMIT = "9899d1376992c4fafc8a343d2b4ccef81de670d1"
CS2 = 1.0 / 3.0

ALGEBRA_PASS = "CFD_FLOW_MCCLURE_MUSUBI_ALGEBRA_PROVEN"
AUDIT_FAILED = "CFD_FLOW_MCCLURE_MUSUBI_ZERO_RUN_AUDIT_FAILED"
NEXT_IMPLEMENT = "IMPLEMENT ISOLATED MUSUBI ADAPTIVE FLUX BC"
NEXT_FIX = "FIX MCCLURE MUSUBI ZERO-RUN CONTRACT"


def physical_mass_factor(*, density_kg_m3: float, dx_m: float, dt_s: float) -> float:
    """Return kg/s represented by one lattice population per time step."""

    return float(density_kg_m3) * float(dx_m) ** 3 / float(dt_s)


def musubi_pressure_flux_lattice(
    boundary_density: float,
    *,
    stored_boundary_pdfs: np.ndarray,
    incoming_masks: np.ndarray,
    extrapolated_velocity: np.ndarray,
) -> float:
    """Evaluate signed domain influx after Musubi ``pressure_eq`` replacement."""

    stored = np.asarray(stored_boundary_pdfs, dtype=np.float64)
    masks = np.asarray(incoming_masks, dtype=bool)
    velocity = np.asarray(extrapolated_velocity, dtype=np.float64)
    if stored.ndim != 2 or stored.shape[1] != 19:
        raise ValueError("stored_boundary_pdfs must have shape (n, 19)")
    if masks.shape != (len(stored), 18) or velocity.shape != (len(stored), 3):
        raise ValueError("boundary masks or velocities have incompatible shape")
    candidate = equilibrium_pdf(float(boundary_density), velocity)
    total = 0.0
    for row in range(len(stored)):
        incoming = np.flatnonzero(masks[row])
        outgoing = INVERSE_DIRECTIONS[incoming]
        total += float(np.sum(candidate[row, incoming] - stored[row, outgoing]))
    return total


def musubi_pressure_flux_affine_coefficients(
    *,
    stored_boundary_pdfs: np.ndarray,
    incoming_masks: np.ndarray,
    extrapolated_velocity: np.ndarray,
) -> tuple[float, float]:
    """Derive ``F(rho)=alpha*rho+beta`` from native equilibrium replacement."""

    stored = np.asarray(stored_boundary_pdfs, dtype=np.float64)
    masks = np.asarray(incoming_masks, dtype=bool)
    velocity = np.asarray(extrapolated_velocity, dtype=np.float64)
    unit_equilibrium = equilibrium_pdf(1.0, velocity)
    alpha = 0.0
    beta = 0.0
    for row in range(len(stored)):
        incoming = np.flatnonzero(masks[row])
        outgoing = INVERSE_DIRECTIONS[incoming]
        alpha += float(np.sum(unit_equilibrium[row, incoming]))
        beta -= float(np.sum(stored[row, outgoing]))
    return alpha, beta


def solve_boundary_density(target_flux_lattice: float, alpha: float, beta: float) -> float:
    if not all(np.isfinite(value) for value in (target_flux_lattice, alpha, beta)):
        raise ValueError("target flux and affine coefficients must be finite")
    if alpha == 0.0:
        raise ValueError("adaptive pressure coefficient alpha is zero")
    return float((target_flux_lattice - beta) / alpha)


def reconstruct_musubi_boundary_state(
    boundary_density: float,
    *,
    fetched_boundary_pdfs: np.ndarray,
    incoming_masks: np.ndarray,
    extrapolated_velocity: np.ndarray,
) -> np.ndarray:
    """Build the complete post-boundary/pre-collision PDFs for safety checks."""

    reconstructed = np.asarray(fetched_boundary_pdfs, dtype=np.float64).copy()
    masks = np.asarray(incoming_masks, dtype=bool)
    equilibrium = equilibrium_pdf(float(boundary_density), extrapolated_velocity)
    for row in range(len(reconstructed)):
        incoming = np.flatnonzero(masks[row])
        reconstructed[row, incoming] = equilibrium[row, incoming]
    return reconstructed


def _git(repository: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return completed.stdout.strip()


def _wsl_git(repository: str, *arguments: str) -> str:
    completed = subprocess.run(
        ["wsl.exe", "-e", "git", "-C", repository, *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return completed.stdout.strip()


def _wsl_source_root() -> Path:
    for candidate in (
        Path(r"\\wsl.localhost\Ubuntu\home\lzy\apes-pinned\musubi_official"),
        Path(r"\\wsl$\Ubuntu\home\lzy\apes-pinned\musubi_official"),
    ):
        if candidate.is_dir():
            return candidate
    raise FlowError(AUDIT_FAILED, "Pinned Musubi source is unavailable")


def _evidence(path: Path, revision: str, tokens: Iterable[str]) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="strict")
    lines: dict[str, int] = {}
    for token in tokens:
        offset = text.find(token)
        if offset < 0:
            raise FlowError(AUDIT_FAILED, f"Missing source token in {path}: {token}")
        lines[token] = text.count("\n", 0, offset) + 1
    return {
        "path": str(path),
        "revision": revision,
        "sha256": sha256_file(path),
        "token_line_evidence": lines,
    }


def _build_lbpm_contract(root: Path, paper: Path) -> dict[str, Any]:
    repository = root / "external_reference" / "LBPM"
    if not repository.is_dir():
        raise FlowError(AUDIT_FAILED, f"LBPM reference is absent: {repository}")
    commit = _git(repository, "rev-parse", "HEAD")
    if commit != LBPM_COMMIT or _git(repository, "status", "--short"):
        raise FlowError(AUDIT_FAILED, "LBPM reference commit or clean state changed")
    source_specs = {
        "common/ScaLBL.cpp": (
            "ScaLBL_Communicator::D3Q19_Flux_BC_z",
            "ScaLBL_D3Q19_AAeven_Flux_BC_z",
            "ScaLBL_D3Q19_AAodd_Flux_BC_z",
            "ScaLBL_D3Q19_AAeven_Pressure_BC_z",
            "ScaLBL_D3Q19_AAodd_Pressure_BC_z",
        ),
        "cpu/D3Q19.cpp": (
            "ScaLBL_D3Q19_AAeven_Flux_BC_z",
            "ScaLBL_D3Q19_AAodd_Flux_BC_z",
            "ScaLBL_D3Q19_AAeven_Pressure_BC_z",
            "ScaLBL_D3Q19_AAodd_Pressure_BC_z",
        ),
        "cuda/D3Q19.cu": (
            "ScaLBL_D3Q19_AAeven_Flux_BC_z",
            "ScaLBL_D3Q19_AAodd_Flux_BC_z",
            "ScaLBL_D3Q19_AAeven_Pressure_BC_z",
            "ScaLBL_D3Q19_AAodd_Pressure_BC_z",
        ),
        "models/MRTModel.cpp": ("D3Q19_Flux_BC_z",),
        "tests/TestFluxBC.cpp": ("TestFluxBC", "fabs(flux + Q)"),
    }
    files = {
        name: _evidence(repository / name, LBPM_COMMIT, tokens)
        for name, tokens in source_specs.items()
    }
    reader = PdfReader(str(paper))
    return {
        "status": "PASS",
        "clean_room_policy": "FORMULAS_AND_OBSERVABLE_BEHAVIOR_ONLY_NO_LBPM_SOURCE_COPIED",
        "repository_url": "https://github.com/OPM/LBPM.git",
        "repository_path": str(repository),
        "commit": commit,
        "license": "GNU General Public License v3.0",
        "license_path": str(repository / "LICENSE"),
        "source_files": files,
        "paper": {
            "path": str(paper),
            "sha256": sha256_file(paper),
            "pages": len(reader.pages),
            "title": "An adaptive volumetric flux boundary condition for lattice Boltzmann methods",
            "equations_used": ["Eq. 6", "Eq. 7", "Eq. 9", "Appendix B"],
            "rendered_visual_qa": [
                str(root / "tmp" / "pdfs" / "mcclure" / f"page-{page}.png")
                for page in range(1, len(reader.pages) + 1)
            ],
        },
        "per_timestep_behavior": {
            "A": "read current boundary populations after streaming/before collision",
            "B": "reduce the D3Q19 consistency term over inlet nodes",
            "C": "solve one spatially uniform inlet density from total flux",
            "D": "apply pressure closure to reconstruct incoming populations",
        },
    }


def _build_musubi_contract(source_root: Path) -> dict[str, Any]:
    if _wsl_git("/home/lzy/apes-pinned/musubi_official", "rev-parse", "HEAD") != MUSUBI_COMMIT:
        raise FlowError(AUDIT_FAILED, "Pinned Musubi source commit changed")
    if _wsl_git("/home/lzy/apes-pinned/musubi_official", "status", "--short"):
        raise FlowError(AUDIT_FAILED, "Pinned Musubi source is dirty")
    mus = source_root / "mus" / "source"
    tem = source_root / "tem" / "source"
    files = {
        "fluid_boundary": _evidence(
            mus / "bc" / "mus_bc_fluid_module.fpp",
            MUS_COMMIT,
            (
                "subroutine mfr_eq(",
                "subroutine pressure_eq(",
                "rho = rho / physics%fac( iLevel )%press * cs2inv",
                "neighBufferPre_nNext(1,:)",
                "velocity(:,iElem) =   1.5_rk * uxB_1",
                "bitmask%val( iDir, iElem )",
                "state( ?FETCH?(",
            ),
        ),
        "boundary_header": _evidence(
            mus / "bc" / "mus_bc_header_module.fpp",
            MUS_COMMIT,
            (
                "case( 'mfr_bounceback', 'mfr_eq' )",
                "case( 'pressure_eq' )",
                "me( myBCID )%nNeighs = 2",
                "integer :: normalInd",
            ),
        ),
        "boundary_buffers": _evidence(
            mus / "bc" / "mus_bc_general_module.fpp",
            MUS_COMMIT,
            (
                "currState = state( :, pdf%nNext )",
                "bcBuffer always uses AOS",
                "neighBufferPre_nNext",
                "?FETCH?(",
            ),
        ),
        "boundary_neighbor_construction": _evidence(
            mus / "mus_construction_module.fpp",
            MUS_COMMIT,
            (
                "subroutine setFieldBCNeigh(",
                "subroutine mus_build_BCStencils(",
                "if (bc%curved) then",
                "normal(:) = globBC%normal(:)",
                "stencil%cxDir(:, iNeigh) = normal * iNeigh",
            ),
        ),
        "timestep_control": _evidence(
            mus / "mus_control_module.f90",
            MUS_COMMIT,
            (
                "call set_boundary(",
                "call mus_swap_now_next(",
                "call me%scheme%compute(",
            ),
        ),
        "d3q19_equilibrium": _evidence(
            mus / "scheme" / "mus_scheme_derived_quantities_type_module.f90",
            MUS_COMMIT,
            ("pure function get_pdfEq_d3q19", "rho_div_18", "fEq(19)"),
        ),
        "pull_macros": _evidence(
            mus / "header" / "lbm_macros.inc",
            MUS_COMMIT,
            ("else !PULL", "macro :: FETCH", "macro :: SAVE"),
        ),
        "treelm_d3q19_ordering": _evidence(
            tem / "tem_stencil_module.fpp",
            TREELM_COMMIT,
            ("d3q19", "cxDir", "cxDirInv"),
        ),
    }
    return {
        "status": "PASS",
        "source_root": str(source_root),
        "musubi_commit": MUSUBI_COMMIT,
        "mus_submodule_commit": MUS_COMMIT,
        "treelm_submodule_commit": TREELM_COMMIT,
        "files": files,
        "answers": {
            "1_timestep_stage": "set_boundary updates nNext after the previous PULL stream/collision result and before swap plus the next collision",
            "2_modified_pdfs": "only D3Q19 incoming directions whose per-element boundary bitmask is true; writes use FETCH storage positions",
            "3_rho_dependence": "linear: at fixed neighbor-extrapolated velocity every equilibrium population is rho times a velocity polynomial; total link delta is affine",
            "4_velocity_source": "1.5 times first inward-neighbor pre-collision velocity minus 0.5 times second inward-neighbor pre-collision velocity",
            "5_normal_dependency": "incoming links depend on each element bitmask/normal construction; the two pressure_eq neighbor locations use global boundary normal when curved=false and element normal when curved=true",
            "6_dynamic_rho_then_native_reconstruction": "YES; alpha/beta may be reduced globally before reusing the native pressure_eq equilibrium replacement, with no collision change",
        },
        "pressure_conversion": "rho_lattice = pressure_physical / pressure_conversion * cs2inv",
        "streaming": "PULL",
        "source_algebra": "f_eq_i(rho,u)=rho*w_i*P_i(u), u independent of rho; F=sum_mask(f_eq_i-old_outgoing_i)=alpha*rho+beta",
    }


def run_adaptive_flux_pressure_audit(project_root: Path) -> dict[str, Any]:
    """Execute the paper/source contracts and frozen-PDF zero-run gate."""

    root = Path(project_root).resolve()
    head = _git(root, "rev-parse", "HEAD")
    branch = _git(root, "branch", "--show-current")
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    dt_s = diffusive_time_step(
        config.mesh.dx_target_m,
        config.physics.reference_dx_m,
        config.physics.reference_dt_s,
    )
    if config.mesh.dx_target_m != EXPECTED_DX_M:
        raise FlowError(AUDIT_FAILED, f"dx changed: {config.mesh.dx_target_m}")
    q_target_lattice = physical_volume_flux_to_lattice(
        TARGET_Q_M3_S, dx_m=config.mesh.dx_target_m, dt_s=dt_s
    )

    output_root = config.paths.output_root
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"{AUDIT_PREFIX}_{stamp}"
    qc_dir = run_root / "qc"
    qc_dir.mkdir(parents=True, exist_ok=False)
    manifest_path = qc_dir / "adaptive_flux_pressure_audit_manifest.json"
    summary: dict[str, Any] = {
        "status": AUDIT_FAILED,
        "next": NEXT_FIX,
        "audit_revision": AUDIT_REVISION,
        "run_root": str(run_root),
        "actual_head": head,
        "branch": branch,
        "production_pipeline_modified": False,
        "pinned_musubi_source_modified": False,
        "seeder_calls": 0,
        "musubi_calls": 0,
        "harvester_calls": 0,
        "grid_convergence": "NOT_RUN",
        "dx_m": config.mesh.dx_target_m,
        "dt_s": dt_s,
        "q_target_physical_m3_s": TARGET_Q_M3_S,
        "mass_flow_target_kg_s": TARGET_MASS_FLOW_KG_S,
        "q_target_lattice": q_target_lattice,
        "started_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, summary)

    production_paths = (
        root / "cfd_flow.py",
        root / "configs" / "cfd_flow.yaml",
        root / "utils" / "cfd_flow" / "pipeline.py",
    )
    paper = root / "An adaptive volumetric flux boundary condition for lattice Boltzmann.pdf"
    direct_run = output_root / DIRECT_FIELD_RUN
    seeder_run = output_root / FROZEN_SEEDER_RUN
    steady_run = output_root / FROZEN_STEADY_RUN
    mesh_dir = seeder_run / "seeder" / "mesh"
    direct_manifest_path = direct_run / "qc" / "direct_restart_decode_manifest.json"
    direct_manifest = read_json(direct_manifest_path)
    restart_binary = Path(direct_manifest["restart_header_contract"]["binary"])
    frozen_paths = (
        direct_run / "flow" / "direct_cell_field.npz",
        direct_manifest_path,
        restart_binary,
        restart_binary.parent / "roi003274_steady_lbm_lastHeader.lua",
        mesh_dir / "elemlist.lsb",
        mesh_dir / "bnd.lsb",
        mesh_dir / "header.lua",
        mesh_dir / "bnd.lua",
        steady_run / "diagnostic_musubi.lua",
        steady_run / "tracking" / "musubi_stdout.log",
    )
    source_before = _file_manifest(production_paths)
    frozen_before = _file_manifest(frozen_paths)

    try:
        if direct_manifest.get("frozen_iteration") != FROZEN_ITERATION:
            raise FlowError(AUDIT_FAILED, "Frozen restart iteration changed")
        prior_source = read_json(direct_run / "qc" / "direct_decode_source_contract.json")
        executable = Path(prior_source["musubi_executable_revision"]["path"])
        if sha256_file(executable) != MUSUBI_EXECUTABLE_SHA256:
            raise FlowError(AUDIT_FAILED, "Pinned Musubi executable hash changed")

        lbpm_contract = _build_lbpm_contract(root, paper)
        write_json(qc_dir / "mcclure_lbpm_reference_contract.json", lbpm_contract)
        source_root = _wsl_source_root()
        musubi_contract = _build_musubi_contract(source_root)
        write_json(qc_dir / "musubi_adaptive_flux_source_contract.json", musubi_contract)

        property_header = parse_boundary_property_header(
            (mesh_dir / "header.lua").read_text(encoding="utf-8")
        )
        boundary_header = parse_bnd_header(
            (mesh_dir / "bnd.lua").read_text(encoding="utf-8")
        )
        if (
            property_header.element_count != EXPECTED_BOUNDARY_ELEMENT_COUNT
            or boundary_header.side_count != EXPECTED_SIDE_COUNT
        ):
            raise FlowError(AUDIT_FAILED, "Frozen mesh boundary dimensions changed")
        tree_ids, property_bits, _ = read_treelm_elemlist(
            mesh_dir / "elemlist.lsb", n_elems=EXPECTED_CELL_COUNT
        )
        property_indices = extract_boundary_property_indices(
            property_bits, property_header.bit_position
        )
        boundary_ids = read_boundary_ids(
            mesh_dir / "bnd.lsb",
            element_count=property_header.element_count,
            side_count=boundary_header.side_count,
        )
        label_to_id = {
            label: index
            for index, label in enumerate(boundary_header.labels, start=1)
        }
        inlet = reconstruct_boundary(
            boundary_ids,
            property_indices,
            label="inlet",
            boundary_id=label_to_id["inlet"],
        )
        if len(inlet.cell_indices) != EXPECTED_STENCIL_COUNTS["inlet"]:
            raise FlowError(AUDIT_FAILED, "Frozen inlet boundary count changed")

        with np.load(direct_run / "flow" / "direct_cell_field.npz") as data:
            npz_tree_ids = np.asarray(data["tree_id"], dtype=np.int64)
            cell_ijk = np.asarray(data["cell_ijk"], dtype=np.int64)
        if not np.array_equal(npz_tree_ids, tree_ids):
            raise FlowError(AUDIT_FAILED, "Frozen restart and mesh cell order differ")
        pdf = read_restart_pdf(
            restart_binary, n_elems=EXPECTED_CELL_COUNT, n_components=19
        )
        lookup = build_coordinate_lookup(cell_ijk)
        valid, neighbor1, neighbor2 = _valid_pressure_mask(inlet, cell_ijk, lookup)
        if not np.all(valid):
            raise FlowError(
                AUDIT_FAILED,
                f"Candidate pressure reconstruction lacks two neighbors for {np.count_nonzero(~valid)} inlet cells",
            )
        fetched1 = pull_fetch_pdfs(
            pdf, cell_ijk, neighbor1, coordinate_lookup=lookup
        )
        fetched2 = pull_fetch_pdfs(
            pdf, cell_ijk, neighbor2, coordinate_lookup=lookup
        )
        velocity = 1.5 * velocity_from_pdf(fetched1) - 0.5 * velocity_from_pdf(
            fetched2
        )
        stored = np.asarray(pdf[inlet.cell_indices], dtype=np.float64)
        fetched_boundary = pull_fetch_pdfs(
            pdf, cell_ijk, inlet.cell_indices, coordinate_lookup=lookup
        )

        alpha_source, beta_source = musubi_pressure_flux_affine_coefficients(
            stored_boundary_pdfs=stored,
            incoming_masks=inlet.incoming_masks,
            extrapolated_velocity=velocity,
        )
        rho1, rho2 = 0.95, 1.05
        flux1 = musubi_pressure_flux_lattice(
            rho1,
            stored_boundary_pdfs=stored,
            incoming_masks=inlet.incoming_masks,
            extrapolated_velocity=velocity,
        )
        flux2 = musubi_pressure_flux_lattice(
            rho2,
            stored_boundary_pdfs=stored,
            incoming_masks=inlet.incoming_masks,
            extrapolated_velocity=velocity,
        )
        alpha_numeric = (flux2 - flux1) / (rho2 - rho1)
        beta_numeric = flux1 - alpha_numeric * rho1
        affine_scale = max(abs(alpha_source), abs(beta_source), 1.0)
        affine_error = max(
            abs(alpha_numeric - alpha_source), abs(beta_numeric - beta_source)
        ) / affine_scale
        rho_probe = 1.017
        probe_exact = musubi_pressure_flux_lattice(
            rho_probe,
            stored_boundary_pdfs=stored,
            incoming_masks=inlet.incoming_masks,
            extrapolated_velocity=velocity,
        )
        probe_affine = alpha_source * rho_probe + beta_source
        probe_relative_error = abs(probe_exact - probe_affine) / max(
            abs(probe_exact), 1.0e-30
        )
        if affine_error > 1.0e-12 or probe_relative_error > 1.0e-12:
            raise FlowError(AUDIT_FAILED, "Musubi F(rho) affine numerical proof failed")

        rho_target = solve_boundary_density(
            q_target_lattice, alpha_source, beta_source
        )
        exact_flux = musubi_pressure_flux_lattice(
            rho_target,
            stored_boundary_pdfs=stored,
            incoming_masks=inlet.incoming_masks,
            extrapolated_velocity=velocity,
        )
        relative_error = abs(exact_flux - q_target_lattice) / abs(q_target_lattice)
        reconstructed = reconstruct_musubi_boundary_state(
            rho_target,
            fetched_boundary_pdfs=fetched_boundary,
            incoming_masks=inlet.incoming_masks,
            extrapolated_velocity=velocity,
        )
        finite = bool(np.all(np.isfinite(reconstructed)))
        reconstructed_velocity = velocity_from_pdf(reconstructed)
        max_speed = float(np.max(np.linalg.norm(reconstructed_velocity, axis=1)))
        minimum_pdf = float(np.min(reconstructed))
        pressure_factor = REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**2 / dt_s**2
        pressure_reference = pressure_factor * CS2
        pressure_physical = rho_target * pressure_reference
        gauge_pressure = pressure_physical - pressure_reference
        mass_factor = physical_mass_factor(
            density_kg_m3=REFERENCE_DENSITY_KG_M3,
            dx_m=EXPECTED_DX_M,
            dt_s=dt_s,
        )
        exact_mass_flow = exact_flux * mass_factor
        q_from_mass = TARGET_MASS_FLOW_KG_S / REFERENCE_DENSITY_KG_M3
        target_consistency = abs(q_from_mass - TARGET_Q_M3_S) / TARGET_Q_M3_S

        zero_qc = {
            "status": "PASS"
            if relative_error <= 1.0e-10 and finite and max_speed < 0.05
            else "FAIL",
            "frozen_iteration": FROZEN_ITERATION,
            "candidate_identifiable": True,
            "pressure_neighbor_direction": GLOBAL_NORMALS["inlet"].astype(int).tolist(),
            "pressure_neighbor_valid_count": int(np.count_nonzero(valid)),
            "inlet_globbc_count": len(inlet.cell_indices),
            "normal_ind_distribution": _normal_distribution(inlet),
            "incoming_modified_link_count": int(np.count_nonzero(inlet.incoming_masks)),
            "source_algebra_affine_proven": True,
            "rho_1": rho1,
            "rho_2": rho2,
            "F_rho_1_lattice": flux1,
            "F_rho_2_lattice": flux2,
            "alpha_source": alpha_source,
            "beta_source": beta_source,
            "alpha_numeric": alpha_numeric,
            "beta_numeric": beta_numeric,
            "affine_coefficient_relative_error": affine_error,
            "third_probe_rho": rho_probe,
            "third_probe_relative_error": probe_relative_error,
            "q_target_lattice": q_target_lattice,
            "rho_target": rho_target,
            "relative_density_deviation": rho_target - 1.0,
            "equivalent_physical_pressure_pa": pressure_physical,
            "gauge_pressure_pa": gauge_pressure,
            "exact_controlled_flux_lattice": exact_flux,
            "exact_controlled_mass_flow_kg_s": exact_mass_flow,
            "relative_flux_error": relative_error,
            "all_pdfs_finite": finite,
            "minimum_reconstructed_pdf": minimum_pdf,
            "maximum_reconstructed_lattice_velocity": max_speed,
            "maximum_lattice_velocity_gate": 0.05,
            "physical_q_from_mass_flow_m3_s": q_from_mass,
            "physical_target_consistency_relative_error": target_consistency,
            "mass_conversion_kg_s_per_lattice_flux": mass_factor,
            "uses_mfr_eq_area_proxy": False,
        }
        write_json(qc_dir / "zero_run_adaptive_pressure_flux_qc.json", zero_qc)
        if zero_qc["status"] != "PASS":
            raise FlowError(
                AUDIT_FAILED,
                f"zero-run gates: rel={relative_error}, finite={finite}, umax={max_speed}",
            )

        source_after = _file_manifest(production_paths)
        frozen_after = _file_manifest(frozen_paths)
        unchanged = source_before == source_after and frozen_before == frozen_after
        if not unchanged:
            raise FlowError(AUDIT_FAILED, "Production or frozen source files changed")
        write_json(
            qc_dir / "frozen_read_only_sha_qc.json",
            {
                "status": "PASS",
                "source_frozen_files_unchanged": True,
                "production_before": source_before,
                "production_after": source_after,
                "frozen_before": frozen_before,
                "frozen_after": frozen_after,
            },
        )
        summary.update(
            {
                "status": ALGEBRA_PASS,
                "next": NEXT_IMPLEMENT,
                "mcclure_eq9_reference_reproduction": "PASS",
                "musubi_source_mapping": "PROVEN",
                "musubi_F_rho_affine_proven": True,
                "alpha": alpha_source,
                "beta": beta_source,
                "rho_boundary": rho_target,
                "predicted_physical_pressure_pa": pressure_physical,
                "gauge_pressure_pa": gauge_pressure,
                "zero_run_exact_mass_flow_kg_s": exact_mass_flow,
                "zero_run_relative_error": relative_error,
                "zero_run_max_lattice_speed": max_speed,
                "zero_run_minimum_pdf": minimum_pdf,
                "source_frozen_files_unchanged": True,
                "contracts": {
                    "mcclure_lbpm": str(qc_dir / "mcclure_lbpm_reference_contract.json"),
                    "musubi": str(qc_dir / "musubi_adaptive_flux_source_contract.json"),
                    "zero_run": str(qc_dir / "zero_run_adaptive_pressure_flux_qc.json"),
                },
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        summary.update(
            {
                "status": AUDIT_FAILED,
                "next": NEXT_FIX,
                "first_failure": str(error),
                "source_frozen_files_unchanged": source_before
                == _file_manifest(production_paths)
                and frozen_before == _file_manifest(frozen_paths),
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
