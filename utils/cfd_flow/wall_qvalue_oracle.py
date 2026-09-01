"""Zero-run analytic audit of Seeder q-values consumed by Musubi wall_libb.

This research helper reads existing TreElm meshes only.  It does not launch
Seeder, Musubi, Harvester, or any vascular workflow.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .io import sha256_file


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
