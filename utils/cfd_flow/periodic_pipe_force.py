"""Build and audit frozen-geometry periodic Poiseuille Pipe Force cases.

This module is research-only.  It emits transparent Seeder/Musubi inputs and
self-locating launchers; solver execution remains outside Python so every WSL,
MPI, Seeder and Musubi invocation has a durable script and log.
"""

from __future__ import annotations

import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .io import sha256_file, write_json
from .musubi_pressure_bc_benchmark import (
    BULK_NU_M2_S,
    MEAN_VELOCITY_M_S,
    NU_M2_S,
    PIPE_DIAMETER_M,
    PIPE_RADIUS_M,
    PRESSURE_REFERENCE_PA,
    RHO_KG_M3,
    diffusive_time_step,
    generate_pipe_wall,
    pipe_frame,
)
from .musubi_wall_geometry_benchmark import (
    early_stop_decision,
    effective_radius_ratio,
    poiseuille_pressure_gradient_pa_m,
    source_force_vector,
)
from .qvalue_repair_validation import tiny_cylinder_gate


SEEDER_BASE_SHA = "667109df6fafdcb39f4409e3f5d90f04d75cd33c"
TREELM_BASE_SHA = "53f273dbb8e9dcbe7feeb3d9831a35f5ae3cd72c"
PATCH_SHA256 = "c6b005650e08a9fb0d928340050bfe1dbcdf0d52ee51a5fc72831420e894c368"
SEEDER_SHA256 = "d7be681ca90da706559a4fd7e8f769fdb8f4303b8508f751077205f8e00cc7ed"
MUSUBI_SHA256 = "e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588"
MUSUBI_WSL = (
    "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/"
    "build/musubi_adaptive_flux"
)
MPIRUN_WSL = "/home/lzy/.local/bin/mpirun"
ROOT_LEVEL = 8
TRACKING_INTERVAL = 50
MAXIMUM_ITERATIONS = 2000


@dataclass(frozen=True, slots=True)
class PeriodicPipeCase:
    name: str
    cells_across: int
    axis: tuple[float, float, float]
    translation_cells: tuple[int, int, int]

    @property
    def dx_m(self) -> float:
        return PIPE_DIAMETER_M / self.cells_across

    @property
    def dt_s(self) -> float:
        return diffusive_time_step(self.dx_m)

    @property
    def direction(self) -> np.ndarray:
        value = np.asarray(self.axis, dtype=np.float64)
        return value / np.linalg.norm(value)

    @property
    def translation_m(self) -> np.ndarray:
        return self.dx_m * np.asarray(self.translation_cells, dtype=np.float64)

    @property
    def length_m(self) -> float:
        return float(np.linalg.norm(self.translation_m))


CASES = {
    "axis_n16": PeriodicPipeCase("axis_n16", 16, (0.0, 0.0, 1.0), (0, 0, 96)),
    "axis_n20": PeriodicPipeCase("axis_n20", 20, (0.0, 0.0, 1.0), (0, 0, 120)),
    "axis_n27": PeriodicPipeCase("axis_n27", 27, (0.0, 0.0, 1.0), (0, 0, 162)),
    "oblique_n27": PeriodicPipeCase(
        "oblique_n27", 27, (1.0, 1.0, 0.0), (128, 128, 0)
    ),
}


def _lua_vector(values: Iterable[float]) -> str:
    return "{ " + ", ".join(f"{float(value):.17g}" for value in values) + " }"


def generate_periodic_seeder_lua(case: PeriodicPipeCase) -> str:
    axis, u, v = pipe_frame(case.direction)
    cube_length = (2**ROOT_LEVEL) * case.dx_m
    center = np.full(3, 0.5 * cube_length)
    first = center - 0.5 * case.translation_m
    second = center + 0.5 * case.translation_m
    halfwidth = PIPE_RADIUS_M + 2.0 * case.dx_m
    plane1_origin = first - halfwidth * u - halfwidth * v - axis * cube_length / 2**20
    plane2_origin = second - halfwidth * u - halfwidth * v + axis * cube_length / 2**20
    plane1_a = 2.0 * halfwidth * v
    plane1_b = 2.0 * halfwidth * u
    plane2_a = 2.0 * halfwidth * u
    plane2_b = 2.0 * halfwidth * v
    return f"""-- Frozen dimensionless-kernel periodic Pipe Force mesh.
folder = 'mesh/'
timing_file = 'seeder_timing.res'
comment = '{case.name}; wall_libb continuous q; periodic axial direction'
logging = {{ level = 3 }}
minlevel = {ROOT_LEVEL}
bounding_cube = {{ origin = {{ 0.0, 0.0, 0.0 }}, length = {cube_length:.17g} }}
spatial_object = {{
  {{ attribute = {{ kind = 'seed' }},
     geometry = {{ kind = 'canoND', object = {{ origin = {_lua_vector(center)} }} }} }},
  {{ attribute = {{ kind = 'boundary', label = 'wall', level = minlevel,
                     calc_dist = true, flood_diagonal = false }},
     geometry = {{ kind = 'stl', object = {{ filename = 'geometry/wall.stl' }} }} }},
  {{ attribute = {{ kind = 'periodic' }},
     geometry = {{ kind = 'periodic', object = {{
       plane1 = {{ origin = {_lua_vector(plane1_origin)},
                  vec = {{ {_lua_vector(plane1_a)}, {_lua_vector(plane1_b)} }} }},
       plane2 = {{ origin = {_lua_vector(plane2_origin)},
                  vec = {{ {_lua_vector(plane2_a)}, {_lua_vector(plane2_b)} }} }}
     }} }} }}
}}
"""


def generate_periodic_musubi_lua(case: PeriodicPipeCase) -> str:
    axis = case.direction
    force = source_force_vector(
        poiseuille_pressure_gradient_pa_m(
            mean_velocity_m_s=MEAN_VELOCITY_M_S,
            radius_m=PIPE_RADIUS_M,
            rho_kg_m3=RHO_KG_M3,
            nu_m2_s=NU_M2_S,
        ),
        axis,
    )
    cube_length = (2**ROOT_LEVEL) * case.dx_m
    center = np.full(3, 0.5 * cube_length)
    return f"""-- Frozen-geometry, force-driven periodic Poiseuille validation.
simulation_name = 'periodic_pipe_force_{case.name}'
printRuntimeInfo = true
timing_file = 'tracking/musubi_timing.res'
mesh = 'mesh/'
scaling = 'diffusive'
logging = {{ level = 5 }}
maximum_iterations = {MAXIMUM_ITERATIONS}
dx = {case.dx_m:.17g}
dt = {case.dt_s:.17g}
rho0_phy = {RHO_KG_M3:.17g}
nu_phy = {NU_M2_S:.17g}
bulk_viscosity_phy = {BULK_NU_M2_S:.17g}
pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}
radius = {PIPE_RADIUS_M:.17g}
target_mean = {MEAN_VELOCITY_M_S:.17g}
axis = {_lua_vector(axis)}
center = {_lua_vector(center)}

function radial_squared(x,y,z)
  local rx=x-center[1]; local ry=y-center[2]; local rz=z-center[3]
  local axial=rx*axis[1]+ry*axis[2]+rz*axis[3]
  rx=rx-axial*axis[1]; ry=ry-axial*axis[2]; rz=rz-axial*axis[3]
  return rx*rx+ry*ry+rz*rz
end
function vel_analy(x,y,z,t)
  return 2.0*target_mean*math.max(0.0,1.0-radial_squared(x,y,z)/(radius*radius))
end
function initial_velocity_x(x,y,z,t) return axis[1]*vel_analy(x,y,z,t) end
function initial_velocity_y(x,y,z,t) return axis[2]*vel_analy(x,y,z,t) end
function initial_velocity_z(x,y,z,t) return axis[3]*vel_analy(x,y,z,t) end

physics = {{ dt = dt, rho0 = rho0_phy }}
identify = {{ label = 'periodic_pipe_force', kind = 'fluid_incompressible', layout = 'd3q19', relaxation = 'bgk' }}
fluid = {{ kinematic_viscosity = nu_phy, bulk_viscosity = bulk_viscosity_phy }}
initial_condition = {{
  pressure = pressure_reference_phy,
  velocityX = initial_velocity_x, velocityY = initial_velocity_y, velocityZ = initial_velocity_z
}}
boundary_condition = {{ {{ label = 'wall', kind = 'wall_libb' }} }}
glob_source = {{ force = {_lua_vector(force)}, force_order = 2 }}

variable = {{
  {{ name = 'vel_analy', ncomponents = 1, vartype = 'st_fun', st_fun = vel_analy }},
  {{ name = 'vel_error', ncomponents = 1, vartype = 'operation',
     operation = {{ kind = 'difference', input_varname = {{'vel_mag_phy','vel_analy'}} }} }}
}}

sim_control = {{
  time_control = {{ max = {{ iter = maximum_iterations }}, interval = {{ iter = 50 }} }},
  abort_criteria = {{
    stop_file = 'stop', steady_state = true,
    convergence = {{
      variable = {{'vel_mag_phy','vel_error'}}, shape = {{kind='all'}},
      reduction = {{'average','l2norm'}},
      time_control = {{ min = {{iter=200}}, max = {{iter=maximum_iterations}}, interval = {{iter=50}} }},
      norm = 'average', nvals = 5, absolute = false,
      condition = {{
        {{threshold=1.0e-4, operator='<='}}, {{threshold=1.0e-4, operator='<='}}
      }}
    }}
  }}
}}

tracking = {{
  {{ label = 'mean_velocity', folder = 'tracking/', variable = {{'velocity_phy'}},
     shape = {{ kind = 'all' }}, reduction = {{'average'}},
     time_control = {{ min = {{iter=0}}, max = {{iter=maximum_iterations}}, interval = {{iter=50}} }},
     output = {{ format = 'ascii' }} }},
  {{ label = 'profile', folder = 'tracking/', variable = {{'vel_analy','vel_error'}},
     shape = {{ kind = 'all' }}, reduction = {{'l2norm','l2norm'}},
     time_control = {{ min = {{iter=0}}, max = {{iter=maximum_iterations}}, interval = {{iter=50}} }},
     output = {{ format = 'ascii' }} }},
  {{ label = 'safety', folder = 'tracking/', variable = {{'vel_mag','pdf'}},
     shape = {{ kind = 'all' }}, reduction = {{'max','min'}},
     time_control = {{ min = {{iter=0}}, max = {{iter=maximum_iterations}}, interval = {{iter=50}} }},
     output = {{ format = 'ascii' }} }}
}}
"""


def _launcher_text(kind: str) -> str:
    common = f"""#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=\"$(cd -- \"$(dirname -- \"${{BASH_SOURCE[0]}}\")\" && pwd -P)\"
cd \"$SCRIPT_DIR\"
SEEDER='/home/lzy/apes-worktrees/seeder_dimensionless_kernel_20260830/build/seeder'
MUSUBI='{MUSUBI_WSL}'
MPIRUN='{MPIRUN_WSL}'
FREEZE=\"$SCRIPT_DIR/../../../../qc/seeder_geometry_freeze.json\"
"""
    if kind == "preflight_seeder":
        return common + f"""[[ -f seeder.lua && -f geometry/wall.stl && -x \"$SEEDER\" ]]
[[ -f \"$FREEZE\" ]] && grep -q 'CFD_FLOW_SEEDER_GEOMETRY_KERNEL_VALIDATED' \"$FREEZE\"
[[ \"$(sha256sum \"$SEEDER\" | awk '{{print $1}}')\" == '{SEEDER_SHA256}' ]]
mkdir -p mesh
probe=\"mesh/.write_probe_$$\"; : > \"$probe\"; rm -- \"$probe\"
printf 'CONFIG=%s\\nBINARY_SHA256=%s\\nPREFLIGHT=PASS\\n' \"$SCRIPT_DIR/seeder.lua\" '{SEEDER_SHA256}'
"""
    if kind == "run_seeder":
        return common + """/bin/bash "$SCRIPT_DIR/preflight_seeder.sh" > seeder_preflight.log 2>&1
if find mesh -mindepth 1 -print -quit | grep -q .; then
  printf 'Refusing to overwrite non-empty mesh directory\n' >&2; exit 3
fi
"$SEEDER" seeder.lua > seeder_stdout.log 2> seeder_stderr.log
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "mesh/$f" ]]; done
grep -q 'Seeder created mesh successfully!' seeder_stdout.log
printf 'SEMANTIC_SUCCESS=PASS\n' > seeder_semantic_status.log
"""
    if kind == "preflight_musubi":
        return common + f"""[[ -f musubi.lua && -d mesh && -x \"$MUSUBI\" && -x \"$MPIRUN\" ]]
[[ -f \"$FREEZE\" ]] && grep -q 'CFD_FLOW_SEEDER_GEOMETRY_KERNEL_VALIDATED' \"$FREEZE\"
[[ \"$(sha256sum \"$MUSUBI\" | awk '{{print $1}}')\" == '{MUSUBI_SHA256}' ]]
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s \"mesh/$f\" ]]; done
mkdir -p tracking
probe=\"tracking/.write_probe_$$\"; : > \"$probe\"; rm -- \"$probe\"
printf 'CONFIG=%s\\nBINARY_SHA256=%s\\nPREFLIGHT=PASS\\n' \"$SCRIPT_DIR/musubi.lua\" '{MUSUBI_SHA256}'
"""
    if kind == "run_musubi":
        return common + """/bin/bash "$SCRIPT_DIR/preflight_musubi.sh" > musubi_preflight.log 2>&1
if find tracking -mindepth 1 -print -quit | grep -q .; then
  printf 'Refusing to overwrite non-empty tracking directory\n' >&2; exit 3
fi
"$MPIRUN" --bind-to core --map-by core -np 2 "$MUSUBI" musubi.lua > musubi_stdout.log 2> musubi_stderr.log
grep -q 'Initializing simulation' musubi_stdout.log
test "$(find tracking -type f -name '*.res' -size +0c | wc -l)" -ge 3
printf 'SEMANTIC_SUCCESS=PASS\n' > musubi_semantic_status.log
"""
    raise ValueError(kind)


def prepare_case(root: Path, name: str) -> dict[str, Any]:
    case = CASES[name]
    case_dir = root / "cases" / name
    geometry = case_dir / "geometry"
    geometry.mkdir(parents=True, exist_ok=True)
    cube_center = np.full(3, 0.5 * (2**ROOT_LEVEL) * case.dx_m)
    wall = generate_pipe_wall(
        geometry / "wall.stl",
        case.direction,
        radius_m=PIPE_RADIUS_M,
        length_m=case.length_m,
        sections=128,
    )
    # The shared generator is origin-centred; translate without changing faces.
    import trimesh

    surface = trimesh.load_mesh(geometry / "wall.stl", process=False)
    surface.apply_translation(cube_center)
    surface.export(geometry / "wall.stl")
    wall["sha256"] = sha256_file(geometry / "wall.stl")
    wall["center_m"] = cube_center.tolist()
    (case_dir / "seeder.lua").write_text(
        generate_periodic_seeder_lua(case), encoding="utf-8", newline="\n"
    )
    (case_dir / "musubi.lua").write_text(
        generate_periodic_musubi_lua(case), encoding="utf-8", newline="\n"
    )
    for kind in (
        "preflight_seeder",
        "run_seeder",
        "preflight_musubi",
        "run_musubi",
    ):
        (case_dir / f"{kind}.sh").write_text(
            _launcher_text(kind), encoding="utf-8", newline="\n"
        )
    source = {
        "case": name,
        "orientation": "oblique" if name.startswith("oblique") else "axis_aligned",
        "direction": case.direction.tolist(),
        "translation_cells": list(case.translation_cells),
        "cells_across_diameter": case.cells_across,
        "dx_m": case.dx_m,
        "dt_s": case.dt_s,
        "wall": wall,
        "seeder_base_sha": SEEDER_BASE_SHA,
        "treelm_base_sha": TREELM_BASE_SHA,
        "patch_sha256": PATCH_SHA256,
        "seeder_binary_sha256": SEEDER_SHA256,
        "musubi_binary_sha256": MUSUBI_SHA256,
    }
    write_json(case_dir / "source_mesh_summary.json", source)
    return source


def _read_ascii_result(path: Path) -> tuple[list[str], np.ndarray]:
    header: list[str] | None = None
    rows: list[list[float]] = []
    readable = path
    if os.name == "nt" and not str(path).startswith("\\\\?\\"):
        readable = Path("\\\\?\\" + str(path.resolve()))
    for raw in readable.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            tokens = line[1:].split()
            if tokens and tokens[0] == "time":
                header = tokens
            continue
        rows.append([float(value) for value in line.split()])
    if header is None or not rows:
        raise ValueError(f"missing result header/data: {path}")
    return header, np.asarray(rows, dtype=np.float64)


def _result_file(case_dir: Path, label: str) -> Path:
    matches = sorted((case_dir / "tracking").glob(f"*{label}*.res"))
    if len(matches) != 1:
        raise ValueError(f"expected one {label} result, found {matches}")
    return matches[0]


def _columns(header: list[str], values: np.ndarray, pattern: str) -> np.ndarray:
    indices = [index for index, name in enumerate(header) if re.search(pattern, name)]
    if not indices:
        raise ValueError(f"no columns /{pattern}/ in {header}")
    return values[:, indices]


def audit_case(case_dir: Path) -> dict[str, Any]:
    source = json.loads((case_dir / "source_mesh_summary.json").read_text(encoding="utf-8"))
    case = CASES[str(source["case"])]
    q_gate = tiny_cylinder_gate(case_dir / "mesh", case_dir / "source_mesh_summary.json")
    mh, mv = _read_ascii_result(_result_file(case_dir, "mean_velocity"))
    ph, pv = _read_ascii_result(_result_file(case_dir, "profile"))
    sh, sv = _read_ascii_result(_result_file(case_dir, "safety"))
    times = mv[:, 0]
    if not (np.allclose(times, pv[:, 0]) and np.allclose(times, sv[:, 0])):
        raise ValueError("tracking times do not agree")
    velocity = _columns(mh, mv, r"velocity_phy").reshape(len(mv), -1)[:, :3]
    analytic_l2 = _columns(ph, pv, r"vel_analy").reshape(len(pv), -1)[:, 0]
    error_l2 = _columns(ph, pv, r"vel_error").reshape(len(pv), -1)[:, 0]
    max_lattice = _columns(sh, sv, r"vel_mag(?!_phy)").reshape(len(sv), -1)[:, 0]
    pdf_min = np.min(_columns(sh, sv, r"pdf"), axis=1)
    iterations = np.rint(times / case.dt_s).astype(int)
    mean_axial = velocity @ case.direction
    profile_l2 = error_l2 / analytic_l2
    q_values = mean_axial * math.pi * PIPE_RADIUS_M**2
    samples = [
        {
            "iteration": int(iterations[i]),
            "mean_axial_velocity": float(mean_axial[i]),
            "profile_l2_error": float(profile_l2[i]),
            "flow_rate_m3_s": float(q_values[i]),
            "minimum_pdf": float(pdf_min[i]),
            "maximum_lattice_speed": float(max_lattice[i]),
            "all_finite": bool(
                np.all(
                    np.isfinite(
                        [mean_axial[i], profile_l2[i], q_values[i], pdf_min[i], max_lattice[i]]
                    )
                )
            ),
        }
        for i in range(len(times))
    ]
    accepted_index = len(samples) - 1
    steady = early_stop_decision(samples)
    for index in range(1, len(samples)):
        decision = early_stop_decision(samples[: index + 1])
        if decision["stop"]:
            accepted_index = index
            steady = decision
            break
    last = samples[accepted_index]
    target_q = MEAN_VELOCITY_M_S * math.pi * PIPE_RADIUS_M**2
    mean_error = abs(last["mean_axial_velocity"] - MEAN_VELOCITY_M_S) / MEAN_VELOCITY_M_S
    q_error = abs(last["flow_rate_m3_s"] - target_q) / target_q
    radius_ratio = effective_radius_ratio(last["flow_rate_m3_s"], target_q)
    safety = (
        bool(last["all_finite"])
        and float(last["minimum_pdf"]) > 0.0
        and float(last["maximum_lattice_speed"]) < 0.05
    )
    scientific_pass = (
        q_gate["status"] == "PASS"
        and steady["stop"]
        and safety
        and mean_error <= 0.02
        and q_error <= 0.02
        and float(last["profile_l2_error"]) <= 0.02
        and abs(radius_ratio - 1.0) <= 0.01
    )
    result = {
        "status": "PASS" if scientific_pass else "FAIL",
        "case": case.name,
        "configuration": source,
        "q_gate": q_gate,
        "tracking_samples": samples,
        "accepted_iteration": int(last["iteration"]),
        "early_stop": steady,
        "final": {
            **last,
            "target_mean_velocity_m_s": MEAN_VELOCITY_M_S,
            "target_flow_rate_m3_s": target_q,
            "mean_velocity_relative_error": float(mean_error),
            "flow_rate_relative_error": float(q_error),
            "effective_radius_ratio": float(radius_ratio),
            "effective_radius_bias": float(abs(radius_ratio - 1.0)),
        },
        "gates": {
            "continuous_q": q_gate["status"] == "PASS",
            "steady_window": steady["stop"],
            "numerical_safety": safety,
            "mean_velocity_le_2pct": mean_error <= 0.02,
            "flow_rate_le_2pct": q_error <= 0.02,
            "profile_l2_le_2pct": float(last["profile_l2_error"]) <= 0.02,
            "effective_radius_bias_le_1pct": abs(radius_ratio - 1.0) <= 0.01,
        },
    }
    write_json(case_dir / "case_status.json", result)
    return result


def final_decision(root: Path) -> dict[str, Any]:
    loaded = {
        path.parent.name: json.loads(path.read_text(encoding="utf-8"))
        for path in (root / "cases").glob("*/case_status.json")
    }
    required = {"axis_n16", "axis_n27"}
    if not required.issubset(loaded):
        raise ValueError(f"missing required cases: {sorted(required - loaded.keys())}")
    supplement = "oblique_n27" if loaded["axis_n27"]["status"] == "PASS" else "axis_n20"
    if supplement not in loaded:
        return {"status": "PENDING", "next_case": supplement, "cases": loaded}
    if supplement == "oblique_n27":
        axis_q = loaded["axis_n27"]["final"]["flow_rate_m3_s"]
        oblique_q = loaded["oblique_n27"]["final"]["flow_rate_m3_s"]
        spread = abs(oblique_q - axis_q) / abs(axis_q)
        passed = loaded["axis_n27"]["status"] == "PASS" and loaded[supplement]["status"] == "PASS" and spread <= 0.02
    else:
        spread = None
        passed = False
    archived_calls = sum(
        1 for path in (root / "cases").glob("*/attempt_*") if path.is_dir()
    )
    total_calls = len(loaded) + archived_calls
    result = {
        "status": "CFD_FLOW_WALL_BENCHMARK_PASS" if passed else "CFD_FLOW_WALL_BENCHMARK_FAILED",
        "cases": loaded,
        "scientific_call_ledger": {
            "accepted_case_calls": sorted(loaded),
            "archived_noncontract_calls": archived_calls,
        },
        "small_musubi_semantic_calls": total_calls,
        "small_musubi_hard_max": 4,
        "small_musubi_budget_respected": total_calls <= 4,
        "orientation_flow_spread": spread,
        "first_scientific_failure": None if passed else "PERIODIC_PIPE_FORCE_FINAL_GATE",
        "next": "continuous_q_referee_v2" if passed else "SCIENTIFIC_STOP",
    }
    write_json(root / "periodic_pipe_force_status.json", result)
    return result
