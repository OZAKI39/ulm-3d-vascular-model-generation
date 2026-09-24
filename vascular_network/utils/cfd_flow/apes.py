"""Official APES configuration generation and WSL2 process execution."""

from __future__ import annotations

import math
import os
import re
import shlex
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .config import ApesConfig, FlowConfig
from .geometry import BoundingCube, SeedPoint, SurfacePartition
from .io import FlowError, RunLayout, write_json


@dataclass(frozen=True, slots=True)
class BoundaryConditions:
    density_kg_m3: float
    kinematic_viscosity_m2_s: float
    dynamic_viscosity_pa_s: float
    inlet_port_id: str
    inlet_flow_m3_s: float
    inlet_profile_requested: str
    outlet_port_ids: tuple[str, ...]
    outlet_gauge_pressures_pa: tuple[float, ...]
    outlet_expected_1d_flows_m3_s: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class LatticeScaling:
    dx_m: float
    dt_s: float
    nu_lattice: float
    tau: float
    omega: float
    velocity_mean_m_s: float
    velocity_max_expected_m_s: float
    velocity_max_expected_lattice: float
    mach_max_expected: float
    pressure_conversion_pa: float
    pressure_reference_pa: float
    outlet_absolute_pressures_pa: tuple[float, ...]
    outlet_lattice_densities: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class ApesEnvironment:
    status: str
    execution_environment: str
    wsl_distribution: str
    binaries: dict[str, str | None]
    compilers: dict[str, str | None]
    mpi_ranks: int
    cpu_count: int
    available_ram_bytes: int
    commands: tuple[dict[str, Any], ...]


@dataclass(frozen=True, slots=True)
class ExternalRun:
    returncode: int
    command: tuple[str, ...]
    stdout_path: Path
    stderr_path: Path
    wall_time_s: float


def load_boundary_conditions(value: dict[str, Any]) -> BoundaryConditions:
    """Validate and load only the established steady Newtonian BC model."""

    fluid = value.get("fluid")
    model = value.get("flow_model")
    inlet = value.get("inlet")
    outlets = value.get("outlets")
    if not isinstance(fluid, dict) or not isinstance(model, dict):
        raise ValueError("Invalid fluid/flow_model in boundary conditions")
    if not isinstance(inlet, dict) or not isinstance(outlets, list) or len(outlets) != 3:
        raise ValueError("Boundary conditions must contain one inlet and three outlets")
    fixed = {
        "model": "NEWTONIAN",
        "density_kg_m3": 1056.0,
        "kinematic_viscosity_m2_s": 3.27e-6,
        "dynamic_viscosity_pa_s": 0.00345312,
    }
    for key, expected in fixed.items():
        if fluid.get(key) != expected:
            raise ValueError(f"fluid.{key} must remain {expected!r}")
    expected_model = {
        "steady": True,
        "incompressible": True,
        "rigid_wall": True,
        "no_slip_wall": True,
        "turbulence": False,
        "pulsatility": False,
    }
    for key, expected in expected_model.items():
        if model.get(key) is not expected:
            raise ValueError(f"flow_model.{key} must remain {expected!r}")
    if inlet.get("type") != "VOLUMETRIC_FLOW_RATE" or inlet.get("profile") != "PARABOLIC":
        raise ValueError("The frozen inlet must request volumetric flow with PARABOLIC profile")
    pressures: list[float] = []
    expected_flows: list[float] = []
    port_ids: list[str] = []
    for outlet in outlets:
        if not isinstance(outlet, dict) or outlet.get("type") != "PRESSURE_DIRICHLET":
            raise ValueError("Every outlet must be PRESSURE_DIRICHLET")
        if "P_solver_boundary_pa" not in outlet:
            raise ValueError("Outlet is missing corrected P_solver_boundary_pa")
        port_ids.append(str(outlet["port_id"]))
        pressures.append(float(outlet["P_solver_boundary_pa"]))
        expected_flows.append(float(outlet["Q_expected_1D_m3_s"]))
    return BoundaryConditions(
        density_kg_m3=float(fluid["density_kg_m3"]),
        kinematic_viscosity_m2_s=float(fluid["kinematic_viscosity_m2_s"]),
        dynamic_viscosity_pa_s=float(fluid["dynamic_viscosity_pa_s"]),
        inlet_port_id=str(inlet["port_id"]),
        inlet_flow_m3_s=float(inlet["Q_solver_m3_s"]),
        inlet_profile_requested=str(inlet["profile"]),
        outlet_port_ids=tuple(port_ids),
        outlet_gauge_pressures_pa=tuple(pressures),
        outlet_expected_1d_flows_m3_s=tuple(expected_flows),
    )


def diffusive_time_step(dx_m: float, kinematic_viscosity_m2_s: float) -> float:
    """Tau1 scaling: choose nu_lattice=1/6 directly, without a reference pair."""

    return float(dx_m) ** 2 / (6.0 * float(kinematic_viscosity_m2_s))


def bulk_viscosity_from_kinematic(kinematic_viscosity_m2_s: float) -> float:
    """Pinned Musubi D3Q19 baseline: explicit bulk viscosity is 2/3 nu."""

    return (2.0 / 3.0) * kinematic_viscosity_m2_s


def compute_lattice_scaling(
    config: FlowConfig, bc: BoundaryConditions | None, inlet_area_m2: float
) -> LatticeScaling:
    """Compute the exact pinned single-level physical/lattice relationship."""

    dx = config.mesh.dx_m
    dt = diffusive_time_step(dx, config.physics.kinematic_viscosity_m2_s)
    nu_lat = config.physics.kinematic_viscosity_m2_s * dt / dx**2
    tau = nu_lat / config.physics.lattice_cs_squared + 0.5
    omega = 1.0 / tau
    mean_velocity = config.boundary.target_volume_flow_m3_s / inlet_area_m2
    maximum_velocity = 2.0 * mean_velocity
    maximum_lattice_velocity = maximum_velocity * dt / dx
    mach = maximum_lattice_velocity / math.sqrt(config.physics.lattice_cs_squared)
    pressure_factor = config.physics.density_kg_m3 * dx**2 / dt**2
    pressure_reference = pressure_factor * config.physics.lattice_cs_squared
    absolute = tuple(pressure_reference + value for value in config.boundary.outlet_gauge_pressures_pa)
    densities = tuple(value / (pressure_factor * config.physics.lattice_cs_squared) for value in absolute)
    if not 0.0 < omega < 2.0 or maximum_lattice_velocity >= config.solver.maximum_lattice_speed:
        raise FlowError(
            "CFD_FLOW_LBM_SCALING_INVALID",
            f"omega={omega:.12g}, expected_Ma={mach:.12g}",
        )
    if not all(value > 0.0 for value in densities):
        raise FlowError("CFD_FLOW_LBM_SCALING_INVALID", "Pressure offset produced non-positive lattice density")
    source_differences = np.subtract.outer(
        config.boundary.outlet_gauge_pressures_pa,
        config.boundary.outlet_gauge_pressures_pa,
    )
    absolute_differences = np.subtract.outer(absolute, absolute)
    if not np.allclose(source_differences, absolute_differences, rtol=0.0, atol=1.0e-9):
        raise FlowError("CFD_FLOW_LBM_SCALING_INVALID", "Common pressure offset did not preserve differences")
    return LatticeScaling(
        dx_m=dx,
        dt_s=dt,
        nu_lattice=nu_lat,
        tau=tau,
        omega=omega,
        velocity_mean_m_s=mean_velocity,
        velocity_max_expected_m_s=maximum_velocity,
        velocity_max_expected_lattice=maximum_lattice_velocity,
        mach_max_expected=mach,
        pressure_conversion_pa=pressure_factor,
        pressure_reference_pa=pressure_reference,
        outlet_absolute_pressures_pa=absolute,
        outlet_lattice_densities=densities,
    )


def _lua_number(value: float) -> str:
    return f"{float(value):.17g}"


def _lua_vector(values: np.ndarray) -> str:
    return "{ " + ", ".join(_lua_number(float(item)) for item in values) + " }"


def generate_seeder_lua(
    partition: SurfacePartition, seed: SeedPoint, cube: BoundingCube
) -> str:
    """Generate the official single-level, five-STL Seeder configuration."""

    objects: list[str] = []
    for patch in partition.patches:
        calc_dist = ", calc_dist = true" if patch.label == "wall" else ""
        objects.append(
            "  {\n"
            f"    attribute = {{ kind = 'boundary', label = '{patch.label}', level = minlevel{calc_dist} }},\n"
            "    geometry = { kind = 'stl', object = { filename = "
            f"'../geometry/geometry_solver_m/{patch.label}.stl' }}}}\n"
            "  }"
        )
    objects.append(
        "  {\n"
        "    attribute = { kind = 'seed' },\n"
        f"    geometry = {{ kind = 'canoND', object = {{ origin = {_lua_vector(seed.coordinates_m)} }} }}\n"
        "  }"
    )
    return (
        "-- Generated production configuration; official Seeder syntax.\n"
        "folder = 'mesh/'\n"
        "comment = 'ROI003274 uniform 0.20 um CFD lattice'\n"
        "debug = { debugMode = false, debugFiles = false, debugMesh = 'debug/' }\n"
        f"minlevel = {cube.level}\n"
        f"bounding_cube = {{ origin = {_lua_vector(cube.origin_m)}, length = {_lua_number(cube.side_m)} }}\n"
        "spatial_object = {\n"
        + ",\n".join(objects)
        + "\n}\n"
    )


def generate_musubi_lua(
    config: FlowConfig,
    partition: SurfacePartition | None,
    bc: BoundaryConditions | None,
    scaling: LatticeScaling,
    *,
    mesh_path: str = "../seeder/mesh/",
    maximum_iterations: int | None = None,
) -> str:
    """Render the validated Tau1/adaptive/continuous-q production contract."""

    outlet_functions: list[str] = []
    outlet_entries: list[str] = []
    for index, pressure in enumerate(scaling.outlet_absolute_pressures_pa, start=1):
        label = f"outlet_{index:02d}"
        outlet_functions.append(
            f"function {label}_pressure(x, y, z, t) return {_lua_number(pressure)} end"
        )
        outlet_entries.append(
            "  { label = '"
            + label
            + "', kind = 'pressure_eq', pressure = "
            + label
            + "_pressure }"
        )
    iterations = maximum_iterations or config.execution.fresh_maximum_iterations
    target_lattice = (
        config.boundary.target_mass_flow_kg_s
        / config.physics.density_kg_m3
        * scaling.dt_s
        / scaling.dx_m**3
    )
    return f"""-- Validated Tau1 production configuration.
-- pressure_reference_phy is an LBM numerical offset, not physiological pressure.
simulation_name = 'roi003274_production_tau1'
printRuntimeInfo = true
timing_file = 'tracking/timing.res'
mesh = '{mesh_path}'
scaling = 'diffusive'
logging = {{ level = 5 }}

dx = {_lua_number(scaling.dx_m)}
dt = {_lua_number(scaling.dt_s)}
rho0_phy = {_lua_number(config.physics.density_kg_m3)}
nu_phy = {_lua_number(config.physics.kinematic_viscosity_m2_s)}
bulk_viscosity_phy = {_lua_number(config.physics.bulk_viscosity_m2_s)}
pressure_reference_phy = {_lua_number(scaling.pressure_reference_pa)}
target_lattice_flux_expected = {_lua_number(target_lattice)}
maximum_iterations = {iterations}

{chr(10).join(outlet_functions)}

sim_control = {{
  time_control = {{ max = {{ iter = maximum_iterations }}, interval = {{ iter = 100 }} }},
  abort_criteria = {{ stop_file = 'stop' }}
}}

physics = {{ dt = dt, rho0 = rho0_phy }}
identify = {{ label = 'ROI003274', kind = 'fluid', layout = 'd3q19', relaxation = 'bgk' }}
fluid = {{
  kinematic_viscosity = nu_phy,
  bulk_viscosity = bulk_viscosity_phy
}}
initial_condition = {{ pressure = pressure_reference_phy, velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0 }}

boundary_condition = {{
  {{ label = 'wall', kind = 'wall_libb' }},
  {{ label = 'inlet', kind = 'adaptive_flux_pressure', mass_flowrate = {_lua_number(config.boundary.target_mass_flow_kg_s)} }},
{','.join(chr(10) + entry for entry in outlet_entries)}
}}

restart = {{
  write = 'restart/',
  timeformat = {{ use_iter = true }},
  time_control = {{ min = {{ iter = maximum_iterations }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = 1 }} }}
}}
"""


def generate_harvester_lua() -> str:
    return """-- One final full-volume harvest from the terminal restart.
require 'musubi'
restart = { read = 'restart/roi003274_steady_lbm_lastHeader.lua' }
tracking = {
  label = 'flow_field',
  folder = '../flow/',
  variable = { 'pressure_phy', 'velocity_phy' },
  shape = { kind = 'all' },
  output = { format = 'vtk' }
}
"""


def _run_probe(distribution: str, script: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["wsl.exe", "-d", distribution, "--", "bash", "-lc", script],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )


def inspect_apes_environment(config: ApesConfig) -> ApesEnvironment:
    """Inspect and SHA-pin WSL2 tools; never fall back to a same-name binary."""

    commands: list[dict[str, Any]] = []
    names = {
        "seeder": config.seeder_executable,
        "musubi": config.musubi_executable,
        "mus_harvesting": config.harvesting_executable,
        "mpi_launcher": config.mpi_launcher,
        "lua_compiler": "luac",
    }
    binaries: dict[str, str | None] = {}
    for label, name in names.items():
        probe = _run_probe(config.wsl_distribution, f"command -v {shlex.quote(name)}")
        commands.append({"command": ["command", "-v", name], "returncode": probe.returncode, "stdout": probe.stdout.strip(), "stderr": probe.stderr.strip()})
        binaries[label] = probe.stdout.strip() or None
    compilers: dict[str, str | None] = {}
    for name in ("mpifort", "mpif90", "gfortran"):
        probe = _run_probe(config.wsl_distribution, f"command -v {name}")
        commands.append({"command": ["command", "-v", name], "returncode": probe.returncode, "stdout": probe.stdout.strip(), "stderr": probe.stderr.strip()})
        compilers[name] = probe.stdout.strip() or None
    cpu_probe = _run_probe(config.wsl_distribution, "nproc")
    memory_probe = _run_probe(config.wsl_distribution, "grep '^MemAvailable:' /proc/meminfo")
    cpu_count = int(float(cpu_probe.stdout.strip() or 1))
    memory_match = re.search(r"MemAvailable:\s*(\d+)\s*kB", memory_probe.stdout)
    available_ram = int(memory_match.group(1)) * 1024 if memory_match else 0
    ranks = config.mpi_ranks if config.mpi_ranks is not None else min(cpu_count, 8)
    expected_hashes = {
        "seeder": config.seeder_expected_sha256,
        "musubi": config.musubi_expected_sha256,
        "mus_harvesting": config.harvesting_expected_sha256,
    }
    hash_checks: dict[str, dict[str, Any]] = {}
    for label, expected in expected_hashes.items():
        executable = binaries[label]
        probe = _run_probe(
            config.wsl_distribution,
            f"sha256sum -- {shlex.quote(executable or '')}",
        )
        actual = probe.stdout.split(maxsplit=1)[0] if probe.returncode == 0 else None
        hash_checks[label] = {
            "path": executable,
            "expected_sha256": expected,
            "actual_sha256": actual,
            "status": "PASS" if actual == expected else "FAIL",
        }
        commands.append({"command": ["sha256sum", executable], "returncode": probe.returncode, "stdout": probe.stdout.strip(), "stderr": probe.stderr.strip()})
    status = (
        "PASS"
        if all(binaries.values())
        and all(item["status"] == "PASS" for item in hash_checks.values())
        else "CFD_FLOW_ENVIRONMENT_BLOCKED"
    )
    commands.append({"binary_sha256_checks": hash_checks})
    return ApesEnvironment(
        status=status,
        execution_environment="WSL2",
        wsl_distribution=config.wsl_distribution,
        binaries=binaries,
        compilers=compilers,
        mpi_ranks=ranks,
        cpu_count=cpu_count,
        available_ram_bytes=available_ram,
        commands=tuple(commands),
    )


def windows_to_wsl(path: Path, distribution: str) -> str:
    # wsl.exe otherwise treats Windows backslashes as shell escapes before
    # wslpath sees the argument.  Drive-letter paths with forward slashes are
    # accepted by wslpath and preserve every component exactly.
    windows_path = str(Path(path).resolve()).replace("\\", "/")
    probe = subprocess.run(
        ["wsl.exe", "-d", distribution, "--", "wslpath", "-a", windows_path],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    if probe.returncode != 0 or not probe.stdout.strip():
        raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", f"wslpath failed: {probe.stderr.strip()}")
    return probe.stdout.strip()


def run_wsl_tool(
    *,
    distribution: str,
    workdir: Path,
    command: list[str],
    stdout_path: Path,
    stderr_path: Path,
    timeout_s: int,
) -> ExternalRun:
    """Run one official executable and preserve separate stdout/stderr."""

    workdir_wsl = windows_to_wsl(workdir, distribution)
    command_text = " ".join(shlex.quote(part) for part in command)
    shell_text = f"cd {shlex.quote(workdir_wsl)} && {command_text}"
    try:
        started = time.perf_counter()
        process = subprocess.run(
            ["wsl.exe", "-d", distribution, "--", "bash", "-lc", shell_text],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            env=os.environ.copy(),
        )
    except subprocess.TimeoutExpired as error:
        elapsed = time.perf_counter() - started
        stdout_path.write_text(error.stdout or "", encoding="utf-8")
        stderr_path.write_text(error.stderr or "", encoding="utf-8")
        return ExternalRun(124, tuple(command), stdout_path, stderr_path, elapsed)
    elapsed = time.perf_counter() - started
    stdout_path.write_text(process.stdout, encoding="utf-8")
    stderr_path.write_text(process.stderr, encoding="utf-8")
    return ExternalRun(process.returncode, tuple(command), stdout_path, stderr_path, elapsed)


def parse_mesh_header(mesh_dir: Path) -> dict[str, Any]:
    """Parse the official tree header and per-cell boundary IDs.

    Seeder writes the total fluid-cell count to ``header.lua`` and writes one
    row of 26 signed 64-bit boundary IDs for every boundary-adjacent fluid
    cell to ``bnd.lsb`` (or ``bnd.msb``).  Counting labels in ``bnd.lua`` only
    proves that names were registered; it does not prove that any cell uses a
    label.  This parser therefore reads the documented boundary-ID property
    rather than treating label occurrences as boundary-cell counts.
    """

    headers = sorted(mesh_dir.glob("*header*.lua")) + sorted(mesh_dir.glob("*Header*.lua"))
    if not headers:
        raise FlowError("CFD_FLOW_SEEDER_FAILED", "Seeder produced no mesh header")
    text = headers[0].read_text(encoding="utf-8", errors="replace")
    fluid_match = re.search(r"(?m)^\s*nElems\s*=\s*(\d+)\s*,?\s*$", text)
    if fluid_match is None:
        fluid_match = re.search(r"(?m)^\s*nElements\s*=\s*(\d+)\s*,?\s*$", text)
    fluid_count = int(fluid_match.group(1)) if fluid_match else 0
    minimum_level = re.search(r"(?m)^\s*minLevel\s*=\s*(\d+)", text)
    maximum_level = re.search(r"(?m)^\s*maxLevel\s*=\s*(\d+)", text)

    property_blocks = re.findall(r"\{\s*label\s*=\s*['\"]([^'\"]+)['\"]\s*,(.*?)\}", text, re.DOTALL)
    property_counts: dict[str, int] = {}
    for label, block in property_blocks:
        count = re.search(r"nElems\s*=\s*(\d+)", block)
        if count:
            property_counts[label] = int(count.group(1))
    boundary_elements = property_counts.get("has boundaries", 0)

    bnd_header = mesh_dir / "bnd.lua"
    if not bnd_header.is_file():
        raise FlowError("CFD_FLOW_SEEDER_FAILED", "Seeder produced no bnd.lua boundary header")
    bnd_text = bnd_header.read_text(encoding="utf-8", errors="replace")
    sides_match = re.search(r"nSides\s*=\s*(\d+)", bnd_text)
    labels_block = re.search(r"bclabel\s*=\s*\{(.*?)\}", bnd_text, re.DOTALL)
    if sides_match is None or labels_block is None:
        raise FlowError("CFD_FLOW_SEEDER_FAILED", "Could not parse bnd.lua")
    side_count = int(sides_match.group(1))
    boundary_labels = re.findall(r"['\"]([^'\"]+)['\"]", labels_block.group(1))
    expected_labels = ("wall", "inlet", "outlet_01", "outlet_02", "outlet_03")
    if set(boundary_labels) != set(expected_labels):
        raise FlowError(
            "CFD_FLOW_SEEDER_FAILED",
            f"Unexpected boundary labels in bnd.lua: {boundary_labels}",
        )

    boundary_data = mesh_dir / "bnd.lsb"
    byte_order = "<"
    if not boundary_data.is_file():
        boundary_data = mesh_dir / "bnd.msb"
        byte_order = ">"
    if not boundary_data.is_file() or boundary_elements <= 0:
        raise FlowError("CFD_FLOW_SEEDER_FAILED", "Seeder produced no usable boundary-ID property")
    ids = np.fromfile(boundary_data, dtype=np.dtype(f"{byte_order}i8"))
    expected_values = boundary_elements * side_count
    if ids.size != expected_values:
        raise FlowError(
            "CFD_FLOW_SEEDER_FAILED",
            f"Boundary-ID size {ids.size} != {boundary_elements}*{side_count}",
        )
    ids = ids.reshape(boundary_elements, side_count)
    label_counts = {
        label: int(np.count_nonzero(np.any(ids == label_id, axis=1)))
        for label_id, label in enumerate(boundary_labels, start=1)
    }
    return {
        "header": str(headers[0]),
        "fluid_element_count": fluid_count,
        "minimum_level": int(minimum_level.group(1)) if minimum_level else None,
        "maximum_level": int(maximum_level.group(1)) if maximum_level else None,
        "property_element_counts": property_counts,
        "boundary_element_count": boundary_elements,
        "boundary_side_count": side_count,
        "boundary_labels_in_file_order": boundary_labels,
        "boundary_cell_counts": label_counts,
        # Kept as a compatibility alias for existing result readers.
        "boundary_label_occurrences": label_counts,
        "mesh_files": [str(path) for path in sorted(mesh_dir.iterdir()) if path.is_file()],
    }


def save_lua_files(
    layout: RunLayout,
    config: FlowConfig,
    partition: SurfacePartition,
    seed: SeedPoint,
    cube: BoundingCube,
    bc: BoundaryConditions,
    scaling: LatticeScaling,
) -> None:
    (layout.seeder / "seeder.lua").write_text(generate_seeder_lua(partition, seed, cube), encoding="utf-8")
    (layout.musubi / "musubi.lua").write_text(generate_musubi_lua(config, partition, bc, scaling), encoding="utf-8")
    (layout.musubi / "mus_harvester.lua").write_text(generate_harvester_lua(), encoding="utf-8")
    (layout.musubi / "restart").mkdir(exist_ok=True)


def environment_report(environment: ApesEnvironment, config: ApesConfig) -> dict[str, Any]:
    value = asdict(environment)
    value["official_provenance"] = {
        "seeder_repository": config.seeder_repository,
        "seeder_commit": config.seeder_commit,
        "musubi_repository": config.musubi_repository,
        "musubi_commit": config.musubi_commit,
        "proteus_repository": config.proteus_repository,
        "proteus_reference_commit": config.proteus_reference_commit,
        "upstream_source_modified": False,
    }
    return value


def save_environment(layout: RunLayout, environment: ApesEnvironment, config: ApesConfig) -> None:
    write_json(layout.input / "tool_versions.json", environment_report(environment, config))
