"""One-way RBC-induced drift--diffusion reduced-order transport model.

The network calculation converts one root discharge haematocrit into local
discharge and tube haematocrit values.  The particle calculation supplies the
conditional-mean transverse drift, the Ito/Fick correction for the spatially
varying transverse diffusivity, and the local diffusion coefficient used by the
geometric stochastic step.  It does not change the accepted CFD field,
viscosity, particle mobility, torque, or molecular rates.
"""

from dataclasses import dataclass
import math
import numpy as np
try:
    from numba import njit
except ImportError:  
    njit = None
from ..geometry.vessel_bed_topology import build_vessel_bed_topology


RBC_MAJOR_DIAMETER_UM                   = 8.0
RBC_MARGINATION_STRAIN                  = 450.0
CFL_PLATEAU_DIAMETER_UM                 = 8.0
CFL_PLATEAU_THICKNESS_UM                = 1.0
CFL_DIAMETER_MIN_UM                     = 3.0
CFL_DIAMETER_MAX_UM                     = 50.0
TUBE_HEMATOCRIT_QUANTITATIVE_MIN        = 0.15
TUBE_HEMATOCRIT_QUANTITATIVE_MAX        = 0.30
SHEAR_RATE_QUANTITATIVE_MIN_S_INV       = 100.0
SHEAR_RATE_QUANTITATIVE_MAX_S_INV       = 10_000.0
REFERENCE_TUBE_HEMATOCRIT               = 0.24
REFERENCE_SHEAR_RATE_S_INV              = 1_000.0
REFERENCE_TRANSVERSE_DIFFUSIVITY_UM2_S  = 45.0
RBC_RNG_ALGORITHM                       = "splitmix64_box_muller_v1"


@dataclass(frozen=True, slots=True)
class RedBloodCellNetworkState:
    """Per-vessel haematocrit state and branch audit values."""

    enabled: bool
    root_discharge_hematocrit: float
    effective_rbc_diameter_um: float
    vessel_id: np.ndarray
    parent_vessel_id: np.ndarray
    child_vessel_id: np.ndarray
    topology_depth: np.ndarray
    flow_rate_um3_s: np.ndarray
    diameter_um: np.ndarray
    cfl_width_um: np.ndarray
    scale_activation: np.ndarray
    discharge_hematocrit: np.ndarray
    tube_hematocrit: np.ndarray
    phase_flow_fraction_child_1: np.ndarray
    phase_rbc_flux_fraction_child_1: np.ndarray
    phase_coefficient_a: np.ndarray
    phase_coefficient_b: np.ndarray
    phase_cutoff_x0: np.ndarray
    rbc_flux_conservation_residual_um3_s: np.ndarray
    rbc_flux_conservation_relative_error: np.ndarray
    cfl_diameter_in_range: np.ndarray
    tube_hematocrit_in_quantitative_range: np.ndarray
    dense_diameter_um_by_vessel_id: np.ndarray
    dense_cfl_width_um_by_vessel_id: np.ndarray
    dense_scale_activation_by_vessel_id: np.ndarray
    dense_discharge_hematocrit_by_vessel_id: np.ndarray
    dense_tube_hematocrit_by_vessel_id: np.ndarray
    dense_downstream_tangent_xz_by_vessel_id: np.ndarray

    @property
    def maximum_rbc_flux_conservation_relative_error(self):
        finite = np.isfinite(self.rbc_flux_conservation_relative_error)
        return (
            float(np.max(self.rbc_flux_conservation_relative_error[finite]))
            if np.any(finite)
            else 0.0
        )

    def to_npz_payload(self):
        """Return an explicit, unit-bearing persistence payload."""

        return {
            "model_schema": np.asarray(
                ["rbc_induced_drift_diffusion_reduced_order_transport_v2_thg_cfl"]
            ),
            "enabled": np.asarray([self.enabled], dtype=bool),
            "root_discharge_hematocrit": np.asarray(
                [self.root_discharge_hematocrit], dtype=np.float64
            ),
            "effective_rbc_diameter_um": np.asarray(
                [self.effective_rbc_diameter_um], dtype=np.float64
            ),
            "rbc_major_diameter_um": np.asarray(
                [RBC_MAJOR_DIAMETER_UM], dtype=np.float64
            ),
            "cfl_model": np.asarray(
                ["thg_3pef_cortical_mouse_linear_to_plateau_v1"]
            ),
            "cfl_plateau_diameter_um": np.asarray(
                [CFL_PLATEAU_DIAMETER_UM], dtype=np.float64
            ),
            "cfl_plateau_thickness_um": np.asarray(
                [CFL_PLATEAU_THICKNESS_UM], dtype=np.float64
            ),
            "cfl_measured_diameter_range_um": np.asarray(
                [CFL_DIAMETER_MIN_UM, CFL_DIAMETER_MAX_UM], dtype=np.float64
            ),
            "reference_cumulative_shear_strain": np.asarray(
                [RBC_MARGINATION_STRAIN], dtype=np.float64
            ),
            "reference_transverse_diffusivity_um2_s": np.asarray(
                [REFERENCE_TRANSVERSE_DIFFUSIVITY_UM2_S], dtype=np.float64
            ),
            "reference_tube_hematocrit": np.asarray(
                [REFERENCE_TUBE_HEMATOCRIT], dtype=np.float64
            ),
            "reference_shear_rate_s_inv": np.asarray(
                [REFERENCE_SHEAR_RATE_S_INV], dtype=np.float64
            ),
            "vessel_id": self.vessel_id,
            "parent_vessel_id": self.parent_vessel_id,
            "child_vessel_id": self.child_vessel_id,
            "topology_depth": self.topology_depth,
            "flow_rate_um3_s": self.flow_rate_um3_s,
            "diameter_um": self.diameter_um,
            "cfl_width_um": self.cfl_width_um,
            "scale_activation": self.scale_activation,
            "discharge_hematocrit": self.discharge_hematocrit,
            "tube_hematocrit": self.tube_hematocrit,
            "phase_flow_fraction_child_1": self.phase_flow_fraction_child_1,
            "phase_rbc_flux_fraction_child_1": (
                self.phase_rbc_flux_fraction_child_1
            ),
            "phase_coefficient_a": self.phase_coefficient_a,
            "phase_coefficient_b": self.phase_coefficient_b,
            "phase_cutoff_x0": self.phase_cutoff_x0,
            "rbc_flux_conservation_residual_um3_s": (
                self.rbc_flux_conservation_residual_um3_s
            ),
            "rbc_flux_conservation_relative_error": (
                self.rbc_flux_conservation_relative_error
            ),
            "cfl_diameter_in_range": self.cfl_diameter_in_range,
            "tube_hematocrit_in_quantitative_range": (
                self.tube_hematocrit_in_quantitative_range
            ),
            "downstream_tangent_xz": (
                self.dense_downstream_tangent_xz_by_vessel_id[self.vessel_id]
            ),
        }


@dataclass(frozen=True, slots=True)
class RedBloodCellTransportEvaluation:
    """Per-particle coefficients for deterministic and stochastic RBC transport."""

    # Combined conditional-mean drift and Fick correction.
    velocity_xz_um_s: np.ndarray
    drift_velocity_xz_um_s: np.ndarray
    fick_velocity_xz_um_s: np.ndarray
    transverse_direction_xz: np.ndarray
    local_vessel_diameter_um: np.ndarray
    discharge_hematocrit: np.ndarray
    tube_hematocrit: np.ndarray
    shear_rate_s_inv: np.ndarray
    cfl_width_um: np.ndarray
    target_gap_um: np.ndarray
    center_gap_um: np.ndarray
    core_transverse_diffusivity_um2_s: np.ndarray
    transverse_diffusivity_um2_s: np.ndarray
    diffusion_taper_coordinate: np.ndarray
    diffusion_taper: np.ndarray
    margination_length_um: np.ndarray
    margination_time_s: np.ndarray
    scale_activation: np.ndarray
    nearest_wall_unique: np.ndarray
    hematocrit_in_quantitative_range: np.ndarray
    shear_rate_in_quantitative_range: np.ndarray
    quantitative_applicability: np.ndarray
    transverse_space_valid: np.ndarray
    diffusion_enabled: np.ndarray


def build_red_blood_cell_network(vessels, root_discharge_hematocrit):
    """
    Propagate root discharge haematocrit through one binary perfusion tree.
    """

    root_hd = float(root_discharge_hematocrit)
    if not math.isfinite(root_hd) or not 0.0 <= root_hd <= 1.0:
        raise ValueError("root_discharge_hematocrit must be finite and in [0, 1].")
    if not vessels:
        raise ValueError("Red-blood-cell transport requires at least one vessel.")

    ordered         = sorted(vessels, key=lambda vessel: int(vessel.vid))
    vessel_by_id    = {int(vessel.vid): vessel for vessel in ordered}
    if len(vessel_by_id) != len(ordered):
        raise ValueError("Vessel IDs must be unique for red-blood-cell transport.")
    if any(vessel_id < 0 for vessel_id in vessel_by_id):
        raise ValueError("Vessel IDs must be non-negative for categorical lookup.")

    vessel_ids  = np.asarray(sorted(vessel_by_id), dtype=np.int32)
    row_by_id   = {int(vessel_id): row for row, vessel_id in enumerate(vessel_ids)}
    count       = int(vessel_ids.size)
    diameter    = np.asarray([2.0 * float(vessel_by_id[int(vessel_id)].radius) for vessel_id in vessel_ids], dtype=np.float64)
    flow        = np.asarray([abs(float(vessel_by_id[int(vessel_id)].flow_rate)) for vessel_id in vessel_ids], dtype=np.float64)
    downstream_tangent = np.asarray(
        [
            np.asarray(vessel_by_id[int(vessel_id)].x_d, dtype=np.float64)[[0, 2]]
            - np.asarray(vessel_by_id[int(vessel_id)].x_p, dtype=np.float64)[[0, 2]]
            for vessel_id in vessel_ids
        ],
        dtype=np.float64,
    )

    tangent_norm = np.linalg.norm(downstream_tangent, axis=1)
    if np.any(~np.isfinite(diameter)) or np.any(diameter <= 0.0):
        raise ValueError("Every vessel diameter must be finite and positive.")
    if np.any(~np.isfinite(flow)):
        raise ValueError("Every vessel flow rate must be finite.")
    if np.any(~np.isfinite(downstream_tangent)) or np.any(tangent_norm <= 0.0):
        raise ValueError("Every vessel must have a finite non-zero downstream X-Z direction.")
    downstream_tangent /= tangent_norm[:, None]

    diameter_valid          = (diameter >= CFL_DIAMETER_MIN_UM) & (diameter <= CFL_DIAMETER_MAX_UM)
    vessel_cfl_width        = cfl_width_um(diameter)
    vessel_scale_activation = scale_activation(diameter)
    if root_hd == 0.0:
        return _disabled_network_state(
            vessel_ids,
            diameter,
            flow,
            diameter_valid,
            downstream_tangent,
            vessel_cfl_width,
            vessel_scale_activation,
        )
    topology = build_vessel_bed_topology(ordered)
    if len(topology.root_unit_ids) != 1:
        raise ValueError(
            "RBC transport requires exactly one perfused root."
        )

    parent = np.full(count, -1, dtype=np.int32)
    children = np.full((count, 2), -1, dtype=np.int32)
    depth = np.zeros(count, dtype=np.int32)
    discharge = np.full(count, np.nan, dtype=np.float64)
    flow_fraction = np.full(count, np.nan, dtype=np.float64)
    rbc_fraction = np.full(count, np.nan, dtype=np.float64)
    coefficient_a = np.full(count, np.nan, dtype=np.float64)
    coefficient_b = np.full(count, np.nan, dtype=np.float64)
    cutoff_x0 = np.full(count, np.nan, dtype=np.float64)
    flux_residual = np.full(count, np.nan, dtype=np.float64)
    flux_relative_error = np.full(count, np.nan, dtype=np.float64)

    units = topology.units
    root_unit = int(topology.root_unit_ids[0])
    unit_hd = np.full(len(units), np.nan, dtype=np.float64)
    unit_hd[root_unit] = root_hd
    unit_order = sorted(range(len(units)), key=lambda index: units[index].topology_depth)

    for unit_id in unit_order:
        unit = units[unit_id]
        hd_parent = float(unit_hd[unit_id])
        if not math.isfinite(hd_parent):
            raise ValueError("Red-blood-cell topology traversal reached an unassigned unit.")
        chain = tuple(int(value) for value in unit.segment_ids)
        if not chain:
            raise ValueError("Red-blood-cell topology contains an empty vessel unit.")
        unit_flow = float(unit.flow_rate_um3_s)
        if unit_flow <= 0.0:
            raise ValueError("RBC phase separation requires positive branch flow.")
        for chain_index, vessel_id in enumerate(chain):
            row = row_by_id[vessel_id]
            discharge[row] = hd_parent
            depth[row] = int(unit.topology_depth)
            flow_scale = max(
                abs(flow[row]), abs(unit_flow), np.finfo(np.float64).tiny
            )
            if abs(flow[row] - unit_flow) > 1.0e-10 * flow_scale:
                raise ValueError(
                    f"Vessel {vessel_id} changes flow inside one unbranched topology unit."
                )
            if chain_index > 0:
                parent[row] = chain[chain_index - 1]
                children[row_by_id[chain[chain_index - 1]], 0] = vessel_id

        first_vessel = chain[0]
        if unit.parent_unit_id >= 0:
            parent_segment = int(units[unit.parent_unit_id].segment_ids[-1])
            parent[row_by_id[first_vessel]] = parent_segment

        child_units = tuple(int(value) for value in unit.child_unit_ids)
        if len(child_units) not in {0, 2}:
            raise ValueError(
                "RBC phase separation requires terminal or binary topology units."
            )
        if not child_units:
            continue

        parent_segment = chain[-1]
        child_segments = tuple(int(units[value].segment_ids[0]) for value in child_units)
        children[row_by_id[parent_segment], :] = child_segments
        q_parent = unit_flow
        q_1 = float(units[child_units[0]].flow_rate_um3_s)
        q_2 = float(units[child_units[1]].flow_rate_um3_s)
        if min(q_parent, q_1, q_2) <= 0.0:
            raise ValueError(
                f"Branch at vessel {parent_segment} requires strictly positive flows."
            )
        relative_flow_error = abs(q_parent - q_1 - q_2) / q_parent
        if relative_flow_error > 1.0e-8:
            raise ValueError(
                "Flow conservation is inconsistent with RBC phase separation at "
                f"parent vessel {parent_segment}: "
                f"relative_error={relative_flow_error:.17g}."
            )
        d_parent = diameter[row_by_id[parent_segment]]
        d_1 = diameter[row_by_id[child_segments[0]]]
        d_2 = diameter[row_by_id[child_segments[1]]]
        fq, fe, a_value, b_value, x0_value = phase_separation_fraction(
            q_parent, q_1, hd_parent, d_parent, d_1, d_2
        )
        hd_1 = hd_parent * fe / fq
        hd_2 = hd_parent * (1.0 - fe) / (1.0 - fq)
        for child_hd, child_unit in ((hd_1, child_units[0]), (hd_2, child_units[1])):
            if not math.isfinite(child_hd) or not 0.0 <= child_hd <= 1.0:
                raise ValueError(
                    "RBC phase separation produced an invalid discharge "
                    f"haematocrit at parent vessel {parent_segment}: {child_hd}."
                )
            unit_hd[child_unit] = child_hd

        residual = q_parent * hd_parent - q_1 * hd_1 - q_2 * hd_2
        scale = max(abs(q_parent * hd_parent), np.finfo(np.float64).tiny)
        relative = abs(residual) / scale
        if relative > 1.0e-10:
            raise RuntimeError(
                "Red-blood-cell flux conservation failed at parent vessel "
                f"{parent_segment}: relative_error={relative:.17g}."
            )
        row = row_by_id[parent_segment]
        flow_fraction[row] = fq
        rbc_fraction[row] = fe
        coefficient_a[row] = a_value
        coefficient_b[row] = b_value
        cutoff_x0[row] = x0_value
        flux_residual[row] = residual
        flux_relative_error[row] = relative

    if np.any(~np.isfinite(discharge)):
        missing = vessel_ids[~np.isfinite(discharge)].tolist()
        raise ValueError(f"Unassigned discharge haematocrit for vessel_ids={missing}.")
    tube = fahraeus_tube_hematocrit(discharge, diameter)
    if np.any(~np.isfinite(tube)) or np.any((tube < 0.0) | (tube > 1.0)):
        raise ValueError("The Fåhræus relation produced an invalid tube haematocrit.")

    max_vessel_id = int(vessel_ids[-1])
    dense_diameter = np.full(max_vessel_id + 1, np.nan, dtype=np.float64)
    dense_cfl_width = np.full(max_vessel_id + 1, np.nan, dtype=np.float64)
    dense_scale_activation = np.full(max_vessel_id + 1, np.nan, dtype=np.float64)
    dense_discharge = np.full(max_vessel_id + 1, np.nan, dtype=np.float64)
    dense_tube = np.full(max_vessel_id + 1, np.nan, dtype=np.float64)
    dense_tangent = np.full((max_vessel_id + 1, 2), np.nan, dtype=np.float64)
    dense_diameter[vessel_ids] = diameter
    dense_cfl_width[vessel_ids] = vessel_cfl_width
    dense_scale_activation[vessel_ids] = vessel_scale_activation
    dense_discharge[vessel_ids] = discharge
    dense_tube[vessel_ids] = tube
    dense_tangent[vessel_ids] = downstream_tangent
    tube_valid = (tube >= TUBE_HEMATOCRIT_QUANTITATIVE_MIN) & (
        tube <= TUBE_HEMATOCRIT_QUANTITATIVE_MAX
    )

    return RedBloodCellNetworkState(
        enabled=True,
        root_discharge_hematocrit=root_hd,
        effective_rbc_diameter_um=RBC_MAJOR_DIAMETER_UM,
        vessel_id=vessel_ids,
        parent_vessel_id=parent,
        child_vessel_id=children,
        topology_depth=depth,
        flow_rate_um3_s=flow,
        diameter_um=diameter,
        cfl_width_um=vessel_cfl_width,
        scale_activation=vessel_scale_activation,
        discharge_hematocrit=discharge,
        tube_hematocrit=tube,
        phase_flow_fraction_child_1=flow_fraction,
        phase_rbc_flux_fraction_child_1=rbc_fraction,
        phase_coefficient_a=coefficient_a,
        phase_coefficient_b=coefficient_b,
        phase_cutoff_x0=cutoff_x0,
        rbc_flux_conservation_residual_um3_s=flux_residual,
        rbc_flux_conservation_relative_error=flux_relative_error,
        cfl_diameter_in_range=diameter_valid,
        tube_hematocrit_in_quantitative_range=tube_valid,
        dense_diameter_um_by_vessel_id=dense_diameter,
        dense_cfl_width_um_by_vessel_id=dense_cfl_width,
        dense_scale_activation_by_vessel_id=dense_scale_activation,
        dense_discharge_hematocrit_by_vessel_id=dense_discharge,
        dense_tube_hematocrit_by_vessel_id=dense_tube,
        dense_downstream_tangent_xz_by_vessel_id=dense_tangent,
    )


def phase_separation_fraction(
    parent_flow_um3_s,
    child_1_flow_um3_s,
    parent_discharge_hematocrit,
    parent_diameter_um,
    child_1_diameter_um,
    child_2_diameter_um,
):
    """Evaluate the Pries et al. 1990 binary phase-separation relation."""

    q_parent = float(parent_flow_um3_s)
    q_1 = float(child_1_flow_um3_s)
    hd = float(parent_discharge_hematocrit)
    d_parent = float(parent_diameter_um)
    d_1 = float(child_1_diameter_um)
    d_2 = float(child_2_diameter_um)
    if not all(math.isfinite(value) for value in (q_parent, q_1, hd, d_parent, d_1, d_2)):
        raise ValueError("Phase-separation inputs must be finite.")
    if q_parent <= 0.0 or not 0.0 < q_1 < q_parent:
        raise ValueError("A binary phase-separation split requires 0 < Q1 < Qp.")
    if not 0.0 <= hd <= 1.0 or min(d_parent, d_1, d_2) <= 0.0:
        raise ValueError("Phase-separation haematocrit and diameters are invalid.")

    fq = q_1 / q_parent
    a_value = -6.96 * math.log(d_1 / d_2) / d_parent
    b_value = 1.0 + 6.98 * (1.0 - hd) / d_parent
    x0_value = 0.4 / d_parent
    if x0_value >= 0.5:
        raise ValueError("The phase-separation cutoff leaves no interior flow interval.")
    if fq <= x0_value:
        fe = 0.0
    elif fq >= 1.0 - x0_value:
        fe = 1.0
    else:
        z_value = (fq - x0_value) / (1.0 - 2.0 * x0_value)
        log_odds = a_value + b_value * math.log(z_value / (1.0 - z_value))
        if log_odds >= 0.0:
            factor = math.exp(-log_odds)
            fe = 1.0 / (1.0 + factor)
        else:
            factor = math.exp(log_odds)
            fe = factor / (1.0 + factor)
    return fq, fe, a_value, b_value, x0_value


def fahraeus_tube_hematocrit(discharge_hematocrit, diameter_um):
    """Convert discharge to tube haematocrit with the retained relation."""

    hd = np.asarray(discharge_hematocrit, dtype=np.float64)
    diameter = np.asarray(diameter_um, dtype=np.float64)
    ratio = hd + (1.0 - hd) * (
        1.0 + 1.7 * np.exp(-0.415 * diameter) - 0.6 * np.exp(-0.011 * diameter)
    )
    return hd * ratio


def cfl_width_um(diameter_um):
    """Approximate the cortical-vessel CFL trend measured with THG and 3PEF."""

    diameter = np.asarray(diameter_um, dtype=np.float64)
    if np.any(~np.isfinite(diameter)) or np.any(diameter <= 0.0):
        raise ValueError("CFL vessel diameter must be finite and positive.")
    return CFL_PLATEAU_THICKNESS_UM * np.minimum(
        diameter / CFL_PLATEAU_DIAMETER_UM,
        1.0,
    )


def scale_activation(diameter_um):
    """Return the smooth vessel-scale switch using the 8 um RBC major diameter."""

    diameter = np.asarray(diameter_um, dtype=np.float64)
    chi = diameter / RBC_MAJOR_DIAMETER_UM
    scaled = (chi - 1.25) / 0.25
    interior = np.clip(scaled, 0.0, 1.0)
    smooth = 3.0 * interior * interior - 2.0 * interior * interior * interior
    return np.where(chi <= 1.25, 0.0, np.where(chi >= 1.5, 1.0, smooth))


def evaluate_red_blood_cell_transport(
    sampled_trajectory_vessel_id,
    velocity_gradient_s_inv,
    wall_gap_um,
    inward_wall_normal_xz,
    nearest_wall_unique,
    bubble_radius_um,
    active,
    network,
    *,
    use_numba,
):
    """Evaluate RBC drift, Fick correction, and diffusivity for one snapshot."""

    vessel_id = np.ascontiguousarray(sampled_trajectory_vessel_id, dtype=np.int32)
    gradients = np.ascontiguousarray(velocity_gradient_s_inv, dtype=np.float64)
    gap = np.ascontiguousarray(wall_gap_um, dtype=np.float64)
    normal = np.ascontiguousarray(inward_wall_normal_xz, dtype=np.float64)
    unique = np.ascontiguousarray(nearest_wall_unique, dtype=np.bool_)
    radius = np.ascontiguousarray(bubble_radius_um, dtype=np.float64)
    active_mask = np.ascontiguousarray(active, dtype=np.bool_)
    count = int(vessel_id.size)
    if (
        gradients.shape != (count, 2, 2)
        or gap.shape != (count,)
        or normal.shape != (count, 2)
        or unique.shape != (count,)
        or radius.shape != (count,)
        or active_mask.shape != (count,)
    ):
        raise ValueError("Red-blood-cell particle arrays have inconsistent shapes.")
    if not network.enabled:
        raise ValueError("A disabled red-blood-cell network must not enter the RHS.")

    kernel = (
        _evaluate_red_blood_cell_transport_numba
        if use_numba and njit is not None
        else _evaluate_red_blood_cell_transport_kernel
    )
    values = kernel(
        vessel_id,
        gradients,
        gap,
        normal,
        unique,
        radius,
        active_mask,
        network.dense_diameter_um_by_vessel_id,
        network.dense_cfl_width_um_by_vessel_id,
        network.dense_scale_activation_by_vessel_id,
        network.dense_discharge_hematocrit_by_vessel_id,
        network.dense_tube_hematocrit_by_vessel_id,
        network.dense_downstream_tangent_xz_by_vessel_id,
    )
    return RedBloodCellTransportEvaluation(*values)


def _evaluate_red_blood_cell_transport_kernel(
    trajectory_vessel_id,
    gradients,
    gap,
    normal,
    unique,
    radius,
    active,
    diameter_by_id,
    cfl_width_by_id,
    scale_activation_by_id,
    discharge_by_id,
    tube_by_id,
    downstream_tangent_by_id,
):
    count = trajectory_vessel_id.size
    velocity = np.zeros((count, 2), dtype=np.float64)
    drift_velocity = np.zeros((count, 2), dtype=np.float64)
    fick_velocity = np.zeros((count, 2), dtype=np.float64)
    transverse_direction = np.zeros((count, 2), dtype=np.float64)
    diameter = np.full(count, np.nan, dtype=np.float64)
    discharge = np.full(count, np.nan, dtype=np.float64)
    tube = np.full(count, np.nan, dtype=np.float64)
    shear = np.full(count, np.nan, dtype=np.float64)
    cfl = np.full(count, np.nan, dtype=np.float64)
    target_gap = np.full(count, np.nan, dtype=np.float64)
    center_gap = np.full(count, np.nan, dtype=np.float64)
    core_diffusivity = np.full(count, np.nan, dtype=np.float64)
    diffusivity = np.full(count, np.nan, dtype=np.float64)
    taper_coordinate = np.full(count, np.nan, dtype=np.float64)
    taper = np.full(count, np.nan, dtype=np.float64)
    length = np.full(count, np.nan, dtype=np.float64)
    time_scale = np.full(count, np.inf, dtype=np.float64)
    activation = np.zeros(count, dtype=np.float64)
    wall_unique = np.zeros(count, dtype=np.bool_)
    hematocrit_valid = np.zeros(count, dtype=np.bool_)
    shear_valid = np.zeros(count, dtype=np.bool_)
    quantitative_valid = np.zeros(count, dtype=np.bool_)
    transverse_space_valid = np.zeros(count, dtype=np.bool_)
    diffusion_enabled = np.zeros(count, dtype=np.bool_)

    for lane in range(count):
        if not active[lane]:
            continue
        actual_vessel_id = int(trajectory_vessel_id[lane]) - 1
        if actual_vessel_id < 0 or actual_vessel_id >= diameter_by_id.size:
            raise ValueError("An active particle has no categorical vessel owner.")
        local_diameter = diameter_by_id[actual_vessel_id]
        local_cfl = cfl_width_by_id[actual_vessel_id]
        local_activation = scale_activation_by_id[actual_vessel_id]
        local_discharge = discharge_by_id[actual_vessel_id]
        local_tube = tube_by_id[actual_vessel_id]
        tangent_x = downstream_tangent_by_id[actual_vessel_id, 0]
        tangent_z = downstream_tangent_by_id[actual_vessel_id, 1]
        if not (
            math.isfinite(local_diameter)
            and math.isfinite(local_cfl)
            and math.isfinite(local_activation)
            and math.isfinite(local_discharge)
            and math.isfinite(local_tube)
            and math.isfinite(tangent_x)
            and math.isfinite(tangent_z)
        ):
            raise ValueError("An active particle references an unknown vessel ID.")
        transverse_x = -tangent_z
        transverse_z = tangent_x
        transverse_direction[lane, 0] = transverse_x
        transverse_direction[lane, 1] = transverse_z

        diameter[lane] = local_diameter
        discharge[lane] = local_discharge
        tube[lane] = local_tube
        e00 = gradients[lane, 0, 0]
        e11 = gradients[lane, 1, 1]
        e01 = 0.5 * (gradients[lane, 0, 1] + gradients[lane, 1, 0])
        local_shear = math.sqrt(max(2.0 * (e00 * e00 + e11 * e11 + 2.0 * e01 * e01), 0.0))
        shear[lane] = local_shear
        cfl[lane] = local_cfl
        local_target_gap = max(local_cfl - 2.0 * radius[lane], 0.0)
        target_gap[lane] = local_target_gap
        local_center_gap = 0.5 * local_diameter - radius[lane]
        center_gap[lane] = local_center_gap

        activation[lane] = local_activation

        local_core_diffusivity = (
            REFERENCE_TRANSVERSE_DIFFUSIVITY_UM2_S
            * local_tube
            / REFERENCE_TUBE_HEMATOCRIT
            * local_shear
            / REFERENCE_SHEAR_RATE_S_INV
            * local_activation
        )
        core_diffusivity[lane] = local_core_diffusivity
        local_length = local_center_gap - local_target_gap
        length[lane] = local_length
        local_space_valid = (
            math.isfinite(local_center_gap)
            and math.isfinite(local_length)
            and local_center_gap > 0.0
            and local_length > 0.0
        )
        transverse_space_valid[lane] = local_space_valid
        wall_unique[lane] = unique[lane]
        hematocrit_valid[lane] = (
            TUBE_HEMATOCRIT_QUANTITATIVE_MIN
            <= local_tube
            <= TUBE_HEMATOCRIT_QUANTITATIVE_MAX
        )
        shear_valid[lane] = (
            SHEAR_RATE_QUANTITATIVE_MIN_S_INV
            <= local_shear
            <= SHEAR_RATE_QUANTITATIVE_MAX_S_INV
        )
        diameter_valid = (
            CFL_DIAMETER_MIN_UM <= local_diameter <= CFL_DIAMETER_MAX_UM
        )
        quantitative_valid[lane] = (
            diameter_valid and hematocrit_valid[lane] and shear_valid[lane]
        )

        if local_activation > 0.0 and local_shear > 0.0:
            time_scale[lane] = RBC_MARGINATION_STRAIN / (
                local_activation * local_shear
            )

        if not local_space_valid:
            diffusivity[lane] = 0.0
            taper_coordinate[lane] = 0.0
            taper[lane] = 0.0
            continue

        raw_coordinate = (gap[lane] - local_target_gap) / local_length
        local_coordinate = min(max(raw_coordinate, 0.0), 1.0)
        taper_coordinate[lane] = local_coordinate
        local_taper = (
            3.0 * local_coordinate * local_coordinate
            - 2.0 * local_coordinate * local_coordinate * local_coordinate
        )
        taper[lane] = local_taper
        local_diffusivity = (
            local_core_diffusivity * local_taper
            if gap[lane] > local_target_gap and local_activation > 0.0
            else 0.0
        )
        diffusivity[lane] = local_diffusivity
        diffusion_enabled[lane] = local_diffusivity > 0.0

        excess_gap = max(gap[lane] - local_target_gap, 0.0)
        if excess_gap > 0.0 and local_activation > 0.0 and unique[lane]:
            drift_speed = (
                local_activation
                * local_shear
                / RBC_MARGINATION_STRAIN
                * excess_gap
            )
            drift_velocity[lane, 0] = -drift_speed * normal[lane, 0]
            drift_velocity[lane, 1] = -drift_speed * normal[lane, 1]

        if (
            unique[lane]
            and local_core_diffusivity > 0.0
            and local_coordinate > 0.0
            and local_coordinate < 1.0
        ):
            derivative = (
                local_core_diffusivity
                * 6.0
                * local_coordinate
                * (1.0 - local_coordinate)
                / local_length
            )
            projected_normal = (
                normal[lane, 0] * transverse_x
                + normal[lane, 1] * transverse_z
            )
            fick_speed = derivative * projected_normal
            fick_velocity[lane, 0] = fick_speed * transverse_x
            fick_velocity[lane, 1] = fick_speed * transverse_z

        velocity[lane, 0] = (
            drift_velocity[lane, 0] + fick_velocity[lane, 0]
        )
        velocity[lane, 1] = (
            drift_velocity[lane, 1] + fick_velocity[lane, 1]
        )

    return (
        velocity,
        drift_velocity,
        fick_velocity,
        transverse_direction,
        diameter,
        discharge,
        tube,
        shear,
        cfl,
        target_gap,
        center_gap,
        core_diffusivity,
        diffusivity,
        taper_coordinate,
        taper,
        length,
        time_scale,
        activation,
        wall_unique,
        hematocrit_valid,
        shear_valid,
        quantitative_valid,
        transverse_space_valid,
        diffusion_enabled,
    )


def _disabled_network_state(
    vessel_ids,
    diameter_um,
    flow_rate_um3_s,
    diameter_valid,
    downstream_tangent_xz,
    vessel_cfl_width_um,
    vessel_scale_activation,
):
    """Build inert metadata without imposing new topology rules on the baseline."""

    count = int(vessel_ids.size)
    maximum_id = int(vessel_ids[-1])
    dense_diameter = np.full(maximum_id + 1, np.nan, dtype=np.float64)
    dense_cfl_width = np.full(maximum_id + 1, np.nan, dtype=np.float64)
    dense_scale_activation = np.full(maximum_id + 1, np.nan, dtype=np.float64)
    dense_discharge = np.full(maximum_id + 1, np.nan, dtype=np.float64)
    dense_tube = np.full(maximum_id + 1, np.nan, dtype=np.float64)
    dense_tangent = np.full((maximum_id + 1, 2), np.nan, dtype=np.float64)
    dense_diameter[vessel_ids] = diameter_um
    dense_cfl_width[vessel_ids] = vessel_cfl_width_um
    dense_scale_activation[vessel_ids] = vessel_scale_activation
    dense_discharge[vessel_ids] = 0.0
    dense_tube[vessel_ids] = 0.0
    dense_tangent[vessel_ids] = downstream_tangent_xz
    nan_values = np.full(count, np.nan, dtype=np.float64)
    return RedBloodCellNetworkState(
        enabled=False,
        root_discharge_hematocrit=0.0,
        effective_rbc_diameter_um=RBC_MAJOR_DIAMETER_UM,
        vessel_id=vessel_ids,
        parent_vessel_id=np.full(count, -1, dtype=np.int32),
        child_vessel_id=np.full((count, 2), -1, dtype=np.int32),
        topology_depth=np.zeros(count, dtype=np.int32),
        flow_rate_um3_s=flow_rate_um3_s,
        diameter_um=diameter_um,
        cfl_width_um=vessel_cfl_width_um,
        scale_activation=vessel_scale_activation,
        discharge_hematocrit=np.zeros(count, dtype=np.float64),
        tube_hematocrit=np.zeros(count, dtype=np.float64),
        phase_flow_fraction_child_1=nan_values.copy(),
        phase_rbc_flux_fraction_child_1=nan_values.copy(),
        phase_coefficient_a=nan_values.copy(),
        phase_coefficient_b=nan_values.copy(),
        phase_cutoff_x0=nan_values.copy(),
        rbc_flux_conservation_residual_um3_s=nan_values.copy(),
        rbc_flux_conservation_relative_error=nan_values.copy(),
        cfl_diameter_in_range=diameter_valid,
        tube_hematocrit_in_quantitative_range=np.zeros(count, dtype=bool),
        dense_diameter_um_by_vessel_id=dense_diameter,
        dense_cfl_width_um_by_vessel_id=dense_cfl_width,
        dense_scale_activation_by_vessel_id=dense_scale_activation,
        dense_discharge_hematocrit_by_vessel_id=dense_discharge,
        dense_tube_hematocrit_by_vessel_id=dense_tube,
        dense_downstream_tangent_xz_by_vessel_id=dense_tangent,
    )
if njit is not None:
    _evaluate_red_blood_cell_transport_numba = njit(cache=True, nogil=True)(
        _evaluate_red_blood_cell_transport_kernel
    )
else:  # pragma: no cover
    _evaluate_red_blood_cell_transport_numba = _evaluate_red_blood_cell_transport_kernel
