import math
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from ..generation.stage_specs import StageConfig, parse_tinit_spec
try:
    import yaml
except ImportError as exc:  
    raise ImportError(
        "PyYAML is required to read DCCO YAML configuration files. "
        "Install it with `pip install pyyaml` in the Python environment used "
        "to run vessel_generation.py."
    ) from exc


PACKAGE_ROOT        = Path(__file__).resolve().parents[2]
CONFIG_DIR          = PACKAGE_ROOT / "configs"
DEFAULT_CONFIG_PATH = CONFIG_DIR / "dcco_config.yaml"
SWC_OUTPUT_DIR      = PACKAGE_ROOT / "vessel_swc_models"
TOP_LEVEL_KEYS      = {"output", "config", "stages"}


@dataclass(frozen=True)
class RegionMask:
    """
    Represents a region mask in the x-z plane, which can be either a circle or a rectangle. 
    The mask is used to define avascular or carriage regions in the vascular model.
    """
    kind: str
    x_um: float = 0.0
    y_um: float = 0.0
    z_um: float = 0.0
    radius_um: float = 0.0
    x_min_um: float = 0.0
    x_max_um: float = 0.0
    z_min_um: float = 0.0
    z_max_um: float = 0.0
    y_min_um: float = 0.0
    y_max_um: float = 0.0

    def contains(self, p):
        """
        Check if a point `p` (given as a tuple or list of coordinates) is within the region defined by this mask. 
        The point is expected to be in the form (x, y, z), but only the x and z coordinates are used for the check.
        """
        x = float(p[0])
        y = float(p[1])
        z = float(p[2])

        if self.kind == "circle":
            return (
                (x - self.x_um) ** 2 + (z - self.z_um) ** 2
                <= self.radius_um ** 2
            )
        if self.kind == "sphere":
            return (
                (x - self.x_um) ** 2
                + (y - self.y_um) ** 2
                + (z - self.z_um) ** 2
                <= self.radius_um ** 2
            )
        if self.kind == "rect":
            return (
                self.x_min_um <= x <= self.x_max_um
                and self.z_min_um <= z <= self.z_max_um
            )
        if self.kind == "box":
            return (
                self.x_min_um <= x <= self.x_max_um
                and self.y_min_um <= y <= self.y_max_um
                and self.z_min_um <= z <= self.z_max_um
            )
        
        raise ValueError(f"Unsupported region mask kind: {self.kind!r}")

    def signed_distance(self, p):
        """
        Return positive distance outside the mask and negative distance inside.

        This lets the finite-radius tube check keep the vessel wall, rather
        than only its centreline, outside avascular and carriage regions.
        """
        point = np.asarray(p, dtype=float)
        if self.kind == "circle":
            return float(
                np.linalg.norm(
                    point[[0, 2]] - np.asarray([self.x_um, self.z_um])
                )
                - self.radius_um
            )
        if self.kind == "sphere":
            return float(
                np.linalg.norm(
                    point - np.asarray([self.x_um, self.y_um, self.z_um])
                )
                - self.radius_um
            )
        if self.kind == "rect":
            center = np.asarray(
                [
                    0.5 * (self.x_min_um + self.x_max_um),
                    0.5 * (self.z_min_um + self.z_max_um),
                ]
            )
            half = np.asarray(
                [
                    0.5 * (self.x_max_um - self.x_min_um),
                    0.5 * (self.z_max_um - self.z_min_um),
                ]
            )
            q = np.abs(point[[0, 2]] - center) - half
        elif self.kind == "box":
            center = np.asarray(
                [
                    0.5 * (self.x_min_um + self.x_max_um),
                    0.5 * (self.y_min_um + self.y_max_um),
                    0.5 * (self.z_min_um + self.z_max_um),
                ]
            )
            half = np.asarray(
                [
                    0.5 * (self.x_max_um - self.x_min_um),
                    0.5 * (self.y_max_um - self.y_min_um),
                    0.5 * (self.z_max_um - self.z_min_um),
                ]
            )
            q = np.abs(point - center) - half
        else:
            raise ValueError(f"Unsupported region mask kind: {self.kind!r}")
        return float(np.linalg.norm(np.maximum(q, 0.0)) + min(np.max(q), 0.0))


def _read_yaml(path):
    """
    read yaml file and return a dictionary.
    """
    data = yaml.safe_load(path.read_text(encoding="utf-8-sig"))
    if data is None:
        data = {}
    if not isinstance(data, dict):
        raise ValueError(f"The YAML root in {path} must be a mapping.")

    return data


def _default_yaml_config():
    """
    read the default config mapping under 'config:' in the default YAML.
    """
    data = _read_yaml(DEFAULT_CONFIG_PATH)
    raw_config = data.get("config", {})
    if not isinstance(raw_config, dict):
        raise ValueError(f"{DEFAULT_CONFIG_PATH} must contain a mapping under 'config'.")

    return dict(raw_config)


# ============================
# only keys present in the default YAML are considered valid; any unknown keys will raise an error.
# ============================
CONFIG_KEYS = set(_default_yaml_config())
CUBE_UM     = float(_default_yaml_config()["cube_um"])
SEED        = int(_default_yaml_config()["seed"])


def _parse_region_mask(item):
    """
    Parse a single region mask entry from the YAML configuration. 
    The entry can define either a circular or rectangular region in the x-z plane. 
    The function returns a RegionMask object representing the specified region
    """
    kind = str(item.get("kind", "")).lower()
    if kind == "circle":
        return RegionMask(
            kind="circle",
            x_um=float(item["x_um"]),
            z_um=float(item["z_um"]),
            radius_um=float(item["radius_um"]),
        )
    if kind == "rect":
        return RegionMask(
            kind="rect",
            x_min_um=float(item["x_min_um"]),
            x_max_um=float(item["x_max_um"]),
            z_min_um=float(item["z_min_um"]),
            z_max_um=float(item["z_max_um"]),
        )
    if kind == "sphere":
        return RegionMask(
            kind="sphere",
            x_um=float(item["x_um"]),
            y_um=float(item["y_um"]),
            z_um=float(item["z_um"]),
            radius_um=float(item["radius_um"]),
        )
    if kind == "box":
        return RegionMask(
            kind="box",
            x_min_um=float(item["x_min_um"]),
            x_max_um=float(item["x_max_um"]),
            y_min_um=float(item["y_min_um"]),
            y_max_um=float(item["y_max_um"]),
            z_min_um=float(item["z_min_um"]),
            z_max_um=float(item["z_max_um"]),
        )
    raise ValueError(f"Unsupported mask entry: {item!r}")


def _as_region_masks(raw, key):
    """
    Convert a YAML list of region mask entries into a tuple of RegionMask objects.
    Each entry in the list should be a mapping that defines either a circular or rectangular region.
    """
    if raw is None:
        return ()
    if isinstance(raw, tuple) and all(isinstance(item, RegionMask) for item in raw):
        return raw
    if not isinstance(raw, list):
        raise ValueError(f"{key} must be a YAML list of region mask objects.")

    return tuple(_parse_region_mask(item) for item in raw)


def _normalise_config_mapping(raw):
    """
    Normalise a raw configuration mapping read from YAML.
    This function checks for unknown keys, fills in missing keys with default values, and converts region mask entries into RegionMask objects. 
    It returns a dictionary with all expected configuration keys.
    """
    unknown = sorted(set(raw) - CONFIG_KEYS)
    if unknown:
        raise ValueError(f"Unknown DCCO config key(s): {', '.join(unknown)}")

    data = dict(raw)
    for key in (
        "root_x_um",
        "root_y_um",
        "root_z_um",
        "main_trunk_position_um",
        "main_trunk_inlet_um",
        "main_trunk_outlet_um",
        "main_trunk_min_endpoint_distance_um",
        "main_trunk_phase_rad",
    ):
        if data.get(key) is None:
            data[key] = math.nan

    for key in ("avascular_regions", "carriage_regions"):
        if key in data:
            data[key] = _as_region_masks(data[key], key)

    return data


def _default_config_values():
    return _normalise_config_mapping(_default_yaml_config())


class DCCOConfig:
    """
    Represents the configuration for the DCCO vascular generator.
    The configuration is read from a YAML file and contains various parameters that control the generation of the vascular tree. 
    The configuration is immutable and can be overridden with new values to create modified configurations."""
    __slots__ = tuple(sorted(CONFIG_KEYS))

    def __init__(self, **values):
        data = _default_config_values()
        data.update(_normalise_config_mapping(values))
        for key in self.__slots__:
            setattr(self, key, data[key])
        self._validate()

    def as_dict(self):
        return {key: getattr(self, key) for key in self.__slots__}

    def with_overrides(self, **values):
        data = self.as_dict()
        data.update(_normalise_config_mapping(values))
        return DCCOConfig._from_normalised(data)

    @classmethod
    def _from_normalised(cls, data):
        cfg = cls.__new__(cls)
        missing = sorted(CONFIG_KEYS - set(data))
        if missing:
            raise ValueError(f"Missing DCCO config key(s): {', '.join(missing)}")

        for key in cfg.__slots__:
            setattr(cfg, key, data[key])
        cfg._validate()
        return cfg

    def _validate(self):
        """Reject internally inconsistent direct-flow configurations."""
        if self.geometry_mode not in {"planar_2d", "volumetric_3d"}:
            raise ValueError(
                "geometry_mode must be 'volumetric_3d' or 'planar_2d'."
            )
        if not math.isfinite(float(self.cube_um)) or float(self.cube_um) <= 0.0:
            raise ValueError("cube_um must be finite and positive.")
        if not math.isfinite(float(self.plane_y_um)):
            raise ValueError("plane_y_um must be finite.")
        for key in ("domain_size_x_um", "domain_size_y_um", "domain_size_z_um"):
            value = float(getattr(self, key))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{key} must be finite and positive.")
        if int(self.n_terminals) < 1:
            raise ValueError("n_terminals must be at least 1.")
        if self.geometry_mode == "planar_2d":
            self._validate_planar_hemodynamics()
        else:
            self._validate_volumetric_hemodynamics()
        if float(self.min_vessel_clearance_um) < 0.0:
            raise ValueError("min_vessel_clearance_um cannot be negative.")
        if float(self.boundary_surface_clearance_um) < 0.0:
            raise ValueError("boundary_surface_clearance_um cannot be negative.")
        minimum_segment_length = float(self.min_segment_length_um)
        if (
            not math.isfinite(minimum_segment_length)
            or minimum_segment_length <= 0.0
        ):
            raise ValueError("min_segment_length_um must be finite and positive.")
        minimum_opening = float(self.min_daughter_angle_deg)
        maximum_opening = float(self.max_opening_angle_deg)
        if not 0.0 <= minimum_opening < maximum_opening <= 180.0:
            raise ValueError(
                "Daughter opening angles must satisfy "
                "0 <= min_daughter_angle_deg < max_opening_angle_deg <= 180."
            )
        for key in ("sprout_cv", "sprout_cp", "sprout_cd"):
            value = float(getattr(self, key))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{key} must be finite and non-negative.")

    def _validate_planar_hemodynamics(self):
        if self.planar_terminal_flux_allocation_mode != "equal_terminal":
            raise ValueError(
                "planar_terminal_flux_allocation_mode must be 'equal_terminal'."
            )
        if self.planar_width_scale_mode != "fixed_root_half_width":
            raise ValueError(
                "planar_width_scale_mode must be 'fixed_root_half_width'."
            )
        self._validate_dimensional_hemodynamics(
            inlet_value=self.planar_inlet_total_flux_um2_s,
            inlet_key="planar_inlet_total_flux_um2_s",
            root_size=self.planar_root_half_width_um,
            root_key="planar_root_half_width_um",
            minimum_size=self.planar_min_half_width_um,
            minimum_key="planar_min_half_width_um",
            gamma=self.planar_murray_gamma,
            gamma_key="planar_murray_gamma",
            velocity_min=self.planar_velocity_min_mm_s,
            velocity_max=self.planar_velocity_max_mm_s,
            velocity_prefix="planar",
            inlet_velocity_um_s=lambda flux, half_width: flux / (2.0 * half_width),
        )

    def _validate_volumetric_hemodynamics(self):
        if self.volumetric_terminal_flow_allocation_mode != "equal_terminal":
            raise ValueError(
                "volumetric_terminal_flow_allocation_mode must be 'equal_terminal'."
            )
        if self.volumetric_radius_scale_mode != "fixed_root_radius":
            raise ValueError(
                "volumetric_radius_scale_mode must be 'fixed_root_radius'."
            )
        self._validate_dimensional_hemodynamics(
            inlet_value=self.volumetric_inlet_total_flow_um3_s,
            inlet_key="volumetric_inlet_total_flow_um3_s",
            root_size=self.volumetric_root_radius_um,
            root_key="volumetric_root_radius_um",
            minimum_size=self.volumetric_min_radius_um,
            minimum_key="volumetric_min_radius_um",
            gamma=self.volumetric_murray_gamma,
            gamma_key="volumetric_murray_gamma",
            velocity_min=self.volumetric_velocity_min_mm_s,
            velocity_max=self.volumetric_velocity_max_mm_s,
            velocity_prefix="volumetric",
            inlet_velocity_um_s=lambda flow, radius: (
                flow / (math.pi * radius * radius)
            ),
        )

    def _validate_dimensional_hemodynamics(
        self,
        *,
        inlet_value,
        inlet_key,
        root_size,
        root_key,
        minimum_size,
        minimum_key,
        gamma,
        gamma_key,
        velocity_min,
        velocity_max,
        velocity_prefix,
        inlet_velocity_um_s,
    ):
        inlet_value = float(inlet_value)
        root_size = float(root_size)
        minimum_size = float(minimum_size)
        gamma = float(gamma)
        velocity_min = float(velocity_min)
        velocity_max = float(velocity_max)
        if not math.isfinite(inlet_value) or inlet_value <= 0.0:
            raise ValueError(f"{inlet_key} must be finite and positive.")
        if not math.isfinite(root_size) or root_size <= 0.0:
            raise ValueError(f"{root_key} must be finite and positive.")
        if not math.isfinite(minimum_size) or minimum_size <= 0.0:
            raise ValueError(f"{minimum_key} must be finite and positive.")
        if root_size < minimum_size:
            raise ValueError(f"{root_key} must be at least {minimum_key}.")
        if not math.isfinite(gamma) or gamma <= 0.0:
            raise ValueError(f"{gamma_key} must be finite and positive.")
        if velocity_min < 0.0 or velocity_max < 0.0:
            raise ValueError(
                f"{velocity_prefix} velocity validation bounds cannot be negative."
            )
        if velocity_max > 0.0 and velocity_min > velocity_max:
            raise ValueError(
                f"{velocity_prefix}_velocity_min_mm_s cannot exceed "
                f"{velocity_prefix}_velocity_max_mm_s."
            )
        predicted_terminal_size = root_size * (
            1.0 / float(self.n_terminals)
        ) ** (1.0 / gamma)
        if predicted_terminal_size + 1.0e-12 < minimum_size:
            raise ValueError(
                f"{root_key}, n_terminals, and {gamma_key} predict a terminal "
                f"size below {minimum_key}: "
                f"predicted={predicted_terminal_size:.9g} um, "
                f"minimum={minimum_size:.9g} um."
            )
        inlet_velocity_mm_s = (
            inlet_velocity_um_s(inlet_value, root_size) / 1000.0
        )
        if velocity_max > 0.0 and inlet_velocity_mm_s > velocity_max:
            raise ValueError(
                f"{inlet_key} and {root_key} imply an inlet mean velocity of "
                f"{inlet_velocity_mm_s:.9g} mm/s, above "
                f"{velocity_prefix}_velocity_max_mm_s={velocity_max:.9g}."
            )


def load_generator_config(path):
    """
    read the full generator configuration from a YAML file.

    return `(output, base_cfg, stages)`：
    - `output`: final output path for the SWC file;
    - `base_cfg`: basic configuration for the generator, used as the base for any stage overrides;
    - `stages`: optional staged growth configuration. If empty, the generator will run in normal single-stage mode.
    """
    # ============================
    # Read YAML file and check top-level keys
    # ============================
    path = Path(path)
    data = yaml.safe_load(path.read_text(encoding="utf-8-sig"))
    if data is None:
        data = {}
    if not isinstance(data, dict):
        raise ValueError(f"The YAML root in {path} must be a mapping.")

    # TOP_LEVEL_KEYS = {"output", "config", "stages"}
    unknown = sorted(set(data) - TOP_LEVEL_KEYS)
    if unknown:
        raise ValueError(f"Unknown top-level YAML key(s): {', '.join(unknown)}")

    # ============================
    # Read base configuration under 'config:' key
    # ============================
    raw_config = data.get("config")
    if raw_config is None:
        base_cfg = DCCOConfig()
    elif isinstance(raw_config, dict):
        base_cfg = DCCOConfig(**raw_config)
    else:
        raise ValueError("config must be a YAML mapping.")

    # ===========================
    # Check whether there are any staged growth specifications under 'stages:' key.
    # ===========================
    raw_stages = data.get("stages")
    if raw_stages is None:
        raw_stages = []
    if not isinstance(raw_stages, list):
        raise ValueError("stages must be a YAML list of stage mappings.")

    stages = []
    for index, raw_stage in enumerate(raw_stages, start=1):
        if not isinstance(raw_stage, dict):
            raise ValueError(f"stage {index} must be a YAML mapping.")

        stage = dict(raw_stage)
        name  = str(stage.pop("name", f"stage_{index}"))

        # specify the initial tree for this stage, using the shorthand parser to convert YAML into a TInitSpec.
        tinit = parse_tinit_spec(stage.pop("tinit", None), index)

        # Start with the base configuration, then replace only the values written in this stage.
        stage_values    = base_cfg.as_dict()
        stage_overrides = _normalise_config_mapping(stage)
        stage_values.update(stage_overrides)
        stage_cfg       = DCCOConfig._from_normalised(stage_values)
        if stage_cfg.geometry_mode != base_cfg.geometry_mode:
            raise ValueError(
                f"stage {index} changes geometry_mode from "
                f"{base_cfg.geometry_mode!r} to {stage_cfg.geometry_mode!r}. "
                "A staged run cannot mix planar_2d and volumetric_3d "
                "calculations; use separate generator runs."
            )
        stages.append(StageConfig(name=name, cfg=stage_cfg, tinit=tinit))

    # ===========================
    # Determine the output path for the SWC file.
    # ===========================
    configured_output = data.get("output")
    if configured_output in (None, ""):
        if base_cfg.geometry_mode == "planar_2d":
            output_name = (
                f"ulm_xz_planar_dcco_tree_seed_{int(base_cfg.seed)}.swc"
            )
        else:
            output_name = (
                "ulm_xyz_volumetric_3d_dcco_tree_"
                f"seed_{int(base_cfg.seed)}.swc"
            )
    else:
        output_name = Path(str(configured_output)).name
    output      = SWC_OUTPUT_DIR / output_name

    return output, base_cfg, stages
