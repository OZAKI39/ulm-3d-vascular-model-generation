"""One declared hemisphere, one material and one common loading budget."""
from pathlib import Path
import json, math

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / 'verification/lammps_fragmentation_paradigm_screening_v1'
RESULTS = ROOT / 'results/lammps_fragmentation_paradigm_screening_v1'
LMP = ROOT / 'tools/lammps-install/bin/lmp'
PATTERNS = ('LOCAL_GAUSSIAN_SHEAR', 'OPPOSING_SHEAR', 'VORTEX_TANGENTIAL')
BASELINE = dict(
    model_level='LEVEL_4_QUALITATIVE_EXPLICIT_FRAGMENTATION',
    method_id='LAMMPS_BPM_PERI_LPS', scope='EXPLORATORY_ONLY',
    geometry_type='HEMISPHERE', geometry_id='HEMISPHERE_R562P5UM_CUBIC50UM',
    radius_m=0.0005625, dx_m=0.00005, delta_m=3.01*0.00005,
    K_Pa=100000., G_Pa=1040., rho_kg_m3=1000., s00=.20, alpha=0.,
    load_provider='ANALYTIC_STREAMING_PROVIDER_V1',
    load_budget_N=400.*2*math.pi*0.0005625**2,
    nominal_average_traction_Pa=400.,
    ramp_time_s=.004, hold_time_s=.001, dt_s=5e-7, dump_every=40,
    numerical_drag_rate_s_inv=500., mpi_ranks=2,
    normalization='SUM_SURFACE_NODE_FORCE_MAGNITUDES',
    resolved_fragment_min_nodes=3, wall_interface_model='FIXED_BASE_DEMO',
    load_identity='ENGINEERING_FRAGMENTATION_DEMO_LOAD',
    source='User-authorized engineering demonstration; no actual bubble/CFD magnitude',
    units=dict(length='m', force='N', traction='Pa', mass='kg', time='s'))

def dump(path, value):
    def convert(x):
        if hasattr(x, 'tolist'): return x.tolist()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False, default=convert))
