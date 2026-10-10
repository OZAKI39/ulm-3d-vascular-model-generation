"""Exploratory LPS fragmentation through official LAMMPS BPM topology."""
from .geometry_adapter import prepare_geometry
from .material import prepare_material
from .streaming_patterns import build_streaming_load
from .lammps_writer import write_lammps_case
from .runner import run_lammps
from .fragment_analysis import analyze_fragments, export_result

def render(case):
    from .visualization import render_case
    return render_case(case)
