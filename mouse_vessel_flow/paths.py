"""Project-relative paths: resolution never depends on the shell working directory."""
from pathlib import Path
import json

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CASE_NAMES = ('roi_only', 'best_feasible', 'network_h0', 'mean_2mmps', 'stage_q')

def project_path(relative):
    path = (PROJECT_ROOT / relative).resolve()
    if not path.is_relative_to(PROJECT_ROOT):
        raise ValueError('Path escapes this project: ' + str(relative))
    return path

def case_config(name='roi_only'):
    if name not in CASE_NAMES:
        raise ValueError('Unknown case: ' + str(name))
    return json.loads(project_path(f'cases/{name}/case.json').read_text())

def new_output(path):
    """New computations/postprocessing require an explicit, previously absent directory."""
    path = Path(path).expanduser().resolve()
    frozen = [PROJECT_ROOT / n for n in ('inputs','cases','configs','archive','validation/recorded_cases')]
    frozen += [PROJECT_ROOT/'results'/n for n in CASE_NAMES]
    if any(path == p or path.is_relative_to(p) for p in frozen):
        raise ValueError('Output would overlap retained scientific inputs/results')
    path.mkdir(parents=True, exist_ok=False)
    return path
