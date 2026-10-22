"""Reuse the production solver runner and gates, adapting only the port integrals."""
from pathlib import Path
import importlib.util
import json
import sys

CASE_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CASE_ROOT / 'vendor'))
from pipe_measurements import PipeMeasurements

spec = importlib.util.spec_from_file_location('unchanged_production_runner', CASE_ROOT / 'scripts/solve_remote.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)  # Existing --case and --reference-root CLI parser.
runner.SolutionMeasurements = PipeMeasurements
runner.main()
result = json.loads((runner.args.case / 'reports/execution.json').read_text())
if result['status'] != 'PASS':
    raise SystemExit('Native CFD did not pass its production convergence gates')
