"""Verify final remote inputs and imported execution paths without rerunning physics."""
import sys,importlib,importlib.metadata as md
from .common import *

def main():
    verify()
    source={str(p.relative_to(ROOT)):sha(p) for folder in ['scripts/straight_tube_exp1/f10af','configs/straight_tube_exp1/f10af','tests/f10af'] for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in str(p)}
    imports={name:str(Path(importlib.import_module(name).__file__).resolve()) for name in ['mevo_mb.fast_sim.accelerated_water','mevo_mb.fast_sim.accelerated_sequential','scripts.straight_tube_exp1.f10af.slow']}
    assert all(str(ROOT) in path for path in imports.values())
    assert 'f10af_observable_20261019' in str(Path(sys.prefix).resolve())
    write(BASE/'FINAL_EXECUTION_SEAL.json',dict(identity=identity(),execution_root=str(ROOT),python=sys.executable,environment=str(Path(sys.prefix).resolve()),source_sha256=source,imports=imports,versions={name:md.version(name) for name in ['numpy','scipy','warp-lang','smac','botorch']},protected_inputs_pass=True,scope='Final reproducibility source snapshot; does not replace individual scientific execution commits. No new trajectory.'))
if __name__=='__main__':main()
