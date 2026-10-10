"""Run the unchanged F5 renderer in an isolated output directory with the archived Arial font."""
from pathlib import Path
import argparse
import os

ROOT = Path('/workspace/mevo_acoustic_microbubbles')
PROJECT = ROOT/'repo/computation_examples/mevo_acoustic_microbubbles'
F5 = PROJECT/'results/fast_pzt/f5'

ap=argparse.ArgumentParser()
ap.add_argument('--mode', choices=['figures','preview','video'], required=True)
ap.add_argument('--output-label', required=True)
args=ap.parse_args()
assert args.output_label and '/' not in args.output_label and args.output_label not in ['.','..']
out=ROOT/'renders'/args.output_label
out.mkdir(parents=True, exist_ok=False)
for name in ['science_v1','ON_OFF_FLOW_COMPARISON.json']:
    (out/name).symlink_to(os.path.relpath(F5/name,out))
from scripts.fast_pzt_f5 import style
style.FONT=ROOT/'references/fonts/arial.ttf'
assert style.FONT.is_file()
from scripts.fast_pzt_f5 import render
render.BASE=out
render.main('science_v1',args.mode)
print(out/'visualization')
