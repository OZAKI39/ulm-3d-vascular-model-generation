#!/usr/bin/env python3
"""s1-10 compact C system with 2 mm mold and balanced vertical clearance."""
import argparse
from pathlib import Path
import yaml
from vascular_processing.compact_mold_revision import run as run_dimensions
from vascular_processing.compact_c_pipeline import run as run_caps

def run(config=None, out=None, do_slice=True, do_step=True, do_render=True):
    # Explicitly selecting the old cap configuration still reproduces V2.
    legacy = config and yaml.safe_load(Path(config).read_text()).get('method') == 'COMPACT_C_RETAINER_V2'
    return (run_caps if legacy else run_dimensions)(config, out, do_slice, do_step, do_render)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config')
    p.add_argument('--output')
    p.add_argument('--geometry-only',action='store_true',help='Skip real slices; remains NEEDS_ADJUSTMENT')
    p.add_argument('--no-step',action='store_true')
    p.add_argument('--no-render',action='store_true')
    p.add_argument('--print-plate',action='store_true',help='Pack the latest accepted 2 mm mold, core and four C retainers on one plate and slice a native Bambu 3MF')
    a=p.parse_args()
    if a.print_plate:
        if a.config or a.geometry_only or a.no_step or a.no_render:
            p.error('--print-plate uses the existing accepted models and a complete native slice; only --output may be combined with it')
        from vascular_processing.all_abs_print_plate import run as run_plate
        s=run_plate(a.output)
        raise SystemExit(0 if s['status']=='READY_FOR_BAMBU_IMPORT_PROCESS_UNCALIBRATED' else 2)
    s=run(a.config,a.output,not a.geometry_only,not a.no_step,not a.no_render)
    raise SystemExit(0 if not s['failures'] else 2)
