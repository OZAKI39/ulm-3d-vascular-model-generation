"""Render the configured default Y-clot point cloud from accepted data only."""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

def main():
    from pd_clot.lammps_fragmentation.point_cloud_delivery import render
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stills',action='store_true',help='Render representative exact-state PNG/PDFs only')
    args=p.parse_args()
    config=json.loads((ROOT/'references/y_clot_default_render.json').read_text())
    assert config['default_mode']=='POINT_CLOUD' and not config['smooth_shell_default']
    if 'accepted_validation' in config:
        status=json.loads((ROOT/config['accepted_validation']).read_text())
        assert status['DETACHMENT_TRIGGER_CORRECT']=='PASS' and status['HANDOFF_VELOCITY_CONTINUITY']=='PASS'
    render(ROOT/config['data_source'],config['output_basename'],args.stills)

if __name__=='__main__':main()
