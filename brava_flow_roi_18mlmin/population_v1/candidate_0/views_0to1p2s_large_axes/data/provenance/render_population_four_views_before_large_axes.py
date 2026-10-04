"""Render four candidate-0 views of the same recorded 0–1.2 s population."""
from pathlib import Path
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'candidate_0/views_0to1p2s'
VIEWS = [(240, 25.), (330, 25.), (60, 25.), (150, 25.)]


def main():
    records = []
    for azimuth, elevation in VIEWS:
        name = f'az{azimuth:03d}_el{elevation:.1f}'.replace('.', 'p')
        target = OUT/name
        cmd = [sys.executable, '-u', '-B', str(ROOT/'scripts/render_view05_flow_axes.py'),
               '--output-dir', str(target), '--azimuth-offset', str(azimuth),
               '--elevation', str(elevation), '--end-time', '1.2', '--angle-label',
               '--movie-name', f'microbubble_{name}_0to1p2s.mp4']
        print('START_VIEW', name, flush=True)
        subprocess.run(cmd, check=True)
        manifest = json.loads((target/'data/render_manifest.json').read_text())
        records.append(dict(name=name, azimuth_offset_deg=azimuth, elevation_deg=elevation,
                            world_azimuth_deg=manifest['layout']['world_azimuth_deg'],
                            video=f"{name}/{manifest['movie_path']}",
                            frames=manifest['frames'], fps=manifest['fps'],
                            video_duration_s=manifest['video_duration_s'],
                            source_manifest=f'{name}/data/render_manifest.json'))
    (OUT/'data').mkdir(exist_ok=True)
    (OUT/'data/views.json').write_text(json.dumps(records, indent=2)+'\n')
    print('ALL_FOUR_VIEWS_RENDERED', json.dumps(records), flush=True)


if __name__ == '__main__':
    main()
