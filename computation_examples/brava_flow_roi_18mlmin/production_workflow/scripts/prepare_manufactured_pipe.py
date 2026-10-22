"""Prepare one conforming high-fidelity surface; TetGen remains the only volume mesher."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pyvista as pv
from pipe_surface_qa import surface_qa

ROOT = Path(__file__).resolve().parents[1]


def annulus(inner, outer):
    """Join two concentric rings in angular order, sharing their existing IDs."""
    faces, i, j = [], 0, 0
    ni, no = len(inner), len(outer)
    while i < ni or j < no:
        # Compare exact integer angular fractions, avoiding a geometric tolerance.
        left, right = (i + 1) * no, (j + 1) * ni
        a, b = inner[i % ni], outer[j % no]
        if left == right:
            faces.extend([[a, b, outer[(j + 1) % no]], [a, outer[(j + 1) % no], inner[(i + 1) % ni]]])
            i += 1; j += 1
        elif left < right:
            faces.append([a, b, inner[(i + 1) % ni]]); i += 1
        else:
            faces.append([a, b, outer[(j + 1) % no]]); j += 1
    return faces


def main(radius_name, strategy):
    manufactured = ROOT/'manufactured_pipe'
    if radius_name == 'R1' and strategy in ('P1','P2'):
        previous = 'P0' if strategy == 'P1' else 'P1'
        gate = json.loads((manufactured/('R1_'+previous)/'cfd/analysis/CALIBRATION.json').read_text())
        assert gate['status']=='FAIL_ENGINEERING_CALIBRATION', 'Refine only after the previous strategy failed accuracy'
    if radius_name == 'R2':
        gate = json.loads((manufactured/('R1_'+strategy)/'cfd/analysis/CALIBRATION.json').read_text())
        assert gate['status']=='PASS', 'R2 requires a qualified R1 winner'
        assert not list((ROOT/'straight_pipe/manufactured').glob('R2_*')), 'Only one unchanged R2 policy is allowed'
    frozen = json.loads((ROOT / 'straight_pipe/BENCHMARK_PROTOCOL.json').read_text())
    radius = frozen['radii'][next(k for k in frozen['radii'] if k.startswith(radius_name))]
    policy = dict(frozen['policies'][strategy])
    # Fixed across strategies and radii. R2 uses the identical normalized rule.
    nt, nz, nr = 512, 30, 16
    policy.update(angular_vertices=nt, axial_spacing_over_R=12 / nz,
                  cap_radial_rings=nr, surface_remesh=strategy == 'P2')
    if strategy == 'P2':
        capability = json.loads((manufactured/'NATIVE_LAYER_CAPABILITY.json').read_text())
        policy['use_mmg'] = capability['conditional_P2_configuration']['use_mmg']
        policy['surface_configuration_evidence'] = 'manufactured_pipe/NATIVE_LAYER_CAPABILITY.json'
    name = radius_name + '_' + strategy
    directory = ROOT / 'straight_pipe/manufactured' / name
    directory.mkdir(parents=True, exist_ok=False)
    length = 12 * radius
    angles = np.arange(nt) * 2 * np.pi / nt
    points = [[radius * np.cos(t), radius * np.sin(t), z]
              for z in np.linspace(0, length, nz + 1) for t in angles]
    faces, tags = [], []
    for k in range(nz):
        for j in range(nt):
            a, b = k * nt + j, k * nt + (j + 1) % nt
            c, d = b + nt, a + nt
            faces.extend([[a, b, c], [a, c, d]]); tags.extend([1, 1])
    for tag, z, outer_offset in [(2, 0., 0), (3, length, nz * nt)]:
        center = len(points); points.append([0., 0., z])
        rings = []
        for k in range(1, nr):
            count = nt * k // nr
            ring = list(range(len(points), len(points) + count)); rings.append(ring)
            points.extend([[radius * k / nr * np.cos(t), radius * k / nr * np.sin(t), z]
                           for t in np.arange(count) * 2 * np.pi / count])
        rings.append(list(range(outer_offset, outer_offset + nt)))
        cap = [[center, rings[0][j], rings[0][(j + 1) % len(rings[0])]] for j in range(len(rings[0]))]
        for a, b in zip(rings[:-1], rings[1:]):
            cap.extend(annulus(a, b))
        if tag == 2:
            cap = [triangle[::-1] for triangle in cap]
        faces.extend(cap); tags.extend([tag] * len(cap))
    grid = pv.PolyData(np.asarray(points), np.c_[np.full(len(faces), 3), faces].ravel())
    grid.cell_data['ModelFaceID'] = np.asarray(tags, np.int32)
    grid.save(directory / 'model.vtp')
    previous = json.loads((ROOT / ('straight_pipe/' + name + '/INPUT.json')).read_text())
    config = dict(previous, name='manufactured/' + name, policy=policy,
                  h_bulk_m=radius * policy['h_bulk_over_R'],
                  expected_max_PLC_radial_chord_error_m=radius * (1 - np.cos(np.pi / nt)),
                  model_sha256=hashlib.sha256((directory / 'model.vtp').read_bytes()).hexdigest(),
                  surface_design='512 shared angular vertices; 30 axial intervals; caps have 16 radial rings with 32*k vertices. Same wall/cap geometry across P0/P1/P2.',
                  volume_mesher='Installed native SimVascular TetGen only',
                  design_reason='Analytic facet sagitta <=0.0118 um for both fixed radii; cap radial resolution controls quadratic nodal trace error without flux scaling.')
    (directory / 'INPUT.json').write_text(json.dumps(config, indent=2) + '\n')
    qa = surface_qa(directory / 'model.vtp', radius, length)
    (directory / 'SEAM_QA.json').write_text(json.dumps(qa, indent=2) + '\n')
    assert qa['status'] == 'PASS', qa
    print(json.dumps(dict(case=name, points=grid.n_points, triangles=grid.n_cells,
                          expected_chord_um=config['expected_max_PLC_radial_chord_error_m'] * 1e6,
                          seam_qa=qa['status']), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('radius', choices=['R1', 'R2']); parser.add_argument('strategy', choices=['P0', 'P1', 'P2'])
    args = parser.parse_args(); main(args.radius, args.strategy)
