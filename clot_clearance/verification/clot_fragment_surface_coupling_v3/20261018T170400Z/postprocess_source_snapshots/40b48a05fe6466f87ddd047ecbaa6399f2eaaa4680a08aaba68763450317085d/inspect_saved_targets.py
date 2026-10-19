"""Read-only target evidence from an immutable periodic or drive-off state."""
from common import *
from current_surface import CurrentSurface
from analyze import graph, l2_resolution
import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--key', help='Existing frame_NNNNNN.npz or upper_drive_off.npz; defaults to the latest saved periodic frame')
    args = parser.parse_args()
    keys = O / 'cases/biphasic/keys'
    source = keys / args.key if args.key else sorted(keys.glob('frame_*.npz'))[-1]
    assert source.parent == keys and (source.name.startswith('frame_') or source.name == 'upper_drive_off.npz') and source.suffix == '.npz'
    source_hash = sha(source)
    script_hash = sha(Path(__file__))
    script_archive = O / 'postprocess_source_snapshots' / script_hash / Path(__file__).name
    script_archive.parent.mkdir(parents=True, exist_ok=True)
    if script_archive.exists():
        assert sha(script_archive) == script_hash
    else:
        script_archive.write_bytes(Path(__file__).read_bytes())
    s = load_npz(source)
    ref = load_npz(O / 'provenance/reference.npz')
    assert np.array_equal(s['atom_ids'], np.arange(1, len(ref['X']) + 1))
    labels, sizes, attached = graph.topology(ref, s['d'])
    surf = CurrentSurface(ref, 'uniform')
    z = surf.evaluate(s['x'], s['velocity'], s['d'], float(s['time']))
    codes, names, margins, near = l2_resolution(surf, z)
    assert np.max(abs(z['force'] - s['external_force'])) < 1e-15
    b = ref['bonds']
    targets = []
    for label in np.flatnonzero(~attached & (sizes >= 20) & (sizes <= 200)):
        ids = np.flatnonzero(labels == label)
        member = np.zeros(len(ref['X']), bool)
        member[ids] = True
        internal = member[b].all(1) & (s['d'] < 1)
        external = (member[b].sum(1) == 1) & (s['d'] < 1)
        assert not external.any() and internal.sum() >= len(ids) - 1
        mask = member[surf.tri[:, 0]] & z['active']
        dual = member[surf.owner] & (codes >= 0)
        center = s['x'][ids].mean(0)
        force = s['external_force'][ids]
        area = float(np.sum(z['area'][mask] * z['weight'][mask]))
        targets.append(dict(original_atom_ids=ids + 1, nodes=len(ids), internal_alive_bonds=int(internal.sum()),
                            remaining_live_external_connections=0, anchor_reachable=False,
                            L1_active=int(mask[:surf.L].sum()), L2_active=int(mask[surf.L:].sum()),
                            effective_loaded_area_m2=area,
                            loading_status='PARTIAL_APPROXIMATE_SURFACE_LOADING' if area else 'FLUID_LOADING_UNRESOLVED',
                            current_COM_m=center, COM_velocity_m_s=s['velocity'][ids].mean(0),
                            resultant_force_N=force.sum(0),
                            torque_about_COM_Nm=np.cross(s['x'][ids] - center, force).sum(0),
                            nodal_power_W=float(np.sum(s['velocity'][ids] * force)),
                            L2_first_condition_counts={name: int(np.sum(dual & (codes == i))) for i, name in enumerate(names)},
                            rejected_mapping_near_strict_limit=int(np.sum(dual & near))))
    assert sha(source) == source_hash
    destination = O / 'intermediate_target_observations' / (source.stem + '.json')
    dump(destination, dict(status='INTERMEDIATE_ACCEPTED_STATE_OBSERVATION', source=str(source.relative_to(R)),
                           source_sha256=source_hash, script_sha256=script_hash,
                           archived_script=str(script_archive.relative_to(R)),
                           time_s=float(s['time']), accepted_step=int(s['step']), target_count=len(targets), targets=targets,
                           scope='Actual archived coordinates, original surviving bonds and recorded nodal force. Surface replay agrees with the saved force. Not a final count, survival assessment, or additional physical run.'))
    print(json.dumps(dict(path=str(destination), time_ms=float(s['time']) * 1000,
                          targets=[{k: q[k] for k in ['nodes', 'internal_alive_bonds', 'L1_active', 'L2_active',
                                    'effective_loaded_area_m2', 'resultant_force_N', 'nodal_power_W', 'loading_status']} for q in targets]),
                     default=cv), flush=True)


if __name__ == '__main__':
    main()
