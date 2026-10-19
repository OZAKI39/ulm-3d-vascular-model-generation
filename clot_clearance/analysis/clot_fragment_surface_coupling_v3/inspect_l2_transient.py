"""Read-only diagnosis of the first observed, transient L2-loaded component."""
from common import *
from current_surface import CurrentSurface, CLASS_NAMES, OPEN


def main():
    observation = json.loads((O/'INTERMEDIATE_L2_ACCEPTED_STATE_OBSERVATION.json').read_text())
    cohort = np.array(observation['L2_loaded_components'][0]['original_atom_ids'], dtype=int)-1
    ref = load_npz(V2_O/'provenance/reference.npz')
    surf = CurrentSurface(ref, 'uniform')  # The material is biphasic; the imposed flow is uniform.
    keydir = O/'cases/biphasic/keys'
    records = []
    tracked = None
    for name in ['frame_008554.npz', 'frame_008715.npz', 'frame_008890.npz']:
        path = keydir/name
        digest = sha(path)
        s = load_npz(path)
        assert np.array_equal(s['atom_ids'], np.arange(1, surf.N+1))
        z = surf.evaluate(s['x'], s['velocity'], s['d'], float(s['time']))
        top = z['topology']
        faces = np.flatnonzero(np.isin(surf.owner, cohort))
        active = faces[z['active'][surf.L+faces]]
        if tracked is None:
            tracked = active.copy()
            assert len(tracked) > 0
        components = []
        for label in np.unique(top['labels'][cohort]):
            members = np.flatnonzero(top['labels'] == label)
            components.append(dict(original_atom_ids=(members+1), nodes=len(members),
                                   anchor_reachable=bool(top['attached'][label])))
        details = []
        for f in tracked:
            owner = surf.owner[f]; g = surf.L+f
            gradient = top['gradients'][surf.axis[f]].getrow(owner)
            weights = {int(j):float(v*.5*surf.sign[f]) for j,v in zip(gradient.indices,gradient.data)}
            weights[int(owner)] = weights.get(int(owner),0.)+1.
            details.append(dict(L2_patch_index=int(f), owner_atom_id=int(owner+1),
                reference_axis=int(surf.axis[f]), reference_sign=float(surf.sign[f]),
                component_boundary=bool(top['boundary'][f]),
                reference_fit_valid=bool(top['reference_fit_valid'][owner]),
                current_cell_valid=bool(z['l2_cell_valid'][owner]),
                mapping_weights_admissible=bool(top['weight_ok'][f]),
                proposed_mapping_weights_by_original_atom_id={str(j+1):v for j,v in weights.items()},
                proposed_minimum_weight=min(weights.values()),
                proposed_sum_absolute_weights=sum(abs(v) for v in weights.values()),
                alive_crossing_bonds=int(z['l2_bridge_count'][f]),
                L1_overlap_guard=bool(z['l2_overlap_guard'][f]),
                geometric=bool(z['geometric'][g]), aperture_m=float(z['aperture_m'][g]),
                opening_threshold_m=OPEN, face_class=CLASS_NAMES[int(z['face_class'][g])],
                active=bool(z['active'][g]), effective_area_m2=float(z['area'][g]*z['weight'][g]),
                face_force_N=z['face_force'][g],
                face_power_W=float(z['face_force'][g]@z['face_velocity'][g])))
        records.append(dict(source=str(path.relative_to(R)), source_sha256=digest,
            time_s=float(s['time']), step=int(s['step']), cohort_components=components,
            cohort_L2_active=len(active), boundary_face_classes={name:int(np.sum(
                top['boundary'][faces] & (z['face_class'][surf.L+faces] == i)))
                for i,name in enumerate(CLASS_NAMES)}, tracked_patches=details))
        assert sha(path) == digest
        print(name, float(s['time']), len(active), details, flush=True)
    dump(O/'L2_TRANSIENT_LOAD_DIAGNOSIS.json',dict(status='READ_ONLY_ACCEPTED_STATE_DIAGNOSIS',
        cohort_original_atom_ids=cohort+1, snapshots=records,
        limitation='Discrete eligibility and visibility reasons only; absence of a loaded patch is not evidence of zero physical fluid resistance. No dynamics rerun or surface rule changed.'))


if __name__ == '__main__':
    main()
