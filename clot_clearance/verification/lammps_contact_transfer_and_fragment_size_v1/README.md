# LAMMPS contact transfer and fragment size

Start with REVIEW_LAMMPS_CONTACT_TRANSFER_AND_FRAGMENT_SIZE_V1.md and
FINAL_STATUS.json. This round is complete, EXPLORATORY_ONLY, and stopped for
human review. Four new hemisphere solves were executed: two mandatory and two
conditionally authorized 31.6228 Pa points. No historical solve was repeated.

The old 37-case map is the NATIVE_1P35DX model. New points and the reused
Gaussian100 contact candidate belong to MODIFIED_1P0DX. Unknown map regions
remain unknown. The official source/install/default is unchanged; 1.0dx has
conditional exploratory candidate status only.

## Statistics

Every analyzed frame was reread from the native dump. Broken-bond events were
deduplicated and mapped to the initial graph, clot-only components rebuilt, and
surviving interface bonds used for attachment. Independent results agree with
the previously accepted components and the native final surviving topology.

Use SIZE_THRESHOLD_POLICY.md for exact 34/114-node cutoffs, reference-volume
equivalent diameter, PCA thickness flags and persistence semantics.
PER_FRAME_DETACHED_COMPONENTS.csv includes only detached components; the original
per-case analysis/components.csv also includes attached components and is retained
locally. All one/two-node debris remains in the data. FRAGMENT_LIFETIMES.csv
contains exact-membership episodes, including attached episodes; consult
detached_frame_count. FRAGMENT_SPLIT_RELATIONS.csv links every parent and child.
Membership episode counts must not be interpreted as independent embolus counts.

## Reproduction

The interpreter is
`/home/lzy/projects/computation_examples/_shared/python_environment/bin/python`.
The new four Python modules extend pd_clot/lammps_fragmentation. Existing runner,
provider, graph/classifier and geometry definitions are reused and hash protected.
collect.py, render.py, audit_results.py and write_review.py are postprocessing.
run_stage0.py is native-trajectory reanalysis; its output-directory guards prevent
overwriting the saved audit. The mandatory/conditional run scripts document the
executed procedure, not an authorization to repeat it. Existing case directories
and the four-launch ceiling block accidental reruns. No new package was installed.

Figures and the single synchronized movie are under
results/lammps_contact_transfer_and_fragment_size_v1. The movie uses 61 saved
states including the 0.500 ms 49-node example; there is no interpolation, magnified
displacement or enlarged debris. Both local and full-cloud views use shared fixed
cameras across cases. Playback timing is illustrative, not experimental time.

## Evidence package

HASH_MANIFEST.json records absolute paths, sizes and SHA256, plus inclusion in
evidence.zip. Large native trajectory/surviving-bond/restart files and redundant
all-component CSV files remain at their original locations; their hashes are
retained. The archive includes compact accepted states for new cases and the
reused Gaussian100 movie reference, native break-event evidence, graph/lifetime
tables, inputs/logs, plot data, figures, movie, and source snapshots.

The old executable/build tree is not duplicated. Its exact hash, source SHA,
patch and prior static contact/restart checks are included as provenance.
Complete previous two rounds (1312 protected files, including old source and
executables) were verified unchanged. HISTORY_PROTECTION.json and directory
digests document that check. No commit, push or PR was created.

The archive, its checksum and DELIVERY_CHECK.json are excluded from the manifest
to avoid circular hashes. The manifest excludes itself and live packaging logs.
DELIVERY_CHECK.json is written after ZIP CRC and member/disk hashes pass.
