# Stage O final-output patch

The isolated `external/sv13o/svMultiPhysics` source inherits the frozen
Stage N lifecycle adapter and base commit. `source_patch.json` in the Stage O
reports records before/after hashes for all 2,430 tracked files. The only
additional changed file is `Code/Source/solver/main.cpp`.

The normal VTU cadence stays controlled by the existing XML setting. The
write condition additionally accepts the existing native stop condition, so
the completed final step is saved even between cadence points. Native restart
already saves on that stop condition. No solver, assembly, integration,
restart format, or PETSc lifecycle operation changes.

The short official case requested STOP_SIM at step 11 with cadence 10 and
wrote exactly VTUs 10 and 11 plus the final native checkpoint. See
`reports/sv1_3o/OFFICIAL_OUTPUT_STOP_acceptance.json` and
`reports/sv1_3o/output_consistency.json` for reload and timestamp evidence.

Build recovery: the existing native Stage N source package contained 1,564
of the 2,430 tracked WSL source files. The initial setup stopped on the first
missing file before patching/building. All existing files matched their WSL
hashes; the 866 missing files were copied into the isolated Stage O source.
`native_source_completion.json` lists those paths; `build_remote.py` verified
the complete source manifest before applying this patch and performing the
clean build. The original native Stage N tree was unchanged.
