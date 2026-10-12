# F10H reproduction

Use the existing isolated project environment. Offline work needs NumPy/SciPy,
h5py and Matplotlib; GPU trajectories additionally need the existing Warp 1.18.
Set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src:.`.

Completed science (do not repeat to regenerate figures):

1. baseline: F10G complete delivery, history and 13 input hashes.
2. checks_v1: explicit LUB2/LUB4 static, original mode, control/Graph/fail-stop gates.
3. short_v1: two real LUB4/10us windows from original initial states.
4. long_v1: only O/1102 attempted 10ms; stopped at first 6.370ms rejection.
5. dt5_v1: the two corresponding LUB4/5us short windows.
6. geometry_v1: independent two/four, conditional saved32 and F10G false reactions.
7. b2_short_v1: full static environment, control/fail-stop checks and O/1102 5ms.
8. b2_dt5_v1: same geometric model/stiffness, O/1102 5ms at 5us.
9. performance_v1: matched 0.2ms LUB2/LUB4 diagnostic cost, three warm repeats.
10. checks_v2: recheck after the explicit mode registry was generalized for B2.

Nine server processes cover the ten numbered stages because baseline is local.
All server tasks use immutable Git archives and unique Supervisor labels.
`remote <module> --label <unique>` refuses dirty executable source, occupied GPU
and concurrent MEVO jobs. It never restarts a failed attempt or stops another job.
`archive` audits every remote source against its exact Git archive, verifies actual
inputs, copies evidence, then verifies result hashes. Wait for it to finish.

Offline, after successful complete archive:

```
python -m scripts.fast_pzt_f10h.analyze
python -m scripts.fast_pzt_f10h.events
python -m scripts.fast_pzt_f10h.figures
python -m scripts.fast_pzt_f10h.delivery
python -m scripts.fast_pzt_f10h.verify
```

`delivery` updates the external Chinese progress document only by an audited
minimal insertion; it refuses to overwrite an independently changed document.
Actual PNG and PDF review is recorded separately in VISUAL_REVIEW.json.
`verify` stores derived per-run admissions without changing the run-time raw HDF5
or its PENDING analysis fields. Original identities, stages, inputs, errors and
the final accepted state before rejection are preserved.

No production GPU batch, B2 10ms, new independent seed, full OFF/S/O trio or AI
training was run. Original M0 remains the default. Same-model shape and event
gates are distinct from contact legality and from physical calibration.
