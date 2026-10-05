# Stage 3 quickstart

Run from `/workspace/mevo_acoustic_microbubbles` on the current server (`ssh -p 4159 root@50.115.148.16`), or the local mirror `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles` with the recorded dependencies/data installed.

The current instance requires long jobs under supervisor. The finite helper uses autostart=false and autorestart=false; it does not resume earlier-stage jobs.

```bash
.venv/bin/python scripts/stage3/submit.py secondary scripts/stage3/acceptance.py secondary
.venv/bin/python scripts/stage3/submit.py clustering scripts/stage3/acceptance.py clustering
.venv/bin/python scripts/stage3/submit.py flow scripts/stage3/acceptance.py flow
.venv/bin/python scripts/stage3/submit.py history scripts/stage3/history.py
```

Wait for successful STATE.json files. Then:

```bash
.venv/bin/python -m scripts.stage3.history_report
.venv/bin/python -m scripts.stage3.gates
.venv/bin/python scripts/stage3/submit.py cylinder scripts/stage3/cylinder.py
```

The delivered acoustic STATE.json already points to the matching cached basis. Do not rerun `acoustic.py` for changed phase/amplitude. `configs/stage3/cylinder.yaml` reuses the validated Stage 2.5 table, with no new production table construction.

The same tracker can be run directly with `mevo-mb bubbles run --config configs/stage3/cylinder.yaml --basis <basis-path> --pattern shifted`. The six gates are required; a demonstration override is not formal validation. After changing evidence, regenerate the hashed registry rather than editing a PASS label manually.

`flow.backend` selects zero, poiseuille or vtk. VTK requires files and optionally increasing times_s, with explicit coordinate/velocity scales if not SI. Out-of-domain probes reject. `particles.contact_search_backend` selects ckdtree or brute. `particles.near_field_policy: simulator` enables warning-only model regimes, while table bounds and numerical failures remain fatal.

```bash
.venv/bin/pytest -q
.venv/bin/pytest -q -m deprecated_reference
```

The first command excludes archived Zhang literature assertions from CI; the second audits them independently. Result identities, scope and reproduction commands are in the final report. Production has no history-force option.
