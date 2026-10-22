# BraVa development instructions

- Read [ACTIVE_DEVELOPMENT_SCOPE.md](ACTIVE_DEVELOPMENT_SCOPE.md) first. Follow the current user request; the second cleanup stops after tests, documentation and its independent commit.
- Default reading is limited to production particle source, manufactured Poiseuille reference, the future Gmsh validation location, current authoritative reports and tests. Do not recursively read retired implementations/reports/logs unless the user explicitly requests a historical audit.
- Historical experimental implementations are recovered from Git history, not retained in the active tree. See [RETIRED_EXPERIMENTS_INDEX.md](RETIRED_EXPERIMENTS_INDEX.md). Do not create code/log archives or rewrite Git history.
- Preserve production physics, all current shared tests, formal geometry/flow, original population seeds and accepted states, source/fixed-window contracts, necessary SHA identities and compatibility symlinks. Unknown cleanup candidates are KEEP. Shared repositories must not be staged wholesale.
- Scientific display filtering never changes accounting. Old 1565 results include 141 censored and 2 numerical failures; they are not a new fixed-window run. Conservation, CPU/GPU agreement and smooth animation do not prove local accuracy.
- Retain native SimVascular/TetGen and svMultiPhysics adapters. Future validation candidates stay separate until gates pass. Stop at failed gates; no automatic fallback/old queue/full-domain refinement.
- Read `/etc/vast-agents-guide.md` before server actions. Use project supervision for long jobs, check actual CPU/GPU resources, and do not change other projects' services or system settings.
- For a new BraVa boundary design, use centerline-resistance outlet pressures as instructed. Do not alter the historical calibrated field or call it a new non-calibrated result.
