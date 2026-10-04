#!/usr/bin/env bash
# Run from the existing clot_clearance project root. Only new candidate paths are written.
set -euo pipefail
PY=/home/lzy/projects/computation_examples/_shared/python_environment/bin/python
O="$PWD/verification/pfpd_candidate_v1"
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMBA_CACHE_DIR="$O/provenance/runtime_cache/reproduce" MPLCONFIGDIR="$O/provenance/runtime_cache/matplotlib"
# Resource acquisition (already performed and pinned):
# git clone --filter=blob:none https://github.com/kaipartmann/Peridynamics.jl.git "$O/external/Peridynamics.jl"
# git -C "$O/external/Peridynamics.jl" checkout --detach b6e86563845a83b5cb06c805af19a817bde6b106
# (cd "$O/external/poppler" && apt-get download poppler-utils && dpkg-deb -x poppler-utils*.deb local)
# "$O/external/poppler/local/usr/bin/pdftotext" -layout /home/lzy/projects/2604.09215v1.pdf "$O/resources/paper_pdftotext.txt"
# "$PY" -B "$O/resources/fetch_resources.py"   # Official Julia 1.10.12 tarball checked against official SHA256.
# "$PY" -B "$O/resources/audit_author_resources.py"  # Author-source snapshot; network state can change.
export JULIA_DEPOT_PATH="$O/external/julia_depot" JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1
JULIA="$O/external/julia-1.10.12/bin/julia"
# Initial independent environment setup command:
# "$JULIA" --startup-file=no --project="$O/external/julia_env" -e 'using Pkg; Pkg.develop(path="verification/pfpd_candidate_v1/external/Peridynamics.jl"); Pkg.instantiate(); Pkg.precompile()'
# Use the saved comparisons/Manifest.toml and Project.toml in external/julia_env for the exact resolved package set.
"$PY" -B scripts/verify_pfpd_baqp.py
"$JULIA" --startup-file=no --project="$O/external/julia_env" "$O/comparisons/official_oracle.jl"
"$PY" -B scripts/compare_pfpd_julia.py
# Julia has exited before Python process pools start.
"$PY" -B scripts/run_pfpd_compatibility.py
"$PY" -B -m pytest tests -q -p no:cacheprovider --basetemp="$O/tests/reproduce_tmp"
"$PY" -B scripts/visualize_pfpd_candidate.py
