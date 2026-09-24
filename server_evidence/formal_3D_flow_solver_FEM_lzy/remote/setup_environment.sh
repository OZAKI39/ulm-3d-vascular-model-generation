#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
env_prefix="$project_root/remote/.env"
conda_executable="${FEM3D_CONDA:-/opt/miniforge3/bin/conda}"
export CONDA_PKGS_DIRS="$project_root/remote/.cache/conda-pkgs"
export XDG_CACHE_HOME="$project_root/remote/.cache"
export CONDA_CHANNEL_PRIORITY=strict
export CONDA_NO_PLUGINS=true
mkdir -p "$project_root/outputs/stage00/environment"
printf 'Remote project directory: %s\nEnvironment prefix: %s\n' "$project_root" "$env_prefix"
if [[ "${1:-}" == "--dry-run" ]]; then
    exec "$conda_executable" env create --solver=classic --prefix "$env_prefix" --file "$project_root/remote/environment.yml" --dry-run --json
fi
if [[ $# -ne 0 ]]; then printf 'Only --dry-run is supported\n' >&2; exit 2; fi
if [[ -e "$env_prefix" ]]; then
    test -f "$env_prefix/.fem3d-environment-spec-sha256"
    sha256sum --check "$env_prefix/.fem3d-environment-spec-sha256"
else
    "$conda_executable" env create --solver=classic --prefix "$env_prefix" --file "$project_root/remote/environment.yml" --yes
    sha256sum "$project_root/remote/environment.yml" > "$env_prefix/.fem3d-environment-spec-sha256"
fi
"$conda_executable" list --prefix "$env_prefix" --explicit > "$project_root/outputs/stage00/environment/conda-linux-64.explicit.txt"
"$conda_executable" env export --prefix "$env_prefix" > "$project_root/outputs/stage00/environment/resolved_environment.yml"
"$env_prefix/bin/python" -B -c 'import dolfinx; assert dolfinx.__version__.startswith("0.11."); print(dolfinx.__version__)'
