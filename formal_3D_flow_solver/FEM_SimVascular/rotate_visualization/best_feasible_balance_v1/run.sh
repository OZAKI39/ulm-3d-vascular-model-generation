#!/usr/bin/env bash
set -euo pipefail
rotate_package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
rotate_python="${FEM_ROTATE_PYTHON:-}"
if [[ -z "$rotate_python" ]]; then
    for rotate_candidate in \
        "$rotate_package_dir/.venv/bin/python" \
        /home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python \
        "$rotate_package_dir/../.venv/bin/python" \
        python3; do
        if "$rotate_candidate" -c 'import numpy,scipy,pyvista,matplotlib,PIL,imageio_ffmpeg' >/dev/null 2>&1; then
            rotate_python="$rotate_candidate"
            break
        fi
    done
fi
if [[ -z "$rotate_python" ]]; then
    echo 'No Python environment contains the required packages. Install requirements.txt or set FEM_ROTATE_PYTHON.' >&2
    exit 1
fi
rotate_script=render_visualization.py
if [[ "${1:-}" == --validate ]]; then
    rotate_script=validate_results.py
    shift
elif [[ "${1:-}" == --surface-fields ]]; then
    rotate_script=render_surface_fields.py
    shift
elif [[ "${1:-}" == --validate-surface-fields ]]; then
    rotate_script=validate_surface_fields.py
    shift
fi
exec "$rotate_python" "$rotate_package_dir/$rotate_script" "$@"
