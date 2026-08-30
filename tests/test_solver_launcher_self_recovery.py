from __future__ import annotations

from utils.cfd_flow.repaired_topology_forensics import (
    scaled_seeder_preflight_script,
    wsl_script_file_command,
)


def test_solver_launcher_is_self_locating_script_file_transport() -> None:
    script = scaled_seeder_preflight_script()
    command = wsl_script_file_command("/mnt/e/run/launcher_preflight.sh")

    assert "${BASH_SOURCE[0]}" in script
    assert 'cd "$SCRIPT_DIR"' in script
    assert command[-2:] == ["/bin/bash", "/mnt/e/run/launcher_preflight.sh"]
    assert "-lc" not in command
