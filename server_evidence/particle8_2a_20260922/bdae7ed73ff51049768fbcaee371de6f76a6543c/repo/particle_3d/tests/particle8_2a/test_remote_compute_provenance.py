import socket
import pytest
from particle_3d.particle82_provenance import require_remote


def test_unattested_host_cannot_claim_remote_compute():
    with pytest.raises(ValueError):require_remote({'role':'DEVELOPMENT','hostname':socket.gethostname()},socket.gethostname())
    with pytest.raises(ValueError):require_remote({'role':'HEAVY_COMPUTE_PRIMARY','hostname':'another-host'},'another-host')
