import numpy as np
from particle_3d.particle82a_pipeline import atomic_npz


def test_saved_path_array_does_not_collide_with_archive_filename(tmp_path):
    path=np.arange(12,dtype=float).reshape(3,4);target=tmp_path/'scene.npz'
    atomic_npz(target,path=path,attempts=np.zeros((0,11)))
    assert np.array_equal(np.load(target)['path'],path)
    assert not list(tmp_path.glob('*.tmp*'))
