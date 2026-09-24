import numpy as np
from particle_3d.injection_population import PopulationSource
def test_restore(source):
 source.orientation(13); clone=PopulationSource.restore(source.state())
 assert np.array_equal(source.orientation(100),clone.orientation(100))
