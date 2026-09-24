from particle_3d.injection_population import PopulationSource,STREAMS
import numpy as np
def test_full_state(source):
 for i,name in enumerate(STREAMS): source.rng[name].random(i+2)
 clone=PopulationSource.restore(source.state())
 for name in STREAMS: assert np.array_equal(source.rng[name].random(31),clone.rng[name].random(31))
