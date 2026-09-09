"""CPU diagnostics; synthetic cases are never material-measurement evidence."""
import unittest
import numpy as np
from py_scripts.fluid_physics.analysis import block_statistics


class StatisticsRepairTests(unittest.TestCase):
    def test_mean_and_ci_use_same_complete_blocks(self):
        values=np.r_[np.tile([-1.,1.],50),[50.,50.,50.,50.]]
        result=block_statistics(np.arange(len(values)),values,min_duration=10,min_blocks=8)
        self.assertGreaterEqual(result['block_count'],2)
        self.assertGreater(result['discarded_tail_samples'],0)
        self.assertAlmostEqual(result['mean'],float(np.mean(result['block_means'])),places=12)


if __name__=='__main__':unittest.main()
