"""CPU regression for run-boundary side effects; no native execution."""
import unittest

from py_scripts.single_rbc_repair.continuous_protocol import evolve_stage


class RunBoundaryRegression(unittest.TestCase):
    def test_observation_does_not_reenter_native_initialization(self):
        class Coordinator:
            def __init__(self):
                self.calls = []

            def run(self, steps, dt):
                # The pinned Simulation::run calls _execSplitters and clears
                # object forces on every invocation, including run(0).
                self.calls.append((steps, dt))

        coordinator = Coordinator()
        self.assertEqual(evolve_stage(coordinator, 5000, .001, 100), 5000)
        self.assertEqual(coordinator.calls, [(5000, .001)])


if __name__ == '__main__':
    unittest.main()
