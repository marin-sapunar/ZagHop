#!/usr/bin/env python3
''' Program tests. '''
import unittest
from test_zaghop import ZagHopTest


class VcPyrazineTest(ZagHopTest):
    """ Trajectory tests."""

    def test_edrift(self):
        """ Test max_toten_d option.
        
        Tests that the trajectory is stopped when the total energy
        drift exceeds the threshold."""
        self.run_traj()
        self.compare_energy()

    def test_estep(self):
        """ Test max_toten_d_step option.
        
        Tests that the trajectory is stopped when the total energy
        change in a single step exceeds the threshold."""
        self.run_traj()
        self.compare_energy()

    def test_fssh(self):
        """ Test adiabatic basis FSSH on pyrazine model."""
        self.run_traj()
        self.compare_energy()


if __name__ == '__main__':
    unittest.main(verbosity=2)
