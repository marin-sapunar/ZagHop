#!/usr/bin/env python3
''' Program tests. '''
import os
import unittest
from test_zaghop import ZagHopTest


class PyrazineTest(ZagHopTest):
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

    @classmethod
    def setUpClass(cls):
        """ Make test logs directory.

        Creates the shared log directory if it does not exist. If it already
        exists it is reused, so that tests running in quick succession (each
        in their own process) do not conflict with each other."""
        cls.logdir = "test_vc_pyrazine"
        cls.cwd = os.getcwd()
        cls.idir = os.path.join(cls.cwd, "vc_pyrazine")
        cls.common = os.path.join(cls.idir, "common")
        os.makedirs(cls.logdir, exist_ok=True)


if __name__ == '__main__':
    unittest.main(verbosity=2)
