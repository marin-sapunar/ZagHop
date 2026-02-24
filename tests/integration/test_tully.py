#!/usr/bin/env python3
''' Program tests. '''
import unittest
from test_zaghop import ZagHopTest


class TullyTest(ZagHopTest):
    """ Trajectory tests."""

    def test_fssh(self):
        """ Test adiabatic FSSH on Tully model I."""
        self.run_traj()
        self.compare_energy()

    def test_ldsh(self):
        """ Test LD-FSSH on Tully model I."""
        self.run_traj()
        self.compare_energy()

    def test_lzsh(self):
        """ Test LZSH on Tully model I."""
        self.run_traj()
        self.compare_energy()




if __name__ == '__main__':
    unittest.main(verbosity=2)
