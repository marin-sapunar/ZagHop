#!/usr/bin/env python3
''' Program tests. '''
import os
import unittest
import numpy as np
from test_zaghop import ZagHopTest


class TullyTest(ZagHopTest):
    """ Trajectory tests."""

#   def test_fssh(self):
#       """ Test adiabatic FSSH on Tully model I."""
#       self.run_traj()
#       self.compare_energy()

    def test_ldsh(self):
        """ Test LD-FSSH on Tully model I."""
        self.run_traj()
        self.compare_energy()

    def test_surfhop_off(self):
        """ Test that the trajectory stays in the initial state with $surfhop off.

        With the same seed, LD-FSSH with the default options hops to state 2 at step 98."""
        self.run_traj()
        self.compare_energy()
        energy = np.loadtxt(os.path.join(self.rundir, self.trajdir, "Results", "energy.dat"))
        self.assertTrue(np.all(energy[:, 1] == 1), "Hop with surface hopping turned off.")

#   def test_lzsh(self):
#       """ Test LZSH on Tully model I."""
#       self.run_traj()
#       self.compare_energy()




if __name__ == '__main__':
    unittest.main(verbosity=2)
