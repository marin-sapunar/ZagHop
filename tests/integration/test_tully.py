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

    def test_lzsh(self):
        """ Test LZSH on Tully model I."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("geometry")
        self.assert_regular_grid(0.5)

    def test_lzsh_bisect(self):
        """ Test LZSH with points inserted around the gap minimum on Tully model I.

        The hop happens at an inserted point, from where the trajectory returns to the output
        grid. Up to the hop, the output equals that of the old implementation."""
        self.run_traj()
        self.stdout_contains("Adding new substep")
        self.assert_regular_grid(2.0)
        self.compare_energy()
        self.compare_results("geometry")

    def test_lzsh_frustrated(self):
        """ Test a frustrated LZSH hop on Tully model I.

        The kinetic energy is enough to cross the barrier of the lower state, but not the gap at
        the minimum. After the hop is rejected, the same gap minimum is not evaluated again."""
        self.run_traj()
        self.assertEqual(self.stdout_count("Insufficient energy for hop."), 1)
        energy = self.load_results("energy.dat")
        self.assertTrue(np.all(energy[:, 1] == 1), "Frustrated hop changed the state.")
        self.assertGreater(self.load_results("geometry")[-1], 0.0, "Barrier not crossed.")
        self.compare_energy()

    def test_lzsh_frustrated_reverse(self):
        """ Test a frustrated LZSH hop with velocity reversal on Tully model I.

        The trajectory is reflected at the gap minimum, so the point after the minimum differs
        from the one at which the hop was evaluated."""
        self.run_traj()
        self.assertEqual(self.stdout_count("Insufficient energy for hop."), 1)
        energy = self.load_results("energy.dat")
        self.assertTrue(np.all(energy[:, 1] == 1), "Frustrated hop changed the state.")
        self.assertLess(self.load_results("geometry")[-1], 0.0, "Trajectory not reflected.")
        self.compare_energy()

    def test_lzsh_ii(self):
        """ Test LZSH on Tully model II, with hops at both avoided crossings.

        The output equals that of the old implementation (commit 27e1d1a) at all points of the
        output grid."""
        self.run_traj()
        self.assertEqual(self.stdout_count("Hop occurred"), 2)
        self.stdout_contains("Adding new substep: 2")
        self.assert_regular_grid(0.5)
        self.compare_energy()
        self.compare_results("geometry")


if __name__ == '__main__':
    unittest.main(verbosity=2)
