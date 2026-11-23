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
        self.stdout_contains("Drift in total energy too large.")

    def test_estep(self):
        """ Test max_toten_d_step option.
        
        Tests that the trajectory is stopped when the total energy
        change in a single step exceeds the threshold."""
        self.run_traj()
        self.compare_energy()
        self.stdout_contains("Change in total energy too large.")

    def test_ldsh(self):
        """ Test diabatic basis FSSH on pyrazine model."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("cwf.dat")

    def test_ldsh_diabatic(self):
        """ Test istate_basis = diabatic.

        Diabatic state 3 is projected onto the adiabatic states (weights 0.17, 0.82 for adiabatic
        states 2 and 3). The seed is chosen so that the less likely state 2 is drawn, and the
        wave function starts in it."""
        self.run_traj()
        self.stdout_contains("Initial state 3 of the diabatic basis projected onto the "
                             "adiabatic basis.")
        self.stdout_contains("Initial active state: 2.")
        self.compare_energy()
        self.compare_results("cwf.dat")


if __name__ == '__main__':
    unittest.main(verbosity=2)
