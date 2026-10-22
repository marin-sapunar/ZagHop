#!/usr/bin/env python3
''' Program tests. '''
import unittest
import numpy as np
from test_zaghop import ZagHopTest


class VcSo2Test(ZagHopTest):
    """ Trajectory tests."""


    def test_ldsh_isc_nm(self):
        """ Test diabatic basis FSSH on SO2 model."""
        self.run_traj()
        self.compare_energy()

    def test_ldsh_isc_cart(self):
        """ Test diabatic basis FSSH on SO2 model in Cartesian coordinates."""
        self.run_traj()
        self.compare_energy()

    def test_ldsh_isc_cart_rot(self):
        """ Test removal of translation/rotation (orient = 2) in Cartesian coordinates.

        The initial geometry of test_ldsh_isc_cart is rotated and translated, and a centre of
        mass velocity and a rotation about the reference geometry are added to the velocity.
        Aligning the initial conditions with the reference geometry of the normal modes and
        projecting translation/rotation from the velocity reproduces test_ldsh_isc_cart."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("energy.dat", reference="test_ldsh_isc_cart")

    def test_ldsh_isc_eckart(self):
        """ Test FSSH on SO2 model in Cartesian coordinates in the Eckart frame.

        The potential is invariant to rotation, so the angular momentum is conserved (up to the
        precision of the printed geometry and velocity), including at hops."""
        self.run_traj()
        self.compare_energy()
        ang = self.angular_momentum()
        self.assertTrue(np.abs(ang - ang[0]).max() < 1e-3, "Angular momentum not conserved.")

    def test_ldsh_isc_eckart_rot(self):
        """ Test FSSH on SO2 model in the Eckart frame starting from a rotated geometry.

        The initial geometry of test_ldsh_isc_eckart is rotated and translated, and the velocity
        rotated. The potential is invariant to rotation and translation, which reproduces
        test_ldsh_isc_eckart."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("energy.dat", reference="test_ldsh_isc_eckart")


if __name__ == '__main__':
    unittest.main(verbosity=2)
