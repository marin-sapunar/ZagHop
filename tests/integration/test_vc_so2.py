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
        self.compare_results("prob")
        self.compare_results("cstate_diab")

    def test_ldsh_isc_cart(self):
        """ Test diabatic basis FSSH on SO2 model in Cartesian coordinates."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")
        self.compare_results("cstate_diab")

    def test_ldsh_isc_cart_rot(self):
        """ Test removal of translation/rotation (orient = 2) in Cartesian coordinates.

        The initial geometry of test_ldsh_isc_cart is rotated and translated, and a centre of
        mass velocity and a rotation about the reference geometry are added to the velocity.
        Aligning the initial conditions with the reference geometry of the normal modes and
        projecting translation/rotation from the velocity reproduces test_ldsh_isc_cart."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")
        self.compare_results("cstate_diab")
        self.compare_results("energy.dat", reference="test_ldsh_isc_cart")

    def test_ldsh_isc_eckart(self):
        """ Test FSSH on SO2 model in Cartesian coordinates in the Eckart frame.

        The potential is invariant to rotation, so the angular momentum is conserved (up to the
        precision of the printed geometry and velocity), including at hops."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")
        self.compare_results("cstate_diab")
        ang = self.angular_momentum()
        self.assertTrue(np.abs(ang - ang[0]).max() < 1e-3, "Angular momentum not conserved.")

    def test_ldsh_isc_eckart_rot(self):
        """ Test FSSH on SO2 model in the Eckart frame starting from a rotated geometry.

        The initial geometry of test_ldsh_isc_eckart is rotated and translated, and the velocity
        rotated. The potential is invariant to rotation and translation, which reproduces
        test_ldsh_isc_eckart."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")
        self.compare_results("cstate_diab")
        self.compare_results("energy.dat", reference="test_ldsh_isc_eckart")

    def test_ldsh_isc_eff_nm(self):
        """ Test LD-FSSH on effective spin states built from the SO2 model (projection).

        The model has one state for each of the 7 spin-free states, instead of 13 components."""
        self.run_traj()
        self.assertEqual(self.load_results("energy.dat").shape[1], 4 + 7, "Not 7 states.")
        self.stdout_contains("Smallest retained fraction of the singlet-triplet couplings")
        self.compare_energy()
        self.compare_results("prob")
        self.compare_results("cstate_diab")

    def test_ldsh_isc_eff_rss_nm(self):
        """ Test LD-FSSH on effective spin states built from the SO2 model (rss, signed)."""
        self.run_traj()
        self.assertEqual(self.load_results("energy.dat").shape[1], 4 + 7, "Not 7 states.")
        self.compare_energy()
        self.compare_results("prob")
        self.compare_results("cstate_diab")

    def test_lzsh_isc_eff_nm(self):
        """ Test LZSH on effective spin states built from the SO2 model (projection).

        The combination of the triplet components is computed at a given reference geometry.
        The spin-orbit couplings of the model are constant, so it is the same as at the default
        reference."""
        self.run_traj()
        self.stdout_contains("Reading reference geometry file ref_geom.")
        self.stdout_contains("ISC probability evaluated")
        self.compare_energy()


if __name__ == '__main__':
    unittest.main(verbosity=2)
