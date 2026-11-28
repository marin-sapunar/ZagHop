#!/usr/bin/env python3
''' Program tests. '''
import unittest
import numpy as np
from test_zaghop import ZagHopTest


class VcSo2Test(ZagHopTest):
    """ Trajectory tests."""


    def test_ldsh_isc_nm(self):
        """ Test diabatic basis FSSH on SO2 model, with the populations from the coefficients
        (print 41-43).

        The populations of each basis (13 states with their spin components) sum to 1, the
        adiabatic ones are |c|**2 of cwf.dat and the diabatic ones at t = 0 those of the initial
        active state (cstate_diab). The MCH states mix only diabatic states of the same
        multiplicity and spin component, so the singlets and each component of the triplets
        have the same total population in the MCH and diabatic bases. states.json gives the
        number of states of each basis."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")
        self.compare_results("cstate_diab")
        n_row = len(self.load_results("energy.dat"))
        qpop = {b: self.load_results(f"qpop_{b}") for b in ("adiab", "diab", "mch")}
        for b in ("adiab", "diab", "mch"):
            self.assertEqual(qpop[b].shape, (n_row, 13), f"Shape of qpop_{b}.")
            self.assertTrue(np.allclose(qpop[b].sum(axis=1), 1.0), f"qpop_{b} not normalised.")
        cwf = self.load_results("cwf.dat")
        self.assertTrue(np.allclose(qpop["adiab"], cwf[:, :13]**2 + cwf[:, 13:]**2),
                        "Adiabatic populations differ from |c|**2.")
        cstate = self.load_results("cstate_diab")
        self.assertTrue(np.allclose(qpop["diab"][0], cstate[0]**2),
                        "Initial diabatic populations differ from the active state.")
        for first in (0, 4, 7, 10):
            n = 4 if first == 0 else 3
            self.assertTrue(np.allclose(qpop["mch"][:, first:first+n].sum(axis=1),
                                        qpop["diab"][:, first:first+n].sum(axis=1)),
                            f"Populations of the group of states from {first + 1} differ between "
                            "the MCH and diabatic bases.")
        self.assertFalse(np.allclose(qpop["mch"], qpop["diab"]), "MCH = diabatic populations.")
        for b in ("adiab", "diab", "mch"):
            self.compare_results(f"qpop_{b}")
        states = self.load_states()
        self.assertEqual(states["format_version"], 3)
        self.assertEqual(states["dynamics_basis"], "adiabatic")
        self.assertEqual(states["n_state"], {"diabatic": 13, "adiabatic": 13, "mch": 13,
                                             "diabatic_sf": 7, "mch_sf": 7})

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

    def test_ldsh_isc_nstate(self):
        """ Test LD-FSSH on a subset of the states of the SO2 model (3 of 4 singlets, 2 of 3
        triplets).

        The adiabatic states are the eigenstates of the Hamiltonian in the selected MCH states,
        9 states. The diabatic basis keeps the 13 states of the model, and the active state is
        normalised in it. The total energy is conserved as well as with all states
        (test_ldsh_isc_nm)."""
        self.run_traj()
        self.stdout_contains("Multiplicity 1: 3 of 4 states.")
        self.stdout_contains("Multiplicity 3: 2 of 3 states.")
        energy = self.load_results("energy.dat")
        self.assertEqual(energy.shape[1], 4 + 9, "Not 9 states.")
        self.assertEqual(self.load_results("prob").shape[1], 9, "Not 9 probabilities.")
        cstate = self.load_results("cstate_diab")
        self.assertEqual(cstate.shape[1], 13, "Not 13 diabatic states.")
        self.assertTrue(np.allclose((cstate**2).sum(axis=1), 1.0), "Active state not normalised.")
        drift = np.abs(energy[:, 2] - energy[0, 2]).max()
        reference = self.load_reference("test_ldsh_isc_nm", "energy.dat")
        self.assertLess(drift, np.abs(reference[:, 2] - reference[0, 2]).max(),
                        "Total energy drift larger than with all states.")
        self.compare_energy()
        self.compare_results("prob")
        self.compare_results("cstate_diab")

    def test_lzsh_isc_eff_nstate(self):
        """ Test LZSH on effective spin states with a subset of the states of the SO2 model.

        The MCH states of a subset are those of the model, and the trajectory of
        test_lzsh_isc_eff_nm stays in the selected states (3 singlets and 2 triplets), so its
        reference is reproduced. The trajectory is restarted from its backup at 59 fs."""
        self.run_traj()
        self.restart_traj(100.0)
        self.stdout_contains("ISC probability evaluated")
        energy = self.load_results("energy.dat")
        self.assertEqual(energy.shape[1], 4 + 5, "Not 5 states.")
        reference = self.load_reference("test_lzsh_isc_eff_nm", "energy.dat")
        self.assertTrue(np.allclose(energy, reference[:, [0, 1, 2, 3, 4, 5, 6, 8, 9]]),
                        "Difference from the energies of the full model.")
        self.compare_energy()

    def test_lzsh_isc_eff_nm(self):
        """ Test LZSH on effective spin states built from the SO2 model (projection).

        The combination of the triplet components is computed at a given reference geometry.
        The spin-orbit couplings of the model are constant, so it is the same as at the default
        reference. The trajectory is run up to 59 fs and restarted from its backup, shortly
        before a hop. LZSH runs in the mch_sf basis, which is the mch basis for effective spin
        states, and has no coefficients for the populations of print option 43."""
        self.run_traj()
        self.input_error(" noprint  2,3,4,5,6", " noprint  2,3,4,5,6\n print 43",
                         "needs the coefficients of the electronic wave function")
        self.restart_traj(100.0)
        self.stdout_contains("Reading reference geometry file ref_geom.")
        self.stdout_contains("ISC probability evaluated")
        self.compare_energy()
        states = self.load_states()
        self.assertEqual(states["dynamics_basis"], "mch_sf")
        self.assertEqual(states["n_state"], {"diabatic": 7, "adiabatic": 7, "mch": 7,
                                             "diabatic_sf": 7, "mch_sf": 7})


if __name__ == '__main__':
    unittest.main(verbosity=2)
