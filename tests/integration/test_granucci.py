#!/usr/bin/env python3
''' Program tests. '''
import unittest
import numpy as np
from test_zaghop import ZagHopTest


class GranucciTest(ZagHopTest):
    """ Trajectory tests on the one-dimensional models of Granucci, Persico and Spighi,
    J. Chem. Phys. 137, 22A501 (2012).

    granucci-st has one singlet and one triplet, with couplings sigma(r) h that change sign at
    r_sigma. The couplings are of rank 1, so the projection and signed rss models are exact.
    granucci-sst has two singlets, whose couplings to the triplet are at the angle phi. Each
    trajectory starts at r = 9 bohr, before the singlet-triplet crossing at r_c = 10 bohr. The
    velocity is rescaled along the velocity, because the coupling vectors between the bright and
    the dark triplet components are exactly zero."""

    def max_prob(self):
        """ Largest hopping probability of the trajectory. """
        return self.load_results("prob").max()

    def test_ldsh_st_components(self):
        """ Test LD-FSSH on granucci-st with all components, with the sign change at r_c.

        The sign change at the crossing makes it a trivial crossing: the active state switches
        to the top state with probability 1, and keeps its singlet character."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")
        self.assertGreater(self.max_prob(), 0.999)
        self.assertEqual(self.load_results("energy.dat")[-1, 1], 4)

    def test_ldsh_st_projection(self):
        """ Test LD-FSSH on the projection model of granucci-st.

        The couplings are of rank 1, so the trajectory is that of all components."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")
        self.assertEqual(self.load_results("energy.dat").shape[1], 4 + 2, "Not 2 states.")
        new = self.load_results("energy.dat")[:, 2:4]
        ref = self.load_reference("test_ldsh_st_components", "energy.dat")[:, 2:4]
        self.assertTrue(np.allclose(new, ref), "Energies differ from all components.")
        self.assertGreater(self.max_prob(), 0.999)

    def test_ldsh_st_rss(self):
        """ Test LD-FSSH on the signed rss model of granucci-st.

        With signs by continuity, the coupling follows the sign of sigma(r), and the model is
        the projection model."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")
        self.compare_results("energy.dat", reference="test_ldsh_st_projection")
        self.compare_results("prob", reference="test_ldsh_st_projection")

    def test_ldsh_st_rss_unsigned(self):
        """ Test LD-FSSH on the unsigned rss model of granucci-st.

        The coupling |sigma(r)| |h| has a cusp at the crossing, so the crossing is no longer
        trivial: the probability of staying in the lower state, the triplet after the crossing,
        is about 0.16 (the fixed phase error of Granucci 2012)."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")
        self.assertLess(self.max_prob(), 0.9)

    def test_lzsh_st_components(self):
        """ Test LZSH on granucci-st with all components, with an S -> T hop.

        With r_sigma = 8, the coupling is constant around the crossing."""
        self.run_traj()
        self.stdout_contains("Intersystem crossing from state 1 to state 2")
        self.compare_energy()

    def test_lzsh_st_projection(self):
        """ Test LZSH on the projection model of granucci-st, which is that of all components."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("energy.dat", reference="test_lzsh_st_components")

    def test_ldsh_sst_components(self):
        """ Test LD-FSSH on granucci-sst with all components, at phi = 45 degrees, from S2."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")

    def test_ldsh_sst_projection(self):
        """ Test LD-FSSH on the projection model of granucci-sst at phi = 45 degrees.

        The couplings are of rank 2, and the model keeps (1 + cos(phi)) / 2 of them."""
        self.run_traj()
        self.stdout_contains("Retained fraction of the singlet-triplet couplings at the "
                             "reference:  0.853553")
        self.compare_energy()
        self.compare_results("prob")

    def test_ldsh_sst_rank1_components(self):
        """ Test LD-FSSH on granucci-sst with all components, at phi = 0, from S2.

        The bright triplet crosses the two dark components, which are exactly degenerate, in a
        trivial crossing."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")

    def test_ldsh_sst_rank1_rss(self):
        """ Test LD-FSSH on the rss model of granucci-sst at phi = 0.

        The couplings are of rank 1, so the energies are those of all components, up to the
        rescaling of the velocity at the trivial crossing with the dark components."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("prob")
        new = self.load_results("energy.dat")[:, 2:4]
        ref = self.load_reference("test_ldsh_sst_rank1_components", "energy.dat")[:, 2:4]
        self.assertTrue(np.allclose(new, ref), "Energies differ from all components.")


if __name__ == '__main__':
    unittest.main(verbosity=2)
