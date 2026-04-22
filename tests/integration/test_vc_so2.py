#!/usr/bin/env python3
''' Program tests. '''
import unittest
from test_zaghop import ZagHopTest


class VcSo2Test(ZagHopTest):
    """ Trajectory tests."""


    def test_ldsh_isc(self):
        """ Test diabatic basis FSSH on pyrazine model."""
        self.run_traj()
        self.compare_energy()


if __name__ == '__main__':
    unittest.main(verbosity=2)
