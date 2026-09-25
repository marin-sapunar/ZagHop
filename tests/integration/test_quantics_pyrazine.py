#!/usr/bin/env python3
''' Program tests. '''
import os
import unittest
from test_zaghop import ZagHopTest


class QuanticsPyrazineTest(ZagHopTest):
    """ Trajectory tests using the quantics operator library.

    The quantics interface reads the restart, dvr and oper files from two
    directories above the running trajectory, so the full directory
    structure under common is copied and zaghop is run in zagreb_trj/traj."""

    trajdir = os.path.join("zagreb_trj", "traj")

    def test_ldsh(self):
        """ Test diabatic basis FSSH on pyrazine model."""
        self.run_traj()
        self.compare_energy()


if __name__ == '__main__':
    unittest.main(verbosity=2)
