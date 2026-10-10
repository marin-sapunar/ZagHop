#!/usr/bin/env python3
''' Program tests. '''
import os
import unittest
from test_zaghop import ZagHopTest


class QuanticsSo2Test(ZagHopTest):
    """ Trajectory tests using the quantics operator library.

    The quantics interface reads the restart, dvr and oper files from two
    directories above the running trajectory, so the full directory
    structure under common is copied and zaghop is run in zagreb_trj/traj."""

    trajdir = os.path.join("zagreb_trj", "traj")

    def test_ldsh_isc(self):
        """ Test locally diabatic FSSH with SOC on SO2 model.

        Locally diabatic propagation does not depend on the ordering of the
        diabatic states, so the diabatic coefficients are also compared to
        check the reordering of the triplet components."""
        self.run_traj()
        self.compare_energy()
        self.compare_results("cstate_diab")


if __name__ == '__main__':
    unittest.main(verbosity=2)
