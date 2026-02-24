#!/usr/bin/env python3
''' Program tests. '''
import os
import unittest
from test_zaghop import ZagHopTest


class TullyTest(ZagHopTest):
    """ Trajectory tests."""

    def test_fssh(self):
        """ Test adiabatic FSSH on Tully model I."""
        self.run_traj()
        self.compare_energy()

    def test_ldsh(self):
        """ Test LD-FSSH on Tully model I."""
        self.run_traj()
        self.compare_energy()

    def test_lzsh(self):
        """ Test LZSH on Tully model I."""
        self.run_traj()
        self.compare_energy()

    @classmethod
    def setUpClass(cls):
        """ Make test logs directory.

        Creates the shared log directory if it does not exist. If it already
        exists it is reused, so that tests running in quick succession (each
        in their own process) do not conflict with each other."""
        cls.logdir = "test_tully"
        cls.cwd = os.getcwd()
        cls.idir = os.path.join(cls.cwd, "tully")
        cls.common = os.path.join(cls.idir, "common")
        os.makedirs(cls.logdir, exist_ok=True)


if __name__ == '__main__':
    unittest.main(verbosity=2)
