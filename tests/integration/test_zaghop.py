#!/usr/bin/env python3
''' Program tests. '''
import os
import shutil
import subprocess
import time
import unittest
import numpy as np
import yaml


class ZagHopTest(unittest.TestCase):
    """ Trajectory tests."""

    def test_ldsh_tullyI(self):
        """ Tully model I test case for LD-FSSH."""
        self.run_traj()
        self.compare_energy()

    def test_lzsh_tullyI(self):
        """ Tully model I test case for LZSH."""
        self.run_traj()
        self.compare_energy()

    def test_fssh_tullyI(self):
        """ Tully model I test case for A-FSSH."""
        self.run_traj()
        self.compare_energy()

    def test_fssh_vc_pyrazine(self):
        """ Pyrazine model test case for A-FSSH."""
        self.run_traj()
        self.compare_energy()

#   def test_gsmd_nh3(self):
#       """ Ammonia test case for ground state MD with thermostat.

#       Tests the thermostat and that options are correctly set when surface
#       hopping is turned off.
#       """
#       self.run_traj()
#       self.compare_energy()

    @classmethod
    def setUpClass(cls):
        """ Make test logs directory.

        Creates the shared log directory if it does not exist. If it already
        exists it is reused, so that tests running in quick succession (each
        in their own process) do not conflict with each other."""
        cls.logdir = "zaghop_test"
        cls.idir = os.getcwd()
        os.makedirs(cls.logdir, exist_ok=True)

    def setUp(self):
        """ Set up directory for running a specific test.

        If the test's own subdirectory already exists from a previous run,
        rename it based on its modification time before creating a fresh one."""
        self.name = self.id().split(".")[-1][5:]
        self.rundir = self.logdir + "/" + self.name
        self.inpdir = self.idir + "/" + self.name
        if os.path.isdir(self.rundir):
            old_time = time.gmtime(os.path.getmtime(self.rundir))
            old_time = time.strftime("%y.%m.%d.%H.%M.%S", old_time)
            os.rename(self.rundir, self.rundir + "_" + old_time)
        shutil.copytree(self.inpdir, self.rundir)
        if not os.path.isdir(self.inpdir + "/Reference"):
            raise FileNotFoundError("Reference directory not found.")

    def tearDown(self):
        """ Return to initial directory in case of a failed test. """
        os.chdir(self.idir)

    def run_traj(self):
        """ Run zaghop for the current test case. """
        os.chdir(self.rundir)
        with open("nad.stdout", "w") as out, open("nad.stderr", "w") as err:
            subprocess.call("zaghop", stdout=out, stderr=err)
        self.assertTrue(os.path.isdir("Results"), "Error termination.")
        os.chdir(self.idir)

    def compare_energy(self):
        """ Compare full contents of the energy.dat file. """
        new = np.loadtxt(self.rundir + "/Results/energy.dat", comments="#")
        ref = np.loadtxt(self.rundir + "/Reference/energy.dat", comments="#")
        self.assertTrue(np.allclose(new, ref), "Difference in energy.dat file.")


if __name__ == '__main__':
    unittest.main(verbosity=2)
