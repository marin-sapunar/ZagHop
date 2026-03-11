#!/usr/bin/env python3
''' Program tests. '''
import os
import re
import shutil
import subprocess
import time
import unittest
import numpy as np


class ZagHopTest(unittest.TestCase):
    """ Trajectory tests."""

    @classmethod
    def setUpClass(cls):
        name = re.sub(r'Test$', '', cls.__name__)
        name = re.sub(r'(?<!^)(?=[A-Z])', '_', name).lower()
        cls.logdir = f"test_{name}"
        cls.cwd = os.getcwd()
        cls.idir = os.path.join(cls.cwd, name)
        cls.common = os.path.join(cls.idir, "common")
        os.makedirs(cls.logdir, exist_ok=True)

    def setUp(self):
        """ Set up directory for running a specific test.

        If the test's own subdirectory already exists from a previous run,
        rename it based on its modification time before creating a fresh one."""
        self.name = self.id().split(".")[-1]
        self.rundir = os.path.join(self.logdir, self.name)
        self.inpdir = os.path.join(self.idir, self.name)
        if os.path.isdir(self.rundir):
            old_time = time.gmtime(os.path.getmtime(self.rundir))
            old_time = time.strftime("%y.%m.%d.%H.%M.%S", old_time)
            os.rename(self.rundir, self.rundir + "_" + old_time)
        shutil.copytree(self.inpdir, self.rundir)
        shutil.copytree(self.common, self.rundir, dirs_exist_ok=True)
        if not os.path.isdir(os.path.join(self.inpdir, "Reference")):
            raise FileNotFoundError("Reference directory not found.")

    def tearDown(self):
        """ Return to initial directory in case of a failed test. """
        os.chdir(self.cwd)

    def run_traj(self):
        """ Run zaghop for the current test case. """
        os.chdir(self.rundir)
        with open("nad.stdout", "w") as out, open("nad.stderr", "w") as err:
            subprocess.call("zaghop", stdout=out, stderr=err)
        self.assertTrue(os.path.isdir("Results"), "Error termination.")
        os.chdir(self.cwd)

    def compare_energy(self):
        """ Compare full contents of the energy.dat file. """
        new = np.loadtxt(self.rundir + "/Results/energy.dat", comments="#")
        ref = np.loadtxt(self.rundir + "/Reference/energy.dat", comments="#")
        self.assertTrue(np.allclose(new, ref), "Difference in energy.dat file.")

    def stdout_contains(self, string):
        """ Assert that a string appears in the nad.stdout file. """
        with open(os.path.join(self.rundir, "nad.stdout")) as f:
            self.assertIn(string, f.read(), f"String not found in nad.stdout: {string!r}")


if __name__ == '__main__':
    unittest.main(verbosity=2)
