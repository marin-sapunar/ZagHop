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

    # Subdirectory of the run directory in which zaghop is executed.
    trajdir = ""

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
        rename it based on the modification time of the directory in which
        zaghop was run before creating a fresh one. (copytree copies the
        modification time of the source directory, so the rundir itself only
        gets a new one if zaghop writes to it directly.)"""
        self.name = self.id().split(".")[-1]
        self.rundir = os.path.join(self.logdir, self.name)
        self.inpdir = os.path.join(self.idir, self.name)
        if os.path.isdir(self.rundir):
            old_time = time.gmtime(os.path.getmtime(os.path.join(self.rundir, self.trajdir)))
            old_time = time.strftime("%y.%m.%d.%H.%M.%S", old_time)
            backup = self.rundir + "_" + old_time
            i = 1
            while os.path.exists(backup):
                backup = f"{self.rundir}_{old_time}_{i}"
                i += 1
            os.rename(self.rundir, backup)
        shutil.copytree(self.inpdir, self.rundir)
        shutil.copytree(self.common, self.rundir, dirs_exist_ok=True)
        if not os.path.isdir(os.path.join(self.inpdir, "Reference")):
            raise FileNotFoundError("Reference directory not found.")

    def tearDown(self):
        """ Return to initial directory in case of a failed test. """
        os.chdir(self.cwd)

    def run_traj(self):
        """ Run zaghop for the current test case. """
        os.chdir(os.path.join(self.rundir, self.trajdir))
        with open("nad.stdout", "w") as out, open("nad.stderr", "w") as err:
            subprocess.call("zaghop", stdout=out, stderr=err)
        self.assertTrue(os.path.isdir("Results"), "Error termination.")
        os.chdir(self.cwd)

    def compare_energy(self):
        """ Compare full contents of the energy.dat file. """
        self.compare_results("energy.dat")

    def compare_results(self, fname, reference=None):
        """ Compare full contents of a file in the Results directory.

        By default, the file is compared with the Reference directory of the current test. If
        reference is given, it is compared with the Reference directory of that test instead."""
        new = np.loadtxt(os.path.join(self.rundir, self.trajdir, "Results", fname), comments="#")
        if reference is None:
            refdir = os.path.join(self.rundir, "Reference")
        else:
            refdir = os.path.join(self.idir, reference, "Reference")
        ref = np.loadtxt(os.path.join(refdir, fname), comments="#")
        self.assertTrue(np.allclose(new, ref), f"Difference in {fname} file.")

    def angular_momentum(self):
        """ Angular momentum (atomic units) about the centre of mass at each printed step.

        Requires Cartesian coordinates and the geometry and velocity outputs (print 3,4). The
        masses are read from the geom file in the run directory."""
        rdir = os.path.join(self.rundir, self.trajdir)
        mass = np.loadtxt(os.path.join(rdir, "geom"), usecols=1, ndmin=1) * 1822.88849
        natom = len(mass)
        geom = np.loadtxt(os.path.join(rdir, "Results", "geometry")).reshape(-1, natom, 3)
        velo = np.loadtxt(os.path.join(rdir, "Results", "velocity")).reshape(-1, natom, 3)
        com = np.einsum("i,tia->ta", mass, geom) / mass.sum()
        return np.einsum("i,tia->ta", mass, np.cross(geom - com[:, None, :], velo))

    def stdout_contains(self, string):
        """ Assert that a string appears in the nad.stdout file. """
        with open(os.path.join(self.rundir, self.trajdir, "nad.stdout")) as f:
            self.assertIn(string, f.read(), f"String not found in nad.stdout: {string!r}")

    def stdout_count(self, string):
        """ Number of times a string appears in the nad.stdout file. """
        with open(os.path.join(self.rundir, self.trajdir, "nad.stdout")) as f:
            return f.read().count(string)

    def load_results(self, fname):
        """ Load a file from the Results directory. """
        return np.loadtxt(os.path.join(self.rundir, self.trajdir, "Results", fname), comments="#")

    def assert_regular_grid(self, dt):
        """ Assert that the times in energy.dat are exactly the output grid with step dt (fs). """
        time = self.load_results("energy.dat")[:, 0]
        grid = np.round(np.arange(len(time)) * dt, 5)
        self.assertTrue(np.array_equal(time, grid), "Times in energy.dat are not the output grid.")


if __name__ == '__main__':
    unittest.main(verbosity=2)
