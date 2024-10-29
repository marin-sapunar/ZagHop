#!/bin/env python3
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

        If a directory from a previous test run exists, rename that directory
        based on the time it was modified and create a new one."""
        cls.logdir = "zaghop_test"
        cls.idir = os.getcwd()
        try:
            os.mkdir(cls.logdir)
        except FileExistsError:
            old_log_time = time.gmtime(os.path.getmtime(cls.logdir))
            old_log_time = time.strftime("%y.%m.%d.%H.%M.%S", old_log_time)
            os.rename(cls.logdir, cls.logdir + "_" + old_log_time)
            os.mkdir(cls.logdir)

    def setUp(self):
        """ Set up directory for running a specific test."""
        self.name = self.id().split(".")[-1][5:]
        self.rundir = self.logdir + "/" + self.name
        self.inpdir = self.idir + "/" + self.name
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


class InterfaceTest(unittest.TestCase):
    def test_turbomole_adc2(self):
        yaml = self.system + """
        options : 
          method : ricc2
          model : adc(2)
          basis : def2-SVP

        request :
          energy : True
          gradient : 2
        """
        ref_en = np.array((-56.3328740, -56.0564136, -55.9764839, -55.9660628))
        ref_gr = np.array((-0.02986388, 0.03397675, 0.02238019))
        with open("qm.yaml", "w") as infile:
            infile.write(yaml)
        args = self.parser.parse_args(["-init", "turbomole", "."])
        res = run.run(args)
        res_en = np.array(res["energy"])
        res_gr = np.array(res["gradient"])[0]
        self.assertTrue(np.allclose(ref_en, res_en), "Energies do not match")
        self.assertTrue(np.allclose(ref_gr, res_gr), "Gradients do not match")

        
    @classmethod
    def setUpClass(cls):
        """ Make test logs directory.

        If a directory from a previous test run exists, rename that directory
        based on the time it was modified and create a new one."""
        cls.logdir = "interface_test"
        cls.idir = os.getcwd()
        try:
            os.mkdir(cls.logdir)
        except FileExistsError:
            old_log_time = time.gmtime(os.path.getmtime(cls.logdir))
            old_log_time = time.strftime("%y.%m.%d.%H.%M.%S", old_log_time)
            os.rename(cls.logdir, cls.logdir + "_" + old_log_time)
            os.mkdir(cls.logdir)
        system ="""
        system :
          atoms : [ n, h, h, h ]
          geom :
            - [    -0.621219E-02 ,    0.399125E-02 ,    0.129238E+00 ]
            - [     0.182386E+01 ,    0.484298E-02 ,   -0.623031E+00 ]
            - [    -0.841656E+00 ,    0.149843E+01 ,   -0.566004E+00 ]
            - [    -0.895890E+00 ,   -0.155873E+01 ,   -0.606655E+00 ]
          states :
            - nstate : 4
        """
        cls.system = system
        cls.parser = run.cli()

    def setUp(self):
        """ Set up directory for running a specific test."""
        self.name = self.id().split(".")[-1][5:]
        self.rundir = self.logdir + "/" + self.name
        os.mkdir(self.rundir)
        os.chdir(self.rundir)

    def tearDown(self):
        """ Return to initial directory in case of a failed test. """
        os.chdir(self.idir)





if __name__ == '__main__':
    unittest.main(verbosity=2)
