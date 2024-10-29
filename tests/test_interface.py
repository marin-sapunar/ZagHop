#!/bin/env python3
''' Program tests. '''
import os
import shutil
import subprocess
import time
import unittest
import numpy as np
import yaml
import run


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
