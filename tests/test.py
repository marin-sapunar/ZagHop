#!/bin/env python3
''' Program tests. '''
import os
import shutil
import time
import unittest
import subprocess
import numpy as np


class ZagHopTest(unittest.TestCase):


    def setup_test(self, name):
        self.idir = os.getcwd()
        shutil.copytree(name, "run_" + name)
        os.chdir("run_" + name)
        if not os.path.isdir("Reference"):
            raise FileNotFoundError("Reference directory not found.")


    def run_trajectory(self):
        with open("nad.stdout", "w") as out, open("nad.stderr", "w") as err:
            subprocess.call("zaghop", stdout=out, stderr=err)
        self.assertTrue(os.path.isdir("Results"), "Error termination.")


    def compare_energy(self):
        new = np.loadtxt("Results/energy.dat", comments="#")
        ref = np.loadtxt("Reference/energy.dat", comments="#")
        self.assertTrue(np.allclose(new, ref), "Difference in energy.dat file.")


    def end_test(self):
        os.chdir(self.idir)


class Turbomole(ZagHopTest):
    """ Turbomole interface tests. """

    def test_adc2_pyrazine(self):
        self.setup_test("ld-fssh_pyrazine")
        self.run_trajectory()
        self.compare_energy()
        self.end_test()


#def check_paths():
#    ''' Check that REFDIR and PROG exist. '''
#    if not os.path.isdir(REFDIR):
#        raise RuntimeError(REFDIR + ' directory not found.')
#    if not os.path.isfile(PROG[0]):
#        raise RuntimeError(PROG[0]+ ' not found.')
#    try:
#        os.mkdir(LOGDIR)
#    except FileExistsError:
#        old_log_time = time.gmtime(os.path.getmtime(LOGDIR))
#        old_log_time = time.strftime("%y.%m.%d._%H:%M:%S", old_log_time)
#        os.rename(LOGDIR, LOGDIR[:-1] + old_log_time)
#        os.mkdir(LOGDIR)


if __name__ == '__main__':
#   check_paths()
    unittest.main(verbosity=2)
