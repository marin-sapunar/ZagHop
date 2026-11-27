#!/usr/bin/env python3
''' Program tests. '''
import json
import os
import subprocess
import unittest
import numpy as np
from test_zaghop import ZagHopTest


class VcPyrazineTest(ZagHopTest):
    """ Trajectory tests."""

    def test_edrift(self):
        """ Test max_toten_d option.
        
        Tests that the trajectory is stopped when the total energy
        drift exceeds the threshold."""
        self.run_traj()
        self.compare_energy()
        self.stdout_contains("Drift in total energy too large.")

    def test_estep(self):
        """ Test max_toten_d_step option.
        
        Tests that the trajectory is stopped when the total energy
        change in a single step exceeds the threshold."""
        self.run_traj()
        self.compare_energy()
        self.stdout_contains("Change in total energy too large.")

    def test_ldsh(self):
        """ Test diabatic basis FSSH on pyrazine model.

        The trajectory is run up to 50 fs and restarted from its backup."""
        self.run_traj()
        self.restart_traj(100.0)
        self.compare_energy()
        self.compare_results("cwf.dat")

    def test_ldsh_diabatic(self):
        """ Test istate_basis = diabatic.

        Diabatic state 3 is projected onto the adiabatic states (weights 0.17, 0.82 for adiabatic
        states 2 and 3). The seed is chosen so that the less likely state 2 is drawn, and the
        wave function starts in it."""
        self.run_traj()
        self.stdout_contains("Initial state 3 of the diabatic basis projected onto the "
                             "adiabatic basis.")
        self.stdout_contains("Initial active state: 2.")
        self.compare_energy()
        self.compare_results("cwf.dat")


    def test_ldsh_external(self):
        """ test_ldsh with the potential computed by eval_point (qm keyword).

        Includes the restart, with the wave functions of the stored points from the backup."""
        self.run_traj()
        self.restart_traj(100.0)
        self.compare_energy()
        self.compare_results("cwf.dat")

    def test_eval_point(self):
        """ eval_point at the initial geometry of test_ldsh, and zagpy sample state reading it.

        The energies are those of the first point of test_ldsh. The diabatic state 3 has the
        weights 0.17 and 0.82 in adiabatic states 2 and 3 (see test_ldsh_diabatic)."""
        from zagpy.sample.state import read_point
        os.chdir(self.rundir)
        labels = np.loadtxt("geom", usecols=0, dtype=str)
        q = np.loadtxt("geom", usecols=2)
        request = {
            "protocol": "zaghop", "version": 1, "point_id": 1, "calculation": "new_point",
            "coordinates": "normal_modes", "atoms": labels.tolist(),
            "geometry": [[x] for x in q],
            "states": [{"multiplicity": 1, "nstate": 4}],
            "request": {
                "gradients": [{"multiplicity": 1, "state": 4}],
                "nacv": [{"multiplicity": 1, "states": [3, 4]}],
                "wave_function": True,
                "properties": ["diabatic_transformation", "diabatic_hamiltonian"]}}
        with open("qm_in.json", "w", encoding="utf-8") as f:
            json.dump(request, f)
        with open("eval_point.stdout", "w") as out, open("eval_point.stderr", "w") as err:
            status = subprocess.call(["eval_point", "vibronic_coupling"], stdout=out, stderr=err)
        os.chdir(self.cwd)
        self.assertEqual(status, 0, "eval_point failed.")
        with open(os.path.join(self.rundir, "qm_out.json"), encoding="utf-8") as f:
            res = json.load(f)

        en = np.array(res["states"][0]["energies"])
        ref = np.loadtxt(os.path.join(self.rundir, "Reference", "energy.dat"), comments="#")
        self.assertTrue(np.allclose(en, ref[0, 4:8], atol=1.0e-9, rtol=0.0), "Energies differ.")
        self.assertEqual(np.array(res["gradients"][0]["gradient"]).shape, (len(q), 1))
        self.assertEqual(np.array(res["nacv"][0]["vector"]).shape, (len(q), 1))
        c = np.array(res["wave_function"][0]["c"])
        adt = np.array(res["properties"]["diabatic_transformation"])
        h = np.array(res["properties"]["diabatic_hamiltonian"])
        self.assertTrue(np.array_equal(adt, c))
        self.assertTrue(np.allclose(c.T @ c, np.eye(4), atol=1.0e-12))
        self.assertTrue(np.allclose(c.T @ h @ c, np.diag(en), atol=1.0e-12))

        en_z, osc_z, adt_z = read_point(self.rundir)
        self.assertIsNone(osc_z)
        self.assertTrue(np.array_equal(en_z, en) and np.array_equal(adt_z, adt))
        self.assertTrue(np.array_equal(np.round(adt_z[2, 1:3]**2, 2), [0.17, 0.82]))


if __name__ == '__main__':
    unittest.main(verbosity=2)
