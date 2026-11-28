#!/usr/bin/env python3
''' Tests of zagpy commands: --seed of the commands that draw random numbers, and the estimators
of ensemble population. '''
import filecmp
import json
import os
import re
import shutil
import subprocess
import sys
import unittest
import numpy as np
from zagpy.normalmode import NormalModes
from zagpy.sample.wigner import thermal_wigner_q, thermal_wigner_v
from zagpy.units import boltzmann_au


# Geometry, frequencies and normal modes of the SO2 model.
V0 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vc_so2", "common", "V0.txt")


class ZagpyTest(unittest.TestCase):
    """ Reproducibility of zagpy commands with --seed, and the estimators of ensemble
    population. The input is generated in the tests. """

    def setUp(self):
        """ Create an empty run directory for the test (removing the one of a previous run). """
        self.cwd = os.getcwd()
        self.name = self.id().split(".")[-1]
        self.rundir = os.path.join(self.cwd, "test_zagpy", self.name)
        if os.path.isdir(self.rundir):
            shutil.rmtree(self.rundir)
        os.makedirs(self.rundir)
        os.chdir(self.rundir)

    def tearDown(self):
        """ Return to initial directory. """
        os.chdir(self.cwd)

    def zagpy(self, *args, cwd=".", fail=False):
        """ Run zagpy in cwd and return the printed seed (None if no seed was printed).

        The standard output and error are written to zagpy.stdout and zagpy.stderr in cwd. With
        fail=True, zagpy is expected to stop with an error."""
        os.makedirs(cwd, exist_ok=True)
        res = subprocess.run([sys.executable, "-m", "zagpy", *args], cwd=cwd,
                             capture_output=True, text=True, check=False)
        with open(os.path.join(cwd, "zagpy.stdout"), "w") as out:
            out.write(res.stdout)
        with open(os.path.join(cwd, "zagpy.stderr"), "w") as err:
            err.write(res.stderr)
        if fail:
            self.assertNotEqual(res.returncode, 0, "zagpy did not stop with an error.")
            return None
        self.assertEqual(res.returncode, 0, f"zagpy failed:\n{res.stderr}")
        seed = re.findall(r"^Random seed: (\d+)$", res.stdout, flags=re.MULTILINE)
        self.assertLessEqual(len(seed), 1, "More than one seed printed.")
        return int(seed[0]) if seed else None

    def wigner(self, run, *args):
        """ Sample Wigner points of the SO2 model in normal mode coordinates in directory run.

        Returns the printed seed, the displacements and the velocities (n_point, n_mode)."""
        seed = self.zagpy("sample", "wigner", "-f", "sharc_lvc", "-fn", V0,
                          "--coord", "normal_mode", *args, cwd=run)
        points = sorted(d for d in os.listdir(run) if d.startswith("point"))
        q = np.array([np.loadtxt(os.path.join(run, p, "geom"), usecols=2) for p in points])
        v = np.array([np.loadtxt(os.path.join(run, p, "veloc")) for p in points])
        return seed, q, v

    def assert_same_files(self, dir1, dir2, files):
        """ Assert that the files (paths relative to dir1 and dir2) are identical. """
        for f in files:
            self.assertTrue(filecmp.cmp(os.path.join(dir1, f), os.path.join(dir2, f),
                                        shallow=False), f"{dir1}/{f} and {dir2}/{f} differ.")

    def test_wigner_seed(self):
        """ The same seed gives the same points, a different seed different points, and the
        printed seed of a run without --seed repeats it. """
        seed1, q1, v1 = self.wigner("seed1", "-n", "5", "--seed", "11")
        _, q2, v2 = self.wigner("seed1_again", "-n", "5", "--seed", "11")
        _, q3, v3 = self.wigner("seed2", "-n", "5", "--seed", "12")
        self.assertEqual(seed1, 11)
        files = [os.path.join(f"point{i}", f) for i in range(1, 6) for f in ("geom", "veloc")]
        self.assert_same_files("seed1", "seed1_again", files)
        self.assertTrue(np.array_equal(q1, q2) and np.array_equal(v1, v2))
        self.assertFalse(np.any(q1 == q3) or np.any(v1 == v3), "Different seeds, same points.")

        seed, q4, v4 = self.wigner("no_seed", "-n", "5")
        self.assertIsNotNone(seed, "No seed printed.")
        _, q5, v5 = self.wigner("no_seed_again", "-n", "5", "--seed", str(seed))
        self.assert_same_files("no_seed", "no_seed_again", files)
        self.assertTrue(np.array_equal(q4, q5) and np.array_equal(v4, v5))

    def test_wigner_npoint(self):
        """ The first 3 points of a sample of 5 are the sample of 3 with the same seed. """
        self.wigner("n3", "-n", "3", "--seed", "7")
        self.wigner("n5", "-n", "5", "--seed", "7")
        files = [os.path.join(f"point{i}", f) for i in range(1, 4) for f in ("geom", "veloc")]
        self.assert_same_files("n3", "n5", files)

    def test_wigner_modes(self):
        """ The modes are sampled independently: the displacements and velocities divided by
        their widths differ between the modes (one generator for all modes, not a generator
        with the same seed per mode). """
        T = 300
        _, q, v = self.wigner("run", "-n", "4", "-T", str(T), "--seed", "3")
        freq = NormalModes.from_sharc_lvc(V0).freq
        T_au = T * boltzmann_au
        z_q = q / np.array([thermal_wigner_q(o, T_au) * np.sqrt(o) for o in freq])
        z_v = v / np.array([thermal_wigner_v(o, T_au) * np.sqrt(o) for o in freq])
        self.assertEqual(q.shape, (4, len(freq)))
        for i in range(len(freq)):
            for j in range(i):
                self.assertFalse(np.allclose(z_q[:, i], z_q[:, j]), f"q: modes {i}, {j}.")
                self.assertFalse(np.allclose(z_v[:, i], z_v[:, j]), f"v: modes {i}, {j}.")
            self.assertFalse(np.allclose(z_q[:, i], z_v[:, i]), f"q and v of mode {i}.")

    def write_points(self, n_point=12, n_state=3):
        """ Write qm_out.json files of eval_point in directories p01, p02, ... with oscillator
        strengths 0.5, so every state of every point is selected with probability 0.5. """
        dirs = [f"p{i + 1:02d}" for i in range(n_point)]
        for i, d in enumerate(dirs):
            os.makedirs(d)
            data = {"states": [{"multiplicity": 1,
                                "energies": [0.01 * i + 0.1 * s for s in range(n_state)]}],
                    "properties": {"oscillator_strengths": [0.5] * n_state}}
            with open(os.path.join(d, "qm_out.json"), "w", encoding="utf-8") as f:
                json.dump(data, f)
        return dirs

    def select(self, target, dirs, *args):
        """ Run sample state into directory target and return the seed and the selected
        point/state directories. """
        seed = self.zagpy("sample", "state", "--weight-oscill", "--weight-norm", "none",
                          "--target-dir", target, *args, *dirs)
        return seed, sorted(os.listdir(target))

    def test_state_seed(self):
        """ The same seed gives the same selection, a different seed a different one, and the
        printed seed of a run without --seed repeats it. """
        dirs = self.write_points()
        seed1, sel1 = self.select("seed1", dirs, "--seed", "5")
        _, sel2 = self.select("seed1_again", dirs, "--seed", "5")
        _, sel3 = self.select("seed2", dirs, "--seed", "6")
        self.assertEqual(seed1, 5)
        self.assertTrue(0 < len(sel1) < 36, "Expected a partial selection.")
        self.assertEqual(sel1, sel2)
        self.assertNotEqual(sel1, sel3)
        seed, sel4 = self.select("no_seed", dirs)
        self.assertIsNotNone(seed, "No seed printed.")
        _, sel5 = self.select("no_seed_again", dirs, "--seed", str(seed))
        self.assertEqual(sel4, sel5)

    def test_state_dir_order(self):
        """ The selection depends on the set of directories, not on their order. """
        dirs = self.write_points()
        _, sel1 = self.select("sorted", dirs, "--seed", "5")
        _, sel2 = self.select("reversed", dirs[::-1], "--seed", "5")
        self.assertEqual(sel1, sel2)

    def write_trajectories(self, n_traj=10, n_step=6, layout=True):
        """ Write Results/energy.dat (time, active state), Results/status (1) and, with layout,
        Results/states.json (dynamics in the adiabatic basis) of two-state trajectories t01,
        t02, ... with random active states. """
        state = np.random.default_rng(0).integers(1, 3, size=(n_traj, n_step))
        state[0] = 1
        state[1] = 2
        dirs = [f"t{i + 1:02d}" for i in range(n_traj)]
        for d, s in zip(dirs, state):
            self.write_trajectory(d, s, 1, layout=("adiabatic", 2) if layout else None)
        return dirs

    @staticmethod
    def write_trajectory(cdir, state, status, layout=None, files=None, version=3):
        """ Write the Results directory of a trajectory with active states state, printed every
        0.5 fs: energy.dat (time, active state), status, and, if given, states.json with the
        basis of the dynamics and its number of states (layout), in format version, and other
        files (name: array of rows). """
        os.makedirs(os.path.join(cdir, "Results"))
        np.savetxt(os.path.join(cdir, "Results", "energy.dat"),
                   np.column_stack([np.arange(len(state)) * 0.5, state]), fmt="%.2f %d",
                   header="time state")
        with open(os.path.join(cdir, "Results", "status"), "w") as f:
            f.write(f"{status}\n")
        if layout is not None:
            basis, n = layout
            with open(os.path.join(cdir, "Results", "states.json"), "w") as f:
                n_state = {"diabatic": 4, "adiabatic": 4, "mch": 4, "diabatic_sf": 2,
                           "mch_sf": 2}
                if version <= 2:
                    n_state = {"diabatic": 4, "adiabatic": 4, "mch": 2}
                n_state[basis] = n
                json.dump({"format_version": version, "dynamics_basis": basis,
                           "n_state": n_state}, f)
        for fname, rows in (files or {}).items():
            np.savetxt(os.path.join(cdir, "Results", fname), rows)

    def population(self, out, dirs, *args):
        """ Run ensemble population with --bootstrap into directory out. Returns the seed and the
        lower and upper bounds. """
        os.makedirs(out)
        seed = self.zagpy("ensemble", "population", "--bootstrap", "--output",
                          os.path.join(out, "pop"), *args, *dirs)
        bounds = [np.loadtxt(os.path.join(out, f"pop_active_adiabatic_{b}.dat"))
                  for b in ("lbound", "ubound")]
        return seed, bounds

    def test_population_seed(self):
        """ The same seed gives the same bootstrap intervals, a different seed different ones,
        and the printed seed of a run without --seed repeats it. Without --bootstrap, --seed is
        accepted and no seed is printed. """
        dirs = self.write_trajectories()
        files = ["pop_active_adiabatic_lbound.dat", "pop_active_adiabatic_ubound.dat"]
        seed1, b1 = self.population("seed1", dirs, "--seed", "9")
        self.population("seed1_again", dirs, "--seed", "9")
        _, b3 = self.population("seed2", dirs, "--seed", "10")
        self.assertEqual(seed1, 9)
        self.assert_same_files("seed1", "seed1_again", files)
        self.assertFalse(all(np.array_equal(x, y, equal_nan=True) for x, y in zip(b1, b3)),
                         "Different seeds, same intervals.")
        seed, _ = self.population("no_seed", dirs)
        self.assertIsNotNone(seed, "No seed printed.")
        self.population("no_seed_again", dirs, "--seed", str(seed))
        self.assert_same_files("no_seed", "no_seed_again", files)

        os.makedirs("no_bootstrap")
        seed = self.zagpy("ensemble", "population", "--seed", "9", "--output",
                          os.path.join("no_bootstrap", "pop"), *dirs)
        self.assertIsNone(seed, "Seed printed without --bootstrap.")

    def population_file(self, *args, name):
        """ Run ensemble population with args and load the population file name. """
        self.zagpy("ensemble", "population", "--output", "pop", *args)
        return np.loadtxt(name)

    def test_population_active(self):
        """ The 'active' estimator gives the fraction of trajectories in each active state, in
        the basis of the dynamics from states.json, with its number of states. Without
        states.json the values are the same, with a warning and the basis 'dynamics'. """
        dirs = self.write_trajectories()
        pop = self.population_file(*dirs, name="pop_active_adiabatic.dat")
        state = np.array([np.loadtxt(os.path.join(d, "Results", "energy.dat"))[:, 1]
                          for d in dirs])
        expected = np.column_stack([(state == 1).mean(axis=0), (state == 2).mean(axis=0)])
        self.assertTrue(np.allclose(pop[:, 1:], expected), "Wrong active-state populations.")

        os.makedirs("no_layout")
        os.chdir("no_layout")
        no_layout = self.write_trajectories(layout=False)
        pop2 = self.population_file(*no_layout, name="pop_active_dynamics.dat")
        self.assertTrue(np.array_equal(pop, pop2), "Different populations without states.json.")
        with open("zagpy.stderr") as f:
            self.assertIn("states.json not found", f.read())

        # The number of states of the mch_sf basis of the dynamics (3) is taken from
        # states.json, also if the highest state is never active. The basis has one state for
        # each multiplet, so there are no components to sum. In format version 2 of states.json
        # this basis was called mch.
        os.chdir(self.rundir)
        os.makedirs("mch_sf")
        os.chdir("mch_sf")
        self.write_trajectory("t01", [1, 2, 2], 1, layout=("mch_sf", 3))
        self.write_trajectory("t02", [1, 1, 3], 1, layout=("mch", 3), version=2)
        pop = self.population_file("t01", "t02", name="pop_active_mch_sf.dat")
        self.assertEqual(pop.shape, (3, 4), "Not 3 states.")
        self.assertTrue(np.allclose(pop[:, 1:], [[1, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5]]))
        self.zagpy("ensemble", "population", "--sum-mult", "--nstate", "1", "0", "1", "t01",
                   fail=True)

        # The mch basis has the spin components (a singlet and the three components of a
        # triplet), which --sum-mult sums.
        os.chdir(self.rundir)
        os.makedirs("mch")
        os.chdir("mch")
        self.write_trajectory("t01", [1, 2, 4], 1, layout=("mch", 4))
        self.write_trajectory("t02", [1, 3, 1], 1, layout=("mch", 4))
        pop = self.population_file("t01", "t02", name="pop_active_mch.dat")
        self.assertEqual(pop.shape, (3, 5), "Not 4 states.")
        pop = self.population_file("--sum-mult", "--nstate", "1", "0", "1", "t01", "t02",
                                   name="pop_active_mch.dat")
        self.assertTrue(np.allclose(pop[:, 1:], [[1, 0], [0, 1], [0.5, 0.5]]))

    def test_population_quantum(self):
        """ The 'quantum' estimator averages the rows of qpop_<basis>. A trajectory which reached
        the target state (status 2) is continued with its last row, not with its final active
        state. """
        rows = {1: [[1.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.2, 0.6, 0.2], [0.0, 0.5, 0.5]],
                2: [[1.0, 0.0, 0.0], [0.7, 0.2, 0.1], [0.4, 0.4, 0.2], [0.3, 0.3, 0.4]],
                3: [[1.0, 0.0, 0.0], [0.6, 0.3, 0.1]]}
        self.write_trajectory("t1", [1, 1, 2, 2], 1, files={"qpop_mch": rows[1]})
        self.write_trajectory("t2", [1, 1, 1, 3], 1, files={"qpop_mch": rows[2]})
        self.write_trajectory("t3", [1, 3], 2, files={"qpop_mch": rows[3]})
        pop = self.population_file("--estimator", "quantum", "--representation", "mch",
                                   "t1", "t2", "t3", name="pop_quantum_mch.dat")
        continued = rows[3] + [rows[3][-1]] * 2
        expected = (np.array(rows[1]) + np.array(rows[2]) + np.array(continued)) / 3
        self.assertTrue(np.allclose(pop[:, 0], [0.0, 0.5, 1.0, 1.5]), "Wrong times.")
        self.assertTrue(np.allclose(pop[:, 1:], expected), "Wrong quantum populations.")

        # The three states as a singlet and two doublet components: --sum-mult sums the
        # components of the mch basis, as in the other bases.
        pop = self.population_file("--estimator", "quantum", "--representation", "mch",
                                   "--sum-mult", "--nstate", "1", "1", "0", "t1", "t2", "t3",
                                   name="pop_quantum_mch.dat")
        summed = np.column_stack([expected[:, 0], expected[:, 1] + expected[:, 2]])
        self.assertTrue(np.allclose(pop[:, 1:], summed), "Wrong summed quantum populations.")

    def test_population_character(self):
        """ The 'character' estimator averages the squared rows of cstate_<basis>. A trajectory
        which reached the target state (status 2) is continued with its last row, not with 1 in
        the diabatic state with the index of its final active state. """
        rows = {1: [[1.0, 0.0], [0.6, -0.8], [0.0, 1.0]],
                2: [[1.0, 0.0], [-0.8, 0.6]]}
        self.write_trajectory("t1", [1, 1, 2], 1, files={"cstate_diab": rows[1]})
        self.write_trajectory("t2", [1, 1], 2, files={"cstate_diab": rows[2]})
        pop = self.population_file("--estimator", "character", "--representation", "diabatic",
                                   "t1", "t2", name="pop_character_diabatic.dat")
        continued = rows[2] + [rows[2][-1]]
        expected = (np.array(rows[1])**2 + np.array(continued)**2) / 2
        self.assertTrue(np.allclose(pop[:, 1:], expected), "Wrong character populations.")

        # Dynamics with one state for each multiplet (LZSH): cstate_diab is over the diabatic
        # multiplets, the diabatic_sf basis, without components to sum.
        os.makedirs("lz")
        os.chdir("lz")
        self.write_trajectory("t1", [1, 1, 2], 1, layout=("mch_sf", 2),
                              files={"cstate_diab": rows[1]})
        pop = self.population_file("--estimator", "character", "--representation", "diabatic",
                                   "t1", name="pop_character_diabatic_sf.dat")
        self.assertTrue(np.allclose(pop[:, 1:], np.array(rows[1])**2))
        self.zagpy("ensemble", "population", "--estimator", "character", "--representation",
                   "diabatic", "--sum-mult", "--nstate", "2", "0", "0", "t1", fail=True)

    def test_population_errors(self):
        """ Combinations of options and inputs which are errors. """
        self.write_trajectory("t1", [1, 2], 1, layout=("adiabatic", 2),
                              files={"qpop_mch": [[1.0, 0.0], [0.5, 0.5]],
                                     "qpop_diab": [[1.0, 0.0]]})
        self.write_trajectory("t2", [1, 2], 1, layout=("mch_sf", 2))
        cases = [
            (["--representation", "adiabatic", "t1"], "needs --estimator character or quantum"),
            (["--estimator", "character", "t1"], "needs --representation"),
            (["--estimator", "character", "--representation", "mch", "t1"],
             "not available in the mch basis"),
            (["--estimator", "quantum", "--representation", "mch", "--sum-mult", "--nstate",
              "1", "0", "1", "t1"], "gives 4 states"),
            (["--sum-mult", "--nstate", "2", "0", "0", "t2"],
             "--sum-mult is not available in the mch_sf basis"),
            (["--sum-mult", "--nstate", "1", "0", "1", "t1"], "gives 4 states"),
            (["t1", "t2"], "different bases of the dynamics"),
            (["--estimator", "quantum", "--representation", "diabatic", "t1"], "has 1 rows"),
        ]
        for i, (args, message) in enumerate(cases):
            with self.subTest(args=args):
                args = [os.path.abspath(a) if a in ("t1", "t2") else a for a in args]
                self.zagpy("ensemble", "population", *args, cwd=f"case{i}", fail=True)
                with open(os.path.join(f"case{i}", "zagpy.stderr")) as f:
                    self.assertIn(message, f.read())

    def test_seed_argument(self):
        """ A negative seed is an error. """
        self.write_points(n_point=1)
        self.zagpy("sample", "state", "--seed", "-1", "p01", fail=True)
        with open("zagpy.stderr") as f:
            self.assertIn("invalid seed", f.read())


if __name__ == '__main__':
    unittest.main()
