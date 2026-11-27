#!/usr/bin/env python3
''' Tests of the --seed option of the zagpy commands that draw random numbers. '''
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
    """ Reproducibility of zagpy commands with --seed. The input is generated in the tests. """

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

    def write_trajectories(self, n_traj=10, n_step=6):
        """ Write Results/energy.dat (time, active state) and Results/status (1) of two-state
        trajectories t01, t02, ... with random active states. """
        state = np.random.default_rng(0).integers(1, 3, size=(n_traj, n_step))
        state[0] = 1
        state[1] = 2
        dirs = [f"t{i + 1:02d}" for i in range(n_traj)]
        for d, s in zip(dirs, state):
            os.makedirs(os.path.join(d, "Results"))
            np.savetxt(os.path.join(d, "Results", "energy.dat"),
                       np.column_stack([np.arange(n_step) * 0.5, s]), fmt="%.2f %d",
                       header="time state")
            with open(os.path.join(d, "Results", "status"), "w") as f:
                f.write("1\n")
        return dirs

    def population(self, out, dirs, *args):
        """ Run ensemble population with --bootstrap into directory out. Returns the seed and the
        lower and upper bounds. """
        os.makedirs(out)
        seed = self.zagpy("ensemble", "population", "--bootstrap", "--output",
                          os.path.join(out, "pop"), *args, *dirs)
        bounds = [np.loadtxt(os.path.join(out, f"pop_adiabatic_{b}.dat"))
                  for b in ("lbound", "ubound")]
        return seed, bounds

    def test_population_seed(self):
        """ The same seed gives the same bootstrap intervals, a different seed different ones,
        and the printed seed of a run without --seed repeats it. Without --bootstrap, --seed is
        accepted and no seed is printed. """
        dirs = self.write_trajectories()
        files = ["pop_adiabatic_lbound.dat", "pop_adiabatic_ubound.dat"]
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

    def test_seed_argument(self):
        """ A negative seed is an error. """
        self.write_points(n_point=1)
        self.zagpy("sample", "state", "--seed", "-1", "p01", fail=True)
        with open("zagpy.stderr") as f:
            self.assertIn("invalid seed", f.read())


if __name__ == '__main__':
    unittest.main()
