""" Random number generators for zagpy commands.

Every command that draws random numbers adds the --seed option with add_seed_argument and
creates its generator with make_rng. Without --seed, a seed is drawn and printed (as in ZagHop),
so that passing it back as --seed repeats the run.
"""
import argparse
import numpy as np


def _seed(value):
    """ Parse a seed from the command line: an integer >= 0. """
    try:
        seed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid seed: '{value}' (integer >= 0)") from None
    if seed < 0:
        raise argparse.ArgumentTypeError(f"invalid seed: {seed} (integer >= 0)")
    return seed


def add_seed_argument(parser):
    """ Add the --seed option to the parser of a command. """
    parser.add_argument(
        "--seed",
        type=_seed,
        default=None,
        metavar="SEED",
        help="""Seed of the random number generator (integer >= 0). Without it, a seed is
                drawn and printed. The same seed gives the same numbers only with the same
                NumPy and SciPy versions.""")


def make_rng(seed=None):
    """ Return a numpy Generator created from seed and print the seed.

    If seed is None, a new seed is drawn from the entropy of the operating system. The printed
    seed passed back as seed gives the same generator."""
    if seed is None:
        seed = np.random.SeedSequence().entropy
    print(f"Random seed: {seed}")
    return np.random.default_rng(seed)
