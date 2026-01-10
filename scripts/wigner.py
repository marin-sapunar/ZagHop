#!/usr/bin/env python3
""" Script for sampling an ensemble of phase space points
    from a Harmonic Wigner distribution.

    Usage: call wigner.py --help for details. """
import argparse
import os
import sys
import numpy as np
from normalmode import NormalModes

def main():
    """ Adds various options for input format and sampling parameters and then calls
    various functions defined below to create initial conditions from Harmonic
    Wigner distribution. """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=
        """Read geometries and Hessians from Orca, Molden, Turbomole or Gaussian and
        generate initial geometries and velocities from Wigner distribution. For Molden
        and Orca, only the file containing Hessian/normal modes needs to be provided.
        For Turbomole, a 'coord' file containing ground state geometry is needed. For
        Gaussian, one needs to specify .log file's name using -l option.""",
        epilog=
        "")
    parser.add_argument(
        "-f",
        "--in_format",
        choices=["molden", "orca", "gaussian", "turbomole"],
        type=str,
        metavar=("file_format"),
        default="molden",
        help="Format of the file. Supported formats: gaussian, molden, turbomole, orca.")
    parser.add_argument(
        "-n",
        "--npoint",
        type=int,
        metavar=("N"),
        default=1000,
        help="Number of points to be generated from Wigner distribution.")
    parser.add_argument(
        "-T",
        "--temperature",
        type=int,
        metavar=("Temperature"),
        default=0,
        help="Temperature in Kelvin.")
    parser.add_argument(
        "-fn",
        "--file_name",
        type=str,
        metavar=("file_name"),
        default="vib.molden",
        help="""Name of the file containing frequencies/Hessians. .hess file for Orca,
                .fch/.fchk for Gaussian, .molden for molden and aoforce output containing $hessian 
                or $hessian (projected)""")
    parser.add_argument(
        "-l",
        "--log_file",
        type=str,
        metavar=("file.log"),
        default=None,
        help="Name of the .log file. Only needed for Gaussian, along with .fchk file.")
    parser.add_argument(
        "-c",
        "--coord_file",
        type=str,
        metavar=("coord"),
        default=None,
        help="""Location of the coord file containing optimized ground
                state geometry. Only needed for Turbomole.""")
    parser.add_argument(
        "-ign",
        "--ignore_modes",
        type=str,
        metavar=("mode_list"),
        default=None,
        help="""Indices of normal modes to be ignored. Example: -ign '1 3'
                will ignore first and third normal mode.""")

    args = parser.parse_args()
    if args.in_format == "molden":
        nm = NormalModes.from_molden(args.file_name)
    elif args.in_format == "turbomole":
        nm = NormalModes.from_turbomole(args.file_name, args.coord_file)
    elif args.in_format == "orca":
        nm = NormalModes.from_orca(args.file_name)
    elif args.in_format == "gaussian":
        if args.log_file is not None:
            nm = NormalModes.from_gaussian(args.file_name, args.log_file)
        else:
            print("ERROR: .log file not provided!")
            sys.exit(1)
    sample = sample_wigner(nm.freq, args.temperature, args.npoint)
    cwd = os.getcwd() # Current directory, from where the program was called
    for i in range(args.npoint):
        os.mkdir("point" + str(i + 1).zfill(len(str(args.npoint))))
        os.chdir(os.path.join(cwd, "point" + str(i + 1).zfill(len(str(args.npoint)))))
        write_veloc("veloc", nm.to_xyz(refined_modes(sample[1][i], args.ignore_modes),
                    reshape=True, displacement=True))
        write_geom("geom", nm.to_xyz(refined_modes(sample[0][i],
            args.ignore_modes)), nm.atoms, nm.at_mass)
        os.mkdir("qmdir")
        os.chdir("qmdir")
        write_coord("coord", nm.to_xyz(refined_modes(sample[0][i], args.ignore_modes)), nm.atoms)
        os.chdir(cwd)


def sample_wigner(omega, T, npoint):
    """ Takes normal mode frequencies, temperature and number of points and
        returns arrays containing normal mode displacements and velocities
        sampled from Harmonic Wigner distribution """
    sample_q = np.zeros((len(omega), npoint))
    sample_v = np.zeros((len(omega), npoint))
    T_au = T * 3.1668116847144e-6 # @todo units module
    for i, o in enumerate(omega):
        sigma_q = thermal_wigner_q(o, T_au)
        sigma_v = thermal_wigner_v(o, T_au)
        rng = np.random.default_rng()
        sample_q[i] = rng.normal(scale=sigma_q, size=npoint) * np.sqrt(o)
        sample_v[i] = rng.normal(scale=sigma_v, size=npoint) * np.sqrt(o)
    return sample_q.T, sample_v.T


def thermal_wigner_q(omega, T):
    if T > 0:
        return 1 / np.sqrt(omega * 2 * np.tanh(omega / 2 / T))
    return 1 / np.sqrt(omega * 2)


def thermal_wigner_v(omega, T):
    if T > 0:
        return np.sqrt(omega / (2 * np.tanh(omega / (2 * T))))
    return np.sqrt(omega / 2)



def write_coord(fname, geom, atoms):
    """ Writes a coord file for Turbomole input. Takes output file name,
        numpy array containing geometry and array containing atoms and
        creates a coord file. """
    natom = int(geom.shape[0] / 3)
    gformat = 3*'{:20.13f} ' + '  {}\n'
    with open(fname, 'w') as ofile:
        ofile.write('$coord\n')
        for at, xyz in zip(atoms, geom.reshape([natom, 3])):
            ofile.write(gformat.format(*xyz, at.lower()))
        ofile.write('$end\n')

def write_veloc(fname, veloc):
    """ Takes output file name and a numpy array containing
        normal mode velocities and saves them in a file. """
    np.savetxt(fname, veloc)

def write_geom(fname, geom, atoms, at_mass):
    """ Takes as input fname (output file name), geom (a numpy array
        containing coordinates of atoms), atoms (an array containing atom
        labels) and at_mass (an array containing relative atomic masses)
        and creates a file named fname suitable for molecular dynamics
        simulations with zaghop. """
    natom = int(geom.shape[0] / 3)
    gformat = '{} ' + '{} ' + 3*'{:.13f} ' + 'q\n'
    with open(fname, 'w') as ofile:
        for at, m_at, xyz in zip(atoms, at_mass, geom.reshape([natom, 3])):
            ofile.write(gformat.format(at.lower(), m_at, *xyz))

def read_coord(fname, natom):
    """ fname is the name of the file to be read, and natom is the number of
        atoms in the system, Reads a coord file named fname and returns arrays
        containing atom symbols and geometry of the system. """
    rgeom = np.zeros(3*natom)
    atoms = []
    with open(fname, 'r') as ifile:
        for line in ifile:
            if line.strip().startswith('$coord'):
                for i in range(natom):
                    line = ifile.readline().split()
                    atoms.append(line[3])
                    xyz = line[0:3]
                    xyz = [np.float(q) for q in xyz]
                    rgeom[3*i:3*i+3] = xyz
    return atoms, rgeom

def refined_modes(omega, ignore_list):
    """ Takes a numpy array of normal mode frequencies and list of mode numbers
    to be ignored (with number ranging from 1 to 3*n_atom - 6) and returns a list
    with specified frequencies set to zero. """
    if ignore_list is not None:
        for i in ignore_list.split():
            omega[int(i)-1] = 0
    return omega
if __name__ == "__main__":
    main()
