import numpy as np
import argparse
from normalmode import NormalModes
import os

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=
        "Read geometries and Hessians from Orca, Molden, Turbomole or Gaussian and generate initial geometries and velocities from Wigner distribution",
        epilog=
        ".")
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
        help="Name of the file containing frequencies/Hessians. .hess file for Orca, .fch/.fchk for Gaussian, .molden for molden and aoforce output containing $hessian or $hessian (projected)")
    parser.add_argument(
        "-l",
        "--log_file",
        type=str,
        metavar=("file.log"),
        default=None,
        help="Name of the .log file. Only needed for Gaussian, along with .fchk file.")
    args = parser.parse_args()
    if args.in_format == "molden":
        nm = NormalModes.from_molden(args.file_name)
        sample = sample_wigner(nm.freq, args.temperature, args.npoint)
        cwd = os.getcwd() # Current directory, from where the program was called
        for i in range(args.npoint):
            os.mkdir("GEOM_" + str(i + 1))
            os.chdir(os.path.join(cwd, "GEOM_" + str(i + 1)))
            write_coord("coord", nm.to_xyz(sample[0][i]), nm.atoms)
            write_veloc("veloc", sample[1][i], len(nm.atoms))
            os.chdir(cwd)
    if args.in_format == "turbomole":
        nm = NormalModes.from_turbomole(args.file_name)
        sample = sample_wigner(nm.freq, args.temperature, args.npoint)
        cwd = os.getcwd() # Current directory, from where the program was called
        for i in range(args.npoint):
            os.mkdir("GEOM_" + str(i + 1))
            os.chdir(os.path.join(cwd, "GEOM_" + str(i + 1)))
            write_coord("coord", nm.to_xyz(sample[0][i]), nm.atoms)
            write_veloc("veloc", sample[1][i], len(nm.atoms))
            os.chdir(cwd)
    if args.in_format == "orca":
        nm = NormalModes.from_orca(args.file_name)
        sample = sample_wigner(nm.freq, args.temperature, args.npoint)
        cwd = os.getcwd() # Current directory, from where the program was called
        for i in range(args.npoint):
            os.mkdir("GEOM_" + str(i + 1))
            os.chdir(os.path.join(cwd, "GEOM_" + str(i + 1)))
            write_coord("coord", nm.to_xyz(sample[0][i]), nm.atoms)
            write_veloc("veloc", sample[1][i], len(nm.atoms))
            os.chdir(cwd)
    if args.in_format == "gaussian":
        if args.log_file != None:
            nm = NormalModes.from_gaussian(args.file_name, args.log_file)
            sample = sample_wigner(nm.freq, args.temperature, args.npoint)
            cwd = os.getcwd() # Current directory, from where the program was called
            for i in range(args.npoint):
                os.mkdir("GEOM_" + str(i + 1))
                os.chdir(os.path.join(cwd, "GEOM_" + str(i + 1)))
                write_coord("coord", nm.to_xyz(sample[0][i]), nm.atoms)
                write_veloc("veloc", sample[1][i], len(nm.atoms))
                os.chdir(cwd)
        else:
            print("ERROR: .log file not provided!")
    return










def sample_wigner(omega, T, npoint):
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
    natom = int(geom.shape[0] / 3)
    gformat = 3*'{:20.13f} ' + '  {}\n'
    with open(fname, 'w') as ofile:
        ofile.write('$coord\n')
        for at, xyz in zip(atoms, geom.reshape([natom, 3])):
            ofile.write(gformat.format(*xyz, at))
        ofile.write('$end\n')

def write_veloc(fname, veloc, natom):
    file = open("veloc", "w")
    np.savetxt(fname, veloc.reshape(natom-2,3))




def read_coord(fname, natom):
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


if __name__ == "__main__":
    main()  
