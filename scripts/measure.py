#!/usr/bin/python3
import argparse
import sys
import numpy as np
import align
import io_xyz
from normalmode import NormalModes


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=
        "Read geometries from xyz files and print requested properties.",
        epilog=
        "One line is printed per input geometry as: file, geometry #, bonds, angles, dihedrals.")
    parser.add_argument(
        "-ref",
        "--reference-structure",
        type=str,
        metavar=("ref_file"),
        help="File containing reference structure for displacement measures.")
    parser.add_argument(
        "-l",
        "--length",
        nargs=2,
        type=int,
        action="append",
        metavar=("a1", "a2"),
        help="Measure length between two atoms.")
    parser.add_argument(
        "-a",
        "--angle",
        nargs=3,
        type=int,
        action="append",
        metavar=("a1", "a2", "a3"),
        help="Measure angle between three atoms.")
    parser.add_argument(
        "-d",
        "--dihedral",
        nargs=4,
        type=int,
        action="append",
        metavar=("a1", "a2", "a3", "a4"),
        help="Measure dihedral angle between four atoms.")
    parser.add_argument(
        "-nm",
        "--normal-mode",
        type=str,
        action="append",
        metavar=("nm-index"),
        help="Measure displacement along given normal mode from reference structure."
             + " To print all displacements use '-nm all'.")
    parser.add_argument(
        "-rmsd",
        "--rmsd",
        action="store_true",
        help="Measure root-mean-square deviation from reference structure.")
    parser.add_argument(
        "-align",
        "--align",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Align to reference structure before calculation.")
    parser.add_argument(
        "-last",
        "--last-only",
        action="store_true",
        help="Only print output for last geometry of each file.")
    parser.add_argument(
        "file",
        nargs="+",
        type=str,
        help="Input files in XYZ format.")
    args = parser.parse_args()

    # Reference structure needs to be in a molden file with normal mode
    # information if any -nm are present
    if args.normal_mode is not None:
        if args.reference_structure is None:
            print("error: the -ref argument is required when -nm is present")
            sys.exit()
        ref_nmode = NormalModes.from_molden(args.reference_structure)
        rgeom = ref_nmode.rgeom
    elif args.reference_structure:
        rgeom = io_xyz.ReadAll(args.reference_structure)[3][0]

    # Convert from one-based indexing in input to zero-based.
    args.length = np.array(args.length or [], dtype=int) - 1
    args.angle = np.array(args.angle or [], dtype=int) - 1
    args.dihedral = np.array(args.dihedral or [], dtype=int) - 1
    if args.normal_mode is None:
        args.normal_mode = []
    elif args.normal_mode[0] == "all":
        args.normal_mode = range(7, len(ref_nmode.freq) + 7)
    args.normal_mode = np.array(args.normal_mode or [], dtype=int) - 7

    longest_name = max([len(file) for file in args.file])
    headers = ["{:{width}s}".format("# File", width=max(longest_name, 6))]
    headers.append("Ind  ")
    headers += ["B {} {}".format(*i+1) for i in args.length]
    headers += ["A {} {} {}".format(*i+1) for i in args.angle]
    headers += ["D {} {} {} {}".format(*i+1) for i in args.dihedral]
    headers += ["NM {}".format(i+7) for i in args.normal_mode]
    if args.rmsd:
        headers.append("RMSD")
    for i, head in enumerate(headers[2:]):
        headers[i+2] = "{:12s}".format(head) # @todo Fix this format
    out_string = ", ".join(headers)
    print(out_string)

    for f in args.file:
        # Read all molecules from file.
        nmol, anat, aatm, axyz, acom = io_xyz.ReadAll(f)
        for mol in range(nmol):
            if args.last_only and mol != nmol - 1:
                continue
            out_string = "{:^{width}s}".format(f, width=longest_name)
            out_string += ", {:<5d}".format(mol+1)
            if args.align and args.reference_structure:
                xyz = align.align(rgeom, axyz[mol])
            else:
                xyz = axyz[mol]
            if args.length.size > 0:
                for i in args.length:
                    val = get_dist(*xyz[i])#* 0.529177249
                    out_string += ", {:<12.7f}".format(val)
            if args.angle.size > 0:
                for i in args.angle:
                    val = get_angle(*xyz[i], degree=True)
                    out_string += ", {:<12.7f}".format(val)
            if args.dihedral.size > 0:
                for i in args.dihedral:
                    val = get_dihedral(*xyz[i], degree=True)
                    out_string += ", {:<12.7f}".format(val)
            if args.normal_mode.size > 0:
                disp_vec = ref_nmode.from_xyz(xyz)
                for i in args.normal_mode:
                    out_string += ", {:<12.7f}".format(disp_vec[i])
            if args.rmsd:
                # prealign=False because xyz is already aligned if the option
                # is selected so no need to repeat
                val = align.rmsd(rgeom, xyz, prealign=False) 
                out_string += ", {:<12.7f}".format(val)
            print(out_string)


def get_dist(point_0, point_1):
    """ Returns distance between two points. """
    return np.linalg.norm(point_0 - point_1)


def get_angle(p0, p1, p2, degree = False):
    ''' Calculate angle for 3 points. '''
    # Get vecotrs from p1 to p0 and p1 to p2.
    bond_0 = p0 - p1
    bond_1 = p2 - p1
    # Normalize the vectors.
    bond_0 /= np.linalg.norm(bond_0)
    bond_1 /= np.linalg.norm(bond_1)
    angle = np.arccos(np.dot(bond_0, bond_1))

    if degree:
        angle = np.rad2deg(angle)

    return angle


def get_dihedral(p0, p1, p2, p3, degree = False):
    ''' Calculate dihedral angle for 4 points.
    Taken from: https://stackoverflow.com/a/34245697
    '''
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    dihedral = np.arctan2(y,x)
    if degree:
        dihedral = np.rad2deg(dihedral)

    return dihedral


if __name__ == "__main__":
    main()