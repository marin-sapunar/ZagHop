""" IO for XYZ files. """
import numpy as np


class toau:
    angs = 1.88972612456


def ReadAll(path):
    """ Read all information from an XYZ file. """
    anat = []
    aatm = []
    axyz = []
    acom = []
    with open(path, 'r') as ifile:
        for molecule in _read_molecule(ifile):
            anat.append(molecule[0])
            aatm.append(molecule[1])
            axyz.append(molecule[2])
            acom.append(molecule[3])
    nmol = len(anat)
    return nmol, anat, aatm, axyz, acom


def _read_molecule(ifile):
    """ Read a single molecule from a XYZ file.

    Args:
        ifile : File object opened for reading.

    Yields:
        nat : Number of atoms.
        atm : List of atom symbols.
        xyz : Array of atom coordinates.
        com : Comment line from the file.
    """
    while True:
        nat = ifile.readline()
        if not nat:
            break
        nat = int(nat)
        com = ifile.readline().strip()
        atm = []
        xyz = np.zeros([nat, 3], dtype='float64')
        for i in range(nat):
            sline = ifile.readline().split()
            atm.append(sline[0].title())
            xyz[i] = sline[1:4]
        xyz = xyz * toau.angs
        yield nat, atm, xyz, com


def _write_molecule(ofile, nat, atm, xyz, com):
    """ Write a single molecule to a file in XYZ format.

    Args:
        ofile : File object opened for writing.
        nat : Number of atoms.
        atm : List of atom symbols.
        xyz : Array of atom coordinates.
        com : Comment line from the file.
    """
    ofile.write('{}\n'.format(nat))
    ofile.write('{}\n'.format(com))
    for i in range(nat):
        axyz = xyz[i] / toau.angs
        ofile.write('{} {} {} {}\n'.format(atm[i], axyz[0], axyz[1], axyz[2]))
