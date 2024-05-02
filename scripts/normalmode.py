import numpy as np


class NormalModes():
    """ Class for handling normal modes.

    Attributes:
        rgeom : Reference geometry at which the normal modes are defined.
        nmode : Array of normal mode vectors in mass-weighted coordinates.
        freq : Vibrational frequencies of the normal modes.
        mass : Mass associated with each coordinate.
        d_to_c : Conversion from dimensionless n. mode to cart displacements.
        c_to_d : Conversion from cart to dimensionless n. mode displacements.
    """

    def __init__(self, rgeom, nmode, freq, atoms, mass_weighted=False):
        """ Initialize data in NormalModes instance.

        Args:
            rgeom : Reference geometry at which the normal modes are defined.
            nmode : Normal mode vectors. Input vectors are assumed to be in
                    Cartesian coordinates unless mass_weighted is True.
            freq : Normal mode frequencies.
            atoms : List of atoms (needed for masses).
            mass_weighted : If normal modes passed in nmode are already
                            weighted skip the mass weighting step.
        """
        if rgeom.ndim > 2:
            raise ValueError('Wrong number of dimensions for rgeom.')
        self.rgeom = rgeom.flatten()
        if self.rgeom.shape[0] % 3 != 0:
            raise ValueError('Wrong dimensions for rgeom.')
        natom = int(self.rgeom.shape[0] / 3)
        if len(atoms) != natom:
            raise ValueError('Dimension mismatch between rgeom and atoms')
        if nmode.shape[1] != 3 * natom:
            raise ValueError('Dimension mismatch between rgeom and nmode')
        nfreq = nmode.shape[0]
        if len(freq) != nfreq:
            raise ValueError('Dimension mismatch between freq and nmode.')
        at_mass = [MASS[atom.title()] for atom in atoms]
        self.mass = np.repeat(at_mass, 3) * MASS_UNIT
        if mass_weighted:
            self.nmode = nmode.copy()
        else:
            self.nmode = nmode * np.sqrt(self.mass)
            norm = np.linalg.norm(self.nmode, axis=1)
            self.nmode = (self.nmode.T / norm).T
        self.freq = freq.copy()
        # Remove rotation/translation
        self.freq = self.freq[6:]
        self.nmode = self.nmode[6:]
        conversion = np.einsum('i,j->ij', self.freq, self.mass)
        conversion = np.sqrt(np.abs(conversion))  # @todo abs here as option?
        self.d_to_c = (self.nmode / conversion).T
        self.c_to_d = self.nmode * conversion

    @classmethod
    def from_molden(cls, infile):
        """ Create instance of NormalMode class from a molden file.

        Molden file needs to contain the fr-norm-coord section.

        Arguments:
            infile : Input molden file.
        Returns:
            nmode : New instance of NormalMode class.
        """
        lines = _get_section('fr-coord', infile, quiet=True)
        if lines:
            natom, atoms, coord = _parse_fr_coord(lines)
        else:
            lines = _get_section('atoms', infile, title=True, quiet=False)
            natom, atoms, coord = _parse_atoms(lines)
        # Normal modes.
        lines = _get_section('fr-norm-coord', infile, quiet=True)
        if lines:
            vib = _parse_fr_norm_coord(lines, natom)
            lines = _get_section('freq', infile)
            freq = _parse_freq(lines)
            nmode = cls(coord, vib, freq, atoms, mass_weighted=False)
        return nmode


    def to_xyz(self, q_vec, reshape=False, displacement=False):
        """ Convert NM displacement vector to Cartesian coordinates.

        Argumnets:
            q_vec : Displacements in dimensionless normal mode coordinates.
            reshape : Reshape output coordinates from a 3*natom vector to a
                shape (natom, 3) array.
            displacement : Output displacements (from reference geometry)
                in Cartesian coordinates instead of absolute geometry.

        Returns:
            cart : Geometry in Cartesian coordinates.
        """
        cart = np.dot(self.d_to_c, q_vec)
        if not displacement:
            cart = cart + self.rgeom
        if reshape:
            cart = cart.reshape((-1, 3))
        return cart


    def from_xyz(self, geom, displacement=False):
        """ Convert Cartesian coordinates to NM displacement vector.

        Arguments:
            cart : Geometry in Cartesian coordinates.
            displacement : Input is in displacements (from reference geometry)
                in Cartesian coordinates instead of absolute geometry.
        Returns:
            q_vec : Displacements in dimensionless normal mode coordinates.
        """
        q_vec = geom.flatten()
        if not displacement:
            q_vec = q_vec - self.rgeom
        q_vec = np.dot(self.c_to_d, q_vec)
        return q_vec


    def harmonic_potential_energy(self, q_vec=None, cart=None):
        """ Calculate PE in harmonic approximation for given geometry.

        

        Arguments:
            q_vec : Displacements in dimensionless normal mode coordinates.
            cart : Geometry in Cartesian coordinates.
        Returns:
            pot_en : Harmonic potential energy along each normal mode.

        """
        if cart is None and q_vec is None:
            raise ValueError("Neither 'q_vec' nor 'cart' given on input")
        if cart is not None:
            q_vec = self.from_xyz(cart)
        pot_en = q_vec**2 * self.freq * 0.5
        return pot_en


MASS_UNIT = 1822.885
MASS = {'H' : 1.008,
        'He': 4.002,
        'Li': 6.941,
        'Be': 9.012,
        'B' : 10.811,
        'C' : 12.011,
        'N' : 14.007,
        'O' : 15.999,
        'F' : 18.998,
        'Ne': 20.18,
        'Na': 22.99,
        'Mg': 24.305,
        'Al': 26.982,
        'Si': 28.086,
        'P' : 30.974,
        'S' : 32.065,
        'Cl': 35.453,
        'Ar': 39.948,
        'K' : 39.098,
        'Ca': 40.078,
        'Sc': 44.956,
        'Ti': 47.867,
        'V' : 50.942,
        'Cr': 51.996,
        'Mn': 54.938,
        'Fe': 55.845,
        'Co': 58.933,
        'Ni': 58.693,
        'Cu': 63.546,
        'Zn': 65.38,
        'Ga': 69.723,
        'Ge': 72.64,
        'As': 74.922,
        'Se': 78.96,
        'Br': 79.904,
        'Kr': 83.798,
        'Rb': 85.468,
        'Sr': 87.62,
        'Y' : 88.906,
        'Zr': 91.224,
        'Nb': 92.906,
        'Mo': 95.96,
        'Tc': 98.,
        'Ru': 101.07,
        'Rh': 102.906,
        'Pd': 106.42,
        'Ag': 107.868,
        'Cd': 112.411,
        'In': 114.818,
        'Sn': 118.71,
        'Sb': 121.76,
        'Te': 127.6,
        'I' : 126.904,
        'Xe': 131.293,
        'Cs': 132.905,
        'Ba': 137.327,
        'La': 138.905,
        'Ce': 140.116,
        'Pr': 140.908,
        'Nd': 144.242,
        'Pm': 145.,
        'Sm': 150.36,
        'Eu': 151.964,
        'Gd': 157.25,
        'Tb': 158.925,
        'Dy': 162.5,
        'Ho': 164.93,
        'Er': 167.259,
        'Tm': 168.934,
        'Yb': 173.054,
        'Lu': 174.967,
        'Hf': 178.49,
        'Ta': 180.948,
        'W' : 183.84,
        'Re': 186.207,
        'Os': 190.23,
        'Ir': 192.217,
        'Pt': 195.084,
        'Au': 196.967,
        'Hg': 200.59,
        'Tl': 204.383,
        'Pb': 207.2,
        'Bi': 208.98,
        'Po': 210.,
        'At': 210.,
        'Rn': 222.,
        'Fr': 223.,
        'Ra': 226.,
        'Ac': 227.,
        'Th': 232.038,
        'Pa': 231.036,
        'U' : 238.029,
        'Np': 237.,
        'Pu': 244.,
        'Am': 243.,
        'Cm': 247.,
        'Bk': 247.,
        'Cf': 251.,
        'Es': 252.,
        'Fm': 257.,
        'Md': 258.,
        'No': 259.,
        'Lr': 262.,
        'Rf': 261.,
        'Db': 262.,
        'Sg': 266.,
        'Bh': 264.,
        'Hs': 267.,
        'Mt': 268.,
        'Ds': 271.,
        'Rg': 272.,
        'Cn': 285.,
        'Nh': 284.,
        'Fl': 289.,
        'Mc': 288.,
        'Lv': 292.,
        'Ts': 295.,
        'Og': 294.}


def _parse_freq(lines):
    """ Parse [FREQ] section. Just converts to float. """
    unit = 4.556335e-6 # @todo Units module
    freq = np.array(lines, dtype=np.float64) * unit
    return freq


def _parse_atoms(lines):
    """ Parse a molden [ATOMS] section.

    Keyword line needs to be included in the lines argument to determine units.
    """
    UNITS = {"ANGS" : 1.88972612456,
             "AU" : 1.0}
    unit_str = lines[0].split()[1].upper()
    unit = UNITS[unit_str]
    natom = len(lines) - 1
    coord = np.zeros((natom, 3))
    atoms = []
    for i, line in enumerate(lines[1:]):
        line = line.split()
        atoms.append(line[0].title())
        coord[i] = np.array(line[3:], dtype=np.float64)
    coord = coord * unit
    return natom, atoms, coord


def _parse_fr_coord(lines):
    """ Parse a molden [FR_COORD] section. """
    natom = len(lines)
    coord = np.zeros((natom, 3))
    atoms = []
    for i, line in enumerate(lines[0:]):
        line = line.split()
        atoms.append(line[0].title())
        coord[i] = np.array(line[1:], dtype=np.float64)
    return natom, atoms, coord


def _parse_fr_norm_coord(lines, natom):
    """ Parse a molden [FR-NORM-COORD] section. """
    nfreq = int(len(lines) / (natom + 1))
    vib = np.zeros((nfreq, natom, 3))
    for i in range(nfreq):
        for j in range(natom):
            xyz = lines[i * (natom + 1) + j + 1].split()
            vib[i, j] = np.array(xyz, dtype=float)
    vib = vib.reshape((nfreq, natom * 3))
    return vib


def _get_section(section, path, title=False, quiet=False):
    """ Get a section from a molden file.

    Args:
        section : Section to find.
        path : Path to the molden file.
        quiet : Don't throw error if section is not found.

    Returns:
        lines : List of lines in the found section.
    """
    lines = []
    with open(path, 'r') as cfile:
        try:
            line = _find_section(cfile, section)
        except ValueError:
            if quiet:
                return lines
            raise
        if title:
            lines.append(line)
        for line in cfile:
            line = line.strip()
            if line.startswith(r'#'):
                pass
            elif line.startswith(r'['):
                break
            lines.append(line)
    return lines


def _find_section(cfile, section):
    """ Go to a specific [section] in an opened molden file.

    Args:
        cfile : File to seach.
        section : Name of section to find.

    Returns:
        line : Line containing the section keyword.
    """
    sec = r'[' + section.lower() + r']'
    for line in cfile:
        line = line.strip().lower()
        if line.startswith(sec):
            return line
    raise ValueError(f'Requested section not found: {sec}.')
