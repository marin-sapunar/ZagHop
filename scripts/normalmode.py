import numpy as np
import scipy.linalg


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
        self.at_mass = at_mass#average atomic mass of elements
        self.mass = np.repeat(at_mass, 3) * MASS_UNIT
        self.atoms = atoms
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
    @classmethod
    def from_turbomole(cls, hess_file_name, coord_file = "coord"):
        """ Create instance of NormalMode class from turbomole files.

        Arguments:
            hessfile : Hessian file that contains $hessian (projected)".
        Returns:
            nmode : New instance of NormalMode class.
        """                                                           
        # Reads elements from Turbomole Hessian file
        try:
            a=open(hess_file_name, 'r')
        except:
            print(hess_file_name + " file not found")
        c=a.readlines()
        a.close()
        first_index = c.index("$hessian (projected)\n") # First line of the projected Hessian is the one after the line that contains $hessian (projected)
        final_index = c.index("$end\n")
        n_atom=int(int(c[final_index-1].split()[0])/3) # Number of atoms, first element of last row divided by 3
        # This block creates a Hessian matrix as a numpy array
        hessian_matrix = []
        for i in c[first_index:final_index]:
            d=i.split()
            for j in range(2,len(d)):
                hessian_matrix.append(float(d[j]))
        hessian_matrix = np.array(hessian_matrix).reshape(3*n_atom,3*n_atom)
        # Reads atoms from "coord" file in order to create a mass weighted matrix
        try:
            a=open(coord_file,"r")
        except:
            print(coord_file + " not found!")
            return exit()
        c=a.readlines()
        a.close()
        coord = []
        atom_list=[]#has each atom repeated 3 times
        atoms=[]#has each atom repeated once, used when creating instance of NormalModes class
        for i in c[1:1 + n_atom]:
            coord_line = i.split()
            upper_atom =coord_line[-1].capitalize()#transforms atom symbols to uppercase
            atoms.append(upper_atom)
            for j in range(3):
                atom_list.append(upper_atom)
                coord.append(float(coord_line[j]))
        coord = np.array(coord)

        # Creates mass weighted Hessian
        weight_matrix=[]# Used for creating mass weighted Hessian
        for i in atom_list:
            for j in atom_list:
                weight_matrix.append(1/np.sqrt(MASS[i]*MASS[j]*MASS_UNIT**2))
        weight_matrix = np.array(weight_matrix).reshape(3*n_atom,3*n_atom)
        mass_weighted_hessian = np.multiply(weight_matrix, hessian_matrix)
        force_constants, mode_vectors = scipy.linalg.eigh(mass_weighted_hessian)
        freq=np.emath.sqrt(force_constants).real
        nmode = cls(coord, mode_vectors.T, freq, atoms, mass_weighted=True)
        return nmode    
    @classmethod
    def from_orca(cls, hess_file_name):
        try:
            a = open(hess_file_name, 'r')
        except:
            print("ERROR: File not found!")
            exit()
        orca_out = a.readlines()# List that contains all the lines from .hess file
        a.close()
        try:
            hessian_index = orca_out.index("$hessian\n")# index of the line that contains $hessian, a reference point for defining ranges
        except:
            print("ERROR: $hessian not found! .hess file containing '$hessian' line needs to be provided")
            exit()
        atom_index = orca_out.index("$atoms\n")# for finding atoms, and the number of atoms
        n_atom = int(orca_out[atom_index+1].split()[0])
        # Number of atoms is contained in a line after "$atom", .split()[0] transforms  
        # it from a list to a string which is converted to an integer
        
        final_index = orca_out.index("$vibrational_frequencies\n")-1 
        # Last line of Hessian is the one that is two before the line that contains 
        # $vibrational_frequencies, but -1 is here since for i in range(X,Y) goes from X to Y-1
        
        raw_hessian = orca_out[hessian_index+2 : final_index]
        for i in range(0,len(raw_hessian)-3*n_atom,3*n_atom): 
            raw_hessian.pop(i)
            # Throws out elements that contain column labels, first it removes first line from the list, and then the next problematic
            # line is guaranteed to be the one that has 3*n_atoms greater index once the previous has been removed. Column label line
            # won't be further than len(raw_hessian)-3*n_atom with how the .hess file is organized.
        
        # this block creates a Hessian matrix as a wnumpy array
        refined_hessian=[]
        for i in range(3*n_atom):
            for j in range(i, len(raw_hessian), 3*n_atom):
                # It adds i-th line and every line that is 3*n_atom further from it, ensuring that all the rows associated with the i-th
                # coordinate are together
                temp_list = raw_hessian[j].split()
                for k in temp_list[1:]:# goes from second element of the row, so that row labels won't be added to refined_hessian list
                    refined_hessian.append(float(k))
        hessian_matrix = np.array(refined_hessian).reshape(3*n_atom, 3*n_atom)
    
        # This block will read all the atom symbols from the .hess file and create matrix that weights the Hessian
        atom_list = [] # List that contains atom symbols X 3
        atoms = [] # List that contains atom symbols once
        geometry = [] # List that contains geometry of molecule at which the Hessian was calculated
        for i in orca_out[atom_index + 2: atom_index + 2 + n_atom]:
            atoms.append(i.split()[0])
            for j in range(3):# It will add each atom 3 times to the atom_list to make creating weight matrix easier
                atom_list.append(i.split()[0])
            for k in i.split()[2:]:
                geometry.append(float(k))
        geometry = np.array(geometry).reshape(n_atom, 3)

        # Creates mass weighted Hessian
        weight_matrix=[]# Used for creating mass weighted Hessian
        for i in atom_list:
            for j in atom_list:
                weight_matrix.append(1/np.sqrt(MASS[i]*MASS[j]*MASS_UNIT**2))
        weight_matrix = np.array(weight_matrix).reshape(3*n_atom,3*n_atom)
        mass_weighted_hessian = np.multiply(weight_matrix, hessian_matrix)
        force_constants, mode_vectors = scipy.linalg.eigh(mass_weighted_hessian)
        freq=np.emath.sqrt(force_constants).real
        nmode = cls(geometry, mode_vectors.T, freq, atoms, mass_weighted=False)
        return nmode
    @classmethod
    def from_gaussian(cls, gauss_file_name, log_file_name):
        # hess_ file_name is fch file name, generated with Right Click > Results > View/Edit file in .chk opened in GaussView on Windows.
        # log_file_name is the standard .log file. Rename them to default values ("gaussian.fch.txt" and "STRUCTURE.LOG") or call the
        # function with from_gaussian("custom_name1", "custom_name2")
        try:
            a = open(gauss_file_name, "r")
        except:
            print("ERROR: Formatted checkpoint file not found!")
        gauss_out = a.readlines()
        a.close()

        hessian_index = gauss_out.index([i for i in gauss_out if "Cartesian Force Constants" in i][0])
        # Finds the line that contains substring "Cartesian Force Constants", which is the line before Hessian matrix is printed

        n_atom = int(gauss_out[gauss_out.index([i for i in gauss_out if "Number of atoms" in i][0])].split()[-1])
        # Number of atoms is the integer of the last word in the line that contains n_atom
        a = open(log_file_name)
        gauss_log = a.readlines()
        a.close()

        Z_mat_index = gauss_log.index(" Symbolic Z-matrix:\n") # Index of Z-matrix line in .log file - used to find all the atoms
        atom_list = []
        atoms = []
        for i in gauss_log[Z_mat_index + 2: Z_mat_index + 2 + n_atom]:
            atoms.append(i.split()[0])
            for j in range(3):
                atom_list.append(i.split()[0])
        # Adds all the atom symbols from Symbolic Z-matrix found in .LOG file

        # This block creates Hessian matrix
        lt_hessian_raw = []# unformatted lower triangular Hessian
        n_rows = int((3 * n_atom + (3 * n_atom * 3 * n_atom - 3 * n_atom)/2)//5)
        if (3 * n_atom + (3 * n_atom * 3 * n_atom - 3 * n_atom)/2)%5 != 0:
            n_rows += 1
        # Gaussian stores only lower triangular half of the Hessian matrix, and always in 5 columns, so the number of elements is 
        # (3n_atoms)^2 - 3n_atoms divided by 2 (number of off-diagonal elements/2) + 3n_atoms (number of diagonal elements) all divided by 5
        # and one row is added if the resulting number isn't divisible by 5

        for i in gauss_out[hessian_index + 1: hessian_index + 1 + n_rows]:
            c = i.split()
            for j in c:
                lt_hessian_raw += [float(j)]

        lt_hessian = []# Formatted lower triangular Hessian

        for i in range(3 * n_atom):
            hess_row = []
            for j in range(3 * n_atom):
                if i >= j:
                    hess_row.append(lt_hessian_raw.pop(0))
                    # To create 3n_atom X 3n_atom numpy array that contains Hessian matrix elements in the lower triangle,
                    # the program adds elements from the beggining of the lt_hessian_raw to hess_row list and removes them,
                    # while i is greater than or equal to j meaning (while column is greater than or equal to row, so that it adds those)
                    # elements up until the diagonal. After the diagonal element, the row is filled with 0. Full Hessian is created by
                    # adding the lower triangular matrix and upper triangular matrix of Hessian (transpose of lt_hessian) and subtracting
                    # the diagonal elements of one of those matrices
                else:
                    hess_row.append(0)
            lt_hessian.append(hess_row)
        lt_hessian = np.array(lt_hessian)# Transforms regular lt_hessian array to a numpy array
        hess_diag = np.diag(np.diagonal(lt_hessian))# Diagonal matrix that contains only the diagonal elements of lt_hessian
        hessian_matrix = lt_hessian + lt_hessian.T - hess_diag
        # Full hessian is made by summing lower triangular Hessian and it's transpose, then subtracting the diagonal elements
        # of lt_hessian since they are added twice

        weight_matrix=[]
        for i in atom_list:
            for j in atom_list:
                weight_matrix.append(1/np.sqrt(MASS[i]*MASS[j]*MASS_UNIT**2))
        weight_matrix = np.array(weight_matrix).reshape(3*n_atom,3*n_atom)
        mass_weighted_hessian = np.multiply(weight_matrix, hessian_matrix)
        force_constants, mode_vectors = scipy.linalg.eigh(mass_weighted_hessian)
        #reads coordinates
        coord_index = gauss_out.index([i for i in gauss_out if "Current cartesian coordinates" in i][0])
        geometry = []
        n_rows = int((3 * n_atom // 5))
        if 3*n_atom%5 != 0:
            n_rows += 1
        for i in range(coord_index + 1, coord_index + 1 + n_rows):
            for j in gauss_out[i].split():
                geometry.append(float(j))
        geometry = np.array(geometry).reshape(n_atom, 3)
        weight_matrix = np.array(weight_matrix).reshape(3*n_atom,3*n_atom)
        mass_weighted_hessian = np.multiply(weight_matrix, hessian_matrix)
        force_constants, mode_vectors = scipy.linalg.eigh(mass_weighted_hessian)
        freq=np.emath.sqrt(force_constants).real
        nmode = cls(geometry, mode_vectors.T, freq, atoms, mass_weighted=False)
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
