import numpy as np
from scipy.special import factorial
import os
def index_to_m(index,l):
    """
    Converts array index to angular momentum projection. Currently uses turbomole
    convention. Later it will be changed so that it is ES program agnostic, and all the
    other MOs will be converted to this format.
    """
# FIX PHASE!
    if l!=1:
        tmp_m = (index+index%2)//2
        m = tmp_m * (-1)**(tmp_m%2 + index%2)
        return m
    else:
        return [1,-1,0][index]

def create_wigner_D_matrix(l,alpha,beta,gamma):
    r"""
    Creates a Wigner D-matrix for rotating spherical harmonics, whose elements [m1,m2] are 
    <l,m1|  exp(-i \alpha L_z) exp(-i \beta L_y) exp(-i \gamma L_z) |l,m2>. \alpha, \beta
    and \gamma are the three Euler angles in the ZYZ convention, while l is the azimuthal
    quantum number.
    """
    wigner_D_matrix = np.zeros([2*l +1,2*l +1],dtype=complex)
    for ind_m1 in range(2*l+1):
        for ind_m2 in range(2*l+1):
            D_m1_m2 = 0
            m1 = index_to_m(ind_m1,l)
            m2 = index_to_m(ind_m2,l)
            s_max = min(l+m2,l-m1)
            s_min = max(0,m2-m1)
            prefactor = np.sqrt(factorial(l-m1) * factorial(l+m1) *
                                factorial(l-m2) * factorial(l+m2)
                               )
            for s in range(s_min,s_max+1):
                D_m1_m2 += ( (-1)**(m1 - m2 + s) * np.cos(beta/2)**(2*l + m2 - m1 - 2*s) * np.sin(beta/2) ** (m1-m2+2*s)/
                        (factorial(l + m2 - s) * factorial(s) * factorial(m1 - m2 + s) * factorial(l - m1 - s))  
                        )
            wigner_D_matrix[ind_m1,ind_m2] = prefactor * D_m1_m2 * np.e**(-1J * alpha * m1 - 1J * gamma * m2)
    return wigner_D_matrix

def read_coord(coord_file):
    with open(coord_file) as cfile:
        coordlines = cfile.readlines()
    atom_list = []
    coord_matrix = []
    for line in coordlines:
        if "$" not in line:
            coord_line = line.split()
            atom_list.append(coord_line[-1].upper())
            coord_matrix.append(coord_line[:-1])
    coord_matrix = np.array(coord_matrix,dtype = float)
    #print(coord_matrix)
    #print(atom_list)
    return atom_list,coord_matrix
def symbol_to_ang_mom(symbol):
    if symbol == "s":
        return 0
    elif symbol == "p":
        return 1
    elif symbol == "d":
        return 2
    elif symbol == "f":
        return 3
    elif symbol == "g":
        return 4
    elif symbol == "h":
        return 5

def read_basis(basis_file):
    with open(basis_file) as bfile:
        basislines = bfile.readlines()
    AO_basis = {}
    basis_index = -1
    for line_index in range(len(basislines)):
        line = basislines[line_index]
        if "#" in line:
            basis_funcs = []
            basis_center = line.split()[1].upper()
            basis_index = line_index + 2
        if line_index == basis_index:
            basisline = basislines[basis_index].split()
            if  "*" not in basisline:
                basis_funcs.append(symbol_to_ang_mom(basisline[1]))
                basis_index += int(basisline[0]) + 1
            else:
                AO_basis[basis_center]=(basis_funcs)
    return AO_basis

def get_azimuthal_q_num_list(AO_basis, atom_list):
    print(AO_basis)
    azimuthal_quantum_number_list = []
    for atom in atom_list:
        azimuthal_quantum_number_list += sorted(AO_basis[atom])
    return azimuthal_quantum_number_list
def split_before_dot(lines):
    out_array = []
    for line in lines:
        start = 0
        line.strip("\n")
        line=line.replace("D","e")
        for i in range(len(line)):
            if line[i] == "." and i>1:
                out_array.append(line[start:i-1])
                start = i-1
        out_array.append(line[start:])
    return out_array

def read_MOs(mos_file):
    # COLUMNS CONTAIN BASIS SET COEFFS
    with open(mos_file) as ifile:
        mos_lines = ifile.readlines()
    coeffs = []
    for line_ind in range(len(mos_lines)):
        if "nsaos" in mos_lines[line_ind]:
            ao_indices = []
            n_AOs = int(mos_lines[line_ind].split("=")[-1])
            n_lines = n_AOs//4 + (n_AOs%4!=0)
            tmp_coeffs = mos_lines[line_ind+1:line_ind+1+n_lines]

            coeffs.append(np.array(split_before_dot(tmp_coeffs),dtype =float).flatten())
    coeffs = np.array(coeffs).T
    return coeffs

def calculate_euler_angles(coords_initial, coords_final):
    # LATER - GENERALIZE TO CHECK FOR LINEAR INDEPENDENCE
    v_initial = []
    v_final = []
    for i in range(1,3):
        #for j in range(i-1):
        v_initial.append(coords_initial[i]- coords_initial[i-1])
        v_final.append(coords_final[i]- coords_final[i-1])
    v_initial.append(np.cross(v_initial[0],v_initial[1]))
    v_final.append(np.cross(v_final[0],v_final[1]))
    v_initial = np.array(v_initial).T
    v_final = np.array(v_final).T
    #print(v_final)
    rotation_matrix = np.matmul(v_final,np.linalg.inv(v_initial))
    beta = np.arctan2(np.sqrt(1-rotation_matrix[2,2]**2),(rotation_matrix[2,2]))
    if np.abs(np.sin(beta))>0.0001:
        gamma = (np.arctan2(rotation_matrix[2,1],-(rotation_matrix[2,0])))
        alpha = np.arctan2(rotation_matrix[1,2],rotation_matrix[0,2])
    else:
        alpha = np.arctan2(rotation_matrix[1,0],rotation_matrix[0,0])
        gamma = 0
    return alpha,beta,gamma

def real_to_complex_spherical_harmonics(coeffs):
    
    l = (len(coeffs)-1)//2
    complex_coeffs = []
    if l==1:
        complex_coeffs.append(-1/np.sqrt(2)*(coeffs[0]+1j*coeffs[1]))
        complex_coeffs.append(1/np.sqrt(2)*(coeffs[0]-1j*coeffs[1]))
        complex_coeffs.append(coeffs[2])
    elif l>1:
        complex_coeffs.append(coeffs[0])
        for m_ind in range(1,2*l+1,2):
            complex_coeffs.append(1/np.sqrt(2)*(coeffs[m_ind]-1j*coeffs[m_ind + 1]))
            complex_coeffs.append(-1/np.sqrt(2)*(coeffs[m_ind]+1j*coeffs[m_ind + 1]))
    else:
        return coeffs
    return np.array(complex_coeffs)

def complex_to_real_spherical_harmonics(coeffs):
    ##### FIXXXXXX GENERALIZE THIS IS JUST FOR d !!!
    real_coeffs = []
    l = (len(coeffs)-1)//2
    if l == 1:
        real_coeffs.append(-np.sqrt(2)*np.real(coeffs[0]))
        real_coeffs.append(np.sqrt(2)*np.real(-1j*coeffs[0]))
        real_coeffs.append(np.real(coeffs[2]))
    elif l>1:
        real_coeffs.append(np.real(coeffs[0]))
        for m_ind in range(1,2*l+1,2):
            real_coeffs.append(-np.sqrt(2)*np.real(coeffs[m_ind]))
            real_coeffs.append(np.sqrt(2)*np.real(-1j*coeffs[m_ind]))
    else:
        return coeffs
    return np.array(real_coeffs)
#read_basis(os.path.join("WATER_TESTS","original","basis"))
def rotate_orbitals(
    initial_MO_coeffs,
    coords_initial,
    coords_final,
    azimuthal_quantum_number_list
):
    l_max = max(azimuthal_quantum_number_list)
    alpha,beta,gamma = calculate_euler_angles(coords_initial, coords_final)
    wigner_D_matrix_dict = {}
    for l in range(0,l_max + 1):
        wigner_D_matrix_dict[l] = create_wigner_D_matrix(l,alpha,beta,gamma)
    MO_array_index = 0
    final_MO_coeffs = np.zeros_like(initial_MO_coeffs)
    for l in azimuthal_quantum_number_list:
        if l == 0:
            final_MO_coeffs[MO_array_index] = initial_MO_coeffs[MO_array_index]
            MO_array_index += 1
        else:
            coeff_block = real_to_complex_spherical_harmonics(
                initial_MO_coeffs[MO_array_index:MO_array_index + 2 * l + 1,:]
            )
            final_MO_coeffs[MO_array_index:MO_array_index + 2 * l + 1,:] = (
                complex_to_real_spherical_harmonics(
                    np.matmul(
                        wigner_D_matrix_dict[l],
                        coeff_block
                        )
                )
            )
            MO_array_index += 2 * l + 1
            #print("INITIAL", initial_MO_coeffs[MO_array_index:MO_array_index + 2 * l + 1,:])
            #print("COEFF_BLOCK",coeff_block)
            #print(complex_to_real_spherical_harmonics(np.matmul(
            #            wigner_D_matrix_dict[l],
            #            coeff_block
            #            )))
            #exit()
    return final_MO_coeffs        
AO_basis = read_basis(os.path.join("WATER_TESTS","original","basis"))
atom_list,coords_initial = read_coord(os.path.join("WATER_TESTS","original","coord"))
atoms,coords_rotated = read_coord(os.path.join("WATER_TESTS","3rd_rotated","coord"))
initial_MO_coeffs = read_MOs(os.path.join("WATER_TESTS","original","mos"))
azimuthal_quantum_number_list = get_azimuthal_q_num_list(AO_basis, atom_list)
print(rotate_orbitals(
    initial_MO_coeffs,
    coords_initial,
    coords_rotated,
    azimuthal_quantum_number_list
)[:,0])
print("41")
print(rotate_orbitals(
    initial_MO_coeffs,
    coords_initial,
    coords_rotated,
    azimuthal_quantum_number_list
)[:,-1])
