import numpy as np
from scipy.special import factorial,binom
import os
#import cartesian_spherical_transformations
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
def powers_to_array(power_string):
    # molden ordering of powers (separated by ', ') converted to numpy array. TEMPORARY SOLUTION???
    arr = power_string.split(", ")
    power_array = []
    for el in arr:
        power_array.append(np.array(list(el),dtype = float))
    return np.array(power_array)
def molden_cartesian_GTO_order():
    """
    Element [l] of returned dictionary shows order of coefficients
    for Cartesian GTOs for the total angular momentum l=l_x+l_y+l_z
    """
    # TEMPORARILY WITHOUT GENERATING FORMULA???
    
    GTO_power_dict = {}
    GTO_power_dict[0] = powers_to_array("000")
    GTO_power_dict[1] = powers_to_array("100, 010, 001")
    GTO_power_dict[2] = powers_to_array("200, 020, 002, 110, 101, 011")
    GTO_power_dict[3] = powers_to_array("300, 030, 003, 120, 210, 201, 102, 012, 021, 111")
    GTO_power_dict[4] = powers_to_array("400, 040, 004, 310, 301, 130, 031, 103, 013, 220, 202, 022, 211, 121, 112")
    return GTO_power_dict
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
        complex_coeffs.append(-1/np.sqrt(2)*(coeffs[0]-1j*coeffs[1]))
        complex_coeffs.append(coeffs[2])
    elif l>1:
        complex_coeffs.append(coeffs[0])
        complex_coeffs.append(-1/np.sqrt(2)*(coeffs[1]-1j*coeffs[2]))
        complex_coeffs.append(1/np.sqrt(2)*(coeffs[1]+1j*coeffs[2]))

        complex_coeffs.append(1/np.sqrt(2)*(coeffs[4]+1j*coeffs[3]))

        complex_coeffs.append(1/np.sqrt(2)*(coeffs[4]-1j*coeffs[3]))
    else:
        return coeffs
    return np.array(complex_coeffs)
def complex_to_real_spherical_harmonics(coeffs):
    real_coeffs = []
    l = (len(coeffs)-1)//2
    if l == 1:
        real_coeffs.append(-np.sqrt(2)*np.real(coeffs[0]))
        real_coeffs.append(-np.sqrt(2)*np.real(-1j*coeffs[0]))
        real_coeffs.append(np.real(coeffs[2]))
    elif l>1:
        real_coeffs.append(np.real(coeffs[0]))
        #print("HIHI",coeffs[0])

        real_coeffs.append(-np.sqrt(2)*np.real(coeffs[1]))
        real_coeffs.append(-np.sqrt(2)*np.real(1j*coeffs[1]))

        real_coeffs.append(np.sqrt(2)*np.real(1j*coeffs[4]))
        real_coeffs.append(np.sqrt(2)*np.real(coeffs[4]))
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
def gaussian_to_spherical_coeff(
l,m,l_x,l_y,l_z
):
    """
    l = azimuthal_q_n
    m = magnetic_q_n
    l_x, l_y, l_z = Cartesian Gaussian
    """
    abs_m = abs(m)
    j = (l_x + l_y - abs_m)/2
    if j != int(j):
        return 0
    else:
        j=int(j)
    coeff_prefactor = np.sqrt(
        factorial(2 * l_x) * factorial(2 * l_y) * factorial(2 * l_z) *
         factorial(l-abs_m) /
        (
        factorial(2 * l) * factorial(l_x) * factorial(l_y) *
        factorial(l_z) * factorial(l + abs_m) * factorial(l)
        )
    ) /(2**l) 
    # Factorial(l) is in denominator instead sqrt(l!) up and
    # 1/l! outside the root - in case something needs to be fixed
    if m>0:
        coeff_prefactor *= (-1)**(abs_m)
    sum_1 = 0
    for i in range((l-abs_m)//2+1):
        sum_1 += (binom(l,i) * binom(i,j) * (-1)**i *
            factorial(2 * l - 2 * i) /
            factorial(l - abs_m - 2 * i)
        )
    sum_2 = 0 
    for k in range(j+1):
        sum_2 += (binom(j,k) * binom(abs_m, l_x - 2*k) *
        (-1 + 0J )**(np.sign(m)*(abs_m - l_x + 2*k)/2)
        )
    return coeff_prefactor * sum_1 * sum_2

def cartesian_to_spherical_transformation_matrix(l):
    molden_GTO_order_dict = molden_cartesian_GTO_order()
    GTO_ordering = molden_GTO_order_dict[l]
    transformation_matrix = np.zeros([2*l+1,(l+1)*(l+2)//2],dtype = complex)
    #spherical_GTO_m_ordering = []
    for m_index in range(2*l + 1):
        for GTO_index in range((l+1)*(l+2)//2):
            m = index_to_m(m_index,l)
            l_x,l_y,l_z = GTO_ordering[GTO_index]
            transformation_matrix[m_index,GTO_index] = (
                gaussian_to_spherical_coeff(l,m,l_x,l_y,l_z)
            )
            
    return transformation_matrix
# TESTING PART

#transformation_matrix = cartesian_to_spherical_transformation_matrix(2)
#print(transformation_matrix)
#
##AO_basis = read_basis(os.path.join("WATER_TESTS","original","basis"))
##atom_list,coords_initial = read_coord(os.path.join("WATER_TESTS","original","coord"))
##initial_MO_coeffs = read_MOs(os.path.join("WATER_TESTS","original","mos"))
#
## CHECKING PHASE
#"""
#AO_basis = read_basis(os.path.join("WATER_TESTS","rotated","basis"))
#atom_list,coords_initial = read_coord(os.path.join("WATER_TESTS","rotated","coord"))
#initial_MO_coeffs = read_MOs(os.path.join("WATER_TESTS","rotated","mos"))
#atoms,coords_rotated = read_coord(os.path.join("WATER_TESTS","3rd_rotated","coord"))
#ROTATED_MO_coeffs = read_MOs(os.path.join("WATER_TESTS","3rd_rotated","mos"))
#"""
#"""
#AO_basis = read_basis(os.path.join("WATER_TESTS","3rd_rotated","basis"))
#atom_list,coords_initial = read_coord(os.path.join("WATER_TESTS","3rd_rotated","coord"))
#initial_MO_coeffs = read_MOs(os.path.join("WATER_TESTS","3rd_rotated","mos"))
#atoms,coords_rotated = read_coord(os.path.join("WATER_TESTS","2nd_rotated","coord"))
#ROTATED_MO_coeffs = read_MOs(os.path.join("WATER_TESTS","2nd_rotated","mos"))
#"""
#AO_basis = read_basis(os.path.join("WATER_TESTS","3rd_rotated","basis"))
#atom_list,coords_initial = read_coord(os.path.join("WATER_TESTS","3rd_rotated","coord"))
#initial_MO_coeffs = read_MOs(os.path.join("WATER_TESTS","3rd_rotated","mos"))
#atoms,coords_rotated = read_coord(os.path.join("WATER_TESTS","ANOTHER_ACTUAL_ROTATION","coord"))
#ROTATED_MO_coeffs = read_MOs(os.path.join("WATER_TESTS","ANOTHER_ACTUAL_ROTATION","mos"))
#azimuthal_quantum_number_list = get_azimuthal_q_num_list(AO_basis, atom_list)
#COEFFS_SPHE_HARM = rotate_orbitals(
#    initial_MO_coeffs,
#    coords_initial,
#    coords_rotated,
#    azimuthal_quantum_number_list
#)[:,0]
##print("41")
##print(rotate_orbitals(
##    initial_MO_coeffs,
##    coords_initial,
##    coords_rotated,
##    azimuthal_quantum_number_list
##)[:,-1])
## 13 0.35208720192995E-03
#d_ORBS_GAUSSIAN_STRING  = np.array("""
#    14 0.36466688672864E-04
#    15 -.14648409405826E-04
#    16 -.21818279267038E-04
#    17 -.68963261488800E-05
#    18 -.29307987727516E-04
#    19 -.25761806553887E-04
#    20 -.11154121028243E-03
#    21 0.48087658203949E-04
#    22 0.63453552078483E-04
#    23 0.17017962033446E-04
#    24 0.88519404928930E-04
#    25 0.68705698372418E-04
#    """.split(),dtype = float).reshape(-1,2)
#d_ORBS_GAUSSIAN_COEFFS = d_ORBS_GAUSSIAN_STRING[:,1]
#d_ORBS_1_GAUSS = d_ORBS_GAUSSIAN_COEFFS[0:6]
#d_ORBS_2_GAUSS = d_ORBS_GAUSSIAN_COEFFS[6:]
##print(d_ORBS_2_GAUSS)
#print("OG",COEFFS_SPHE_HARM[13:18])
#print("OG",ROTATED_MO_coeffs[:,0][13:18])
##print("ROTATED",ROTATED_MO_coeffs[:,0])
##d_ORBS_SPHE_HARM_1 = COEFFS_SPHE_HARM[13:18]
#d_ORBS_SPHE_HARM_1 = COEFFS_SPHE_HARM[13:18]
##print(d_ORBS_1_GAUSS)
##print("OG",d_ORBS_SPHE_HARM_1)
##print(transformation_matrix)
#d_ORBS_SPHE_HARM_1 = ROTATED_MO_coeffs[:,0][13:18]
##print("TRANSFORMATION",transformation_matrix.T[4])
#print("HERE",np.matmul(transformation_matrix.T,real_to_complex_spherical_harmonics(d_ORBS_SPHE_HARM_1))/np.sqrt(3))
#VORW_TRANS = np.matmul(transformation_matrix.T,real_to_complex_spherical_harmonics(d_ORBS_SPHE_HARM_1))/np.sqrt(3)
#rev_tranf = np.linalg.pinv(transformation_matrix.T)
##print(complex_to_real_spherical_harmonics(np.matmul(rev_tranf,d_ORBS_1_GAUSS))*np.sqrt(3))#/np.sqrt(2*2-1))
##print(complex_to_real_spherical_harmonics(np.matmul(rev_tranf,VORW_TRANS))*np.sqrt(3))#/np.sqrt(2*2-1))
#np.savetxt("TRANSFORMATION.dat",transformation_matrix.T)
#hmm_trans = np.linalg.pinv(transformation_matrix)
##print("HERE",np.matmul(hmm_trans,real_to_complex_spherical_harmonics(d_ORBS_SPHE_HARM_1))/np.sqrt(3))
#
