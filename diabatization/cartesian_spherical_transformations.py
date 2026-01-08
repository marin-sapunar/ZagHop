import numpy as np
from scipy.special import factorial,binom
import os

def gaussian_to_spherical_coeff(
l,
m,
l_x,
l_y,
l_z
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

l_tot = 3
m=1
triplets = [(i, j, l_tot - i - j) 
            for i in range(l_tot + 1) 
            for j in range(l_tot + 1 - i)]
for l_mu in triplets:
    print(gaussian_to_spherical_coeff(l_tot,m,*l_mu),l_mu)
