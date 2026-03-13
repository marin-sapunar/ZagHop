""" Physical unit conversion factors from atomic units.

Atomic units are the default throughout ZagHop. Each constant equals
one atomic unit expressed in the target unit, so converting *from*
atomic units is a multiplication::

    energy_eV = energy_au * eV

and converting *to* atomic units is a division::

    energy_au = energy_eV / eV

Values follow CODATA 2018 recommended constants.
"""

# Energy
eV = 27.211386245988        # Hartree -> electronvolt
cm1 = 219474.6313702        # Hartree -> cm^-1
kcalmol = 627.5094740631    # Hartree -> kcal/mol
kelvin = 315775.02480407    # Hartree -> Kelvin  (E_h / k_B)

ENERGY = {
    "Hartree": 1.0,
    "eV": eV,
    "cm^-1": cm1,
    "kcal/mol": kcalmol,
    "K": kelvin,
}

# Distance
angstrom = 0.529177210903   # Bohr -> Angstrom
angstrom_inv = 1 / angstrom # Angstrom -> Bohr

# Mass
dalton = 1822.888486209     # Dalton (u) -> electron mass

# Time
fs = 0.02418884326509       # atomic time unit -> femtosecond
ps = fs * 1e-3              # atomic time unit -> picosecond

# Temperature / Boltzmann
boltzmann_au = 1 / kelvin   # k_B in Hartree / Kelvin
