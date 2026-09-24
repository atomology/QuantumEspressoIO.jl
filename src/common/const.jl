
"""
Bohr radius in Angstrom unit.

This is the default (Physical constants, SI (NIST 2018)) value in QE
`Modules/constants.f90`.
"""
const Bohr_QE::Float64 = 0.529177210903

"""
pw.x input keywords
"""
const PW_KEYWORDS = Set([
    "atomic_species",
    "atomic_positions",
    "k_points",
    "additional_k_points",
    "cell_parameters",
    "constraints",
    "occupations",
    "atomic_velocities",
    "atomic_forces",
    "solvents",
    "hubbard",
])
