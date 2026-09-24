using EzXML

"""
    $(SIGNATURES)

Read atomic structure and band structure from QE's XML output.

# Return
- `lattice`: `3 * 3`, Å, each column is a lattice vector
- `atom_positions`: length-`n_atoms` vector, each element is a fractional position
- `atom_labels`: length-`n_atoms` vector, each element is the label of the corresponding atom
- `recip_lattice`: `3 * 3`, Å⁻¹, each column is a reciprocal lattice vector
- `kpoints`: length-`n_kpts` vector, each element is a fractional kpoint
- `fermi_energy`: eV
- `alat`: the `alat` of QE in Å
- `eigenvalues`: `n_bands × n_kpts` matrix of eigenvalues in eV. For
    spin-polarized but without SOC calculations, return two matrices
    `eigenvalues_up` and `eigenvalues_dn` for the two spin channels.
"""
function read_pw_xml(filename::AbstractString)
    # from qe/Modules/constants.f90
    BOHR_RADIUS_ANGS = Bohr_QE  # Angstrom
    HARTREE_SI = 4.3597447222071e-18  # J
    ELECTRONVOLT_SI = 1.602176634e-19  # J
    AUTOEV = HARTREE_SI / ELECTRONVOLT_SI

    doc = readxml(filename)
    output = findfirst("/qes:espresso/output", root(doc))

    # atoms
    atomic_structure = findfirst("atomic_structure", output)
    alat = parse(Float64, atomic_structure["alat"])
    # from bohr to angstrom
    alat *= BOHR_RADIUS_ANGS
    n_atoms = parse(Int, atomic_structure["nat"])

    # structure info, each column is a vector for position or lattice vector
    atom_positions = Vec3{Float64}[]
    atom_labels = Vector{String}(undef, n_atoms)
    lattice = zeros(3, 3)

    for (i, atom) in enumerate(eachelement(findfirst("atomic_positions", atomic_structure)))
        pos = parse.(Float64, split(atom.content))
        push!(atom_positions, pos)
        atom_labels[i] = atom["name"]
    end
    # lattice
    for i in 1:3
        a = findfirst("cell/a$i", atomic_structure)
        lattice[:, i] = parse.(Float64, split(a.content))
    end
    # from cartesian to fractional
    inv_lattice = inv(lattice)
    atom_positions = map(atom_positions) do pos
        Vec3(inv_lattice * pos)
    end
    # from bohr to angstrom
    lattice *= BOHR_RADIUS_ANGS

    # reciprocal lattice
    recip_lattice = zeros(3, 3)
    for i in 1:3
        b = findfirst("basis_set/reciprocal_lattice/b$i", output)
        recip_lattice[:, i] = parse.(Float64, split(b.content))
    end
    # to 1/angstrom
    recip_lattice *= 2π / alat

    band_structure = findfirst("band_structure", output)
    n_kpts = parse(Int, findfirst("nks", band_structure).content)
    lsda = parse(Bool, findfirst("lsda", band_structure).content)
    # noncolin = parse(Bool, findfirst("noncolin", band_structure).content)
    spinorbit = parse(Bool, findfirst("spinorbit", band_structure).content)
    # check spin-polarized case
    if lsda && !spinorbit
        nbnd_up = parse(Int, findfirst("nbnd_up", band_structure).content)
        nbnd_dn = parse(Int, findfirst("nbnd_dw", band_structure).content)
        # they should be the same in QE
        @assert nbnd_up == nbnd_dn
        n_bands = nbnd_up
        eigenvalues_up = zeros(Float64, n_bands, n_kpts)
        eigenvalues_dn = zeros(Float64, n_bands, n_kpts)
    else
        n_bands = parse(Int, findfirst("nbnd", band_structure).content)
        eigenvalues = zeros(Float64, n_bands, n_kpts)
    end
    kpoints = Vec3{Float64}[]
    kweights = Float64[]

    n_electrons = parse(Float64, findfirst("nelec", band_structure).content)
    fermi_energy = parse(Float64, findfirst("fermi_energy", band_structure).content)
    # Hartree to eV
    fermi_energy *= AUTOEV

    inv_recip = inv(recip_lattice)
    ks_energies = findall("ks_energies", band_structure)
    length(ks_energies) == n_kpts || error("expected $n_kpts ks_energies, got $(length(ks_energies))")
    for (ik, ks_energy) in enumerate(ks_energies)
        k_point = findfirst("k_point", ks_energy)
        wt = parse(Float64, k_point["weight"])
        push!(kweights, wt)
        kpt = parse.(Float64, split(k_point.content))
        # to 1/angstrom
        kpt *= 2π / alat
        # from cartesian to fractional
        kpt = inv_recip * kpt
        push!(kpoints, kpt)

        qe_eigenvalues = findfirst("eigenvalues", ks_energy)
        if lsda && !spinorbit
            e = parse.(Float64, split(qe_eigenvalues.content))
            # Hartree to eV
            e .*= AUTOEV
            eigenvalues_up[:, ik] = e[1:n_bands]
            eigenvalues_dn[:, ik] = e[(n_bands + 1):end]
        else
            e = parse.(Float64, split(qe_eigenvalues.content))
            # Hartree to eV
            e .*= AUTOEV
            eigenvalues[:, ik] = e
        end
    end

    # Convert to `Mat3` matrices (each column is a lattice vector)
    lattice = mat3(lattice)
    recip_lattice = mat3(recip_lattice)

    results = (;
        lattice,
        atom_positions,
        atom_labels,
        recip_lattice,
        kpoints,
        kweights,
        n_electrons,
        fermi_energy,
        alat,
    )
    if lsda && !spinorbit
        return (; results..., eigenvalues_up, eigenvalues_dn)
    end
    return (; results..., eigenvalues)
end
