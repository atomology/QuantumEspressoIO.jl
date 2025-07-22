using OrderedCollections
using DelimitedFiles

"""
    $(SIGNATURES)

Read the `projwfc.x` output `projwfc_up` file.

# Returns
- `params`: Dictionary of parameters
- `orbitals`: the projected orbitals, each element is a NamedTuple with fields:
    - `atom_index`: index of the atom
    - `atom_label`: label of the atom
    - `label`: label of the orbital, e.g., "3S"
    - `n`: principal quantum number
    - `l`: azimuthal quantum number
    - `m`: magnetic quantum number
- `projections`: the projection data, size: `n_kpoint * n_bands * n_orbitals`
    - `n_kpoints`: number of kpoints
    - `n_bands`: number of bands
    - `n_orbitals`: number of projected orbitals
"""
function read_projwfc_up(io::IO)
    splitline() = split(strip(readline(io)))

    # The variable names are strange, as they mirror the names in the QE Fortran code.
    # TODO maybe we should invent better names?
    params = OrderedDict{String, Any}()
    # header
    title = strip(readline(io), '\n')
    params["title"] = title

    nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp = parse.(Int, splitline())
    params["nr1x"] = nr1x
    params["nr2x"] = nr2x
    params["nr3x"] = nr3x
    params["nr1"] = nr1
    params["nr2"] = nr2
    params["nr3"] = nr3

    line = splitline()
    ibrav = parse(Int, line[1])
    params["ibrav"] = ibrav

    celldm = parse.(Float64, line[2:end])
    params["celldm"] = celldm

    line = splitline()
    # some version of projwfc.x output the unit_cell
    if length(line) == 3
        readline(io)
        readline(io)
        line = splitline()
    end
    gcutm, dual, ecutwfc = parse.(Float64, line[1:(end - 1)])
    params["gcutm"] = gcutm
    params["dual"] = dual
    params["ecutwfc"] = ecutwfc

    magicnum = parse(Int, line[end])
    (magicnum == 9) || error("magic number mismatch")

    atm = Vector{String}(undef, ntyp)
    zv = zeros(Float64, ntyp)
    for i in 1:ntyp
        line = splitline()
        nt = parse(Int, line[1])
        (nt == i) || error("atom type mismatch")
        atm[i] = line[2]
        zv[i] = parse(Float64, line[3])
    end
    params["atm"] = atm
    params["zv"] = zv

    tau = zeros(Float64, 3, nat)
    ityp = zeros(Int, nat)
    for i in 1:nat
        line = splitline()
        na = parse(Int, line[1])
        (na == i) || error("atom index mismatch")
        tau[:, i] = parse.(Float64, line[2:4])
        ityp[i] = parse(Int, line[5])
    end
    params["tau"] = tau
    params["ityp"] = ityp

    natomwfc, nkstot, nbnd = parse.(Int, splitline())
    params["natomwfc"] = natomwfc
    params["nkstot"] = nkstot
    params["nbnd"] = nbnd

    noncolin, lspinorb = parse_bool.(splitline())
    # TODO support noncolin and lspinorb
    (!noncolin && !lspinorb) || error("noncolin and lspinorb not supported")
    params["noncolin"] = noncolin
    params["lspinorb"] = lspinorb

    # projection data
    nlmchi = Vector{NamedTuple}()
    proj = zeros(Float64, nkstot, nbnd, natomwfc)
    for iw in 1:natomwfc
        line = splitline()
        nwfc = parse(Int, line[1])
        (nwfc == iw) || error("wfc index mismatch")
        na = parse(Int, line[2])
        atm_name = line[3]
        (atm_name == atm[ityp[na]]) || error("atom name mismatch")
        els = line[4]
        n, l, m = parse.(Int, line[5:end])
        push!(nlmchi, (; na, els, n, l, m))
        for ik in 1:nkstot
            for ib in 1:nbnd
                line = splitline()
                k, b = parse.(Int, line[1:2])
                (k == ik && b == ib) || error("kpt and band index mismatch")
                p = parse(Float64, line[3])
                proj[ik, ib, iw] = p
            end
        end
    end
    params["nlmchi"] = nlmchi

    # Construct a nicer list of orbital names
    orbitals = Vector{NamedTuple}(undef, natomwfc)
    for i in 1:natomwfc
        atom_index = nlmchi[i]["na"]
        atom_label = atm[ityp[atom_index]]
        label = nlmchi[i]["els"]
        n = nlmchi[i]["n"]
        l = nlmchi[i]["l"]
        m = nlmchi[i]["m"]
        orbitals[i] = (; atom_index, atom_label, label, n, l, m)
    end

    # TODO add test
    return params, orbitals, proj
end

function read_projwfc_up(filename::AbstractString)
    return open(filename) do io
        read_projwfc_up(io)
    end
end

"""
    read_projwfc_dos(io)

Read `projwfc.x` total PDOS (e.g. `qe.pdos_tot`) or atom-projected file
(e.g. `qe.pdos_atm#1(Si)_wfc#1(s)`).

# Arguments
- `io` or `filename`: An IO stream or a filename to read the DOS data from.

# Returns
- `columns::Matrix{Float64}`: A matrix containing the DOS data.
- `header::Vector{String}`: A vector containing the first line of the file.
"""
function read_projwfc_dos(io::IO)
    header = readline(io)
    header = strip(chopprefix(header, "#"))
    EeV = "E (eV)"
    if startswith(header, EeV)
        header = chopprefix(header, EeV)
        header = strip(header)
    else
        error("Expected header to start with '$EeV', got '$header'")
    end
    header = String[EeV, split(header)...]
    columns = readdlm(io, Float64; comments=true)
    return columns, header
end

function read_projwfc_dos(filename::AbstractString)
    return open(filename) do io
        read_projwfc_dos(io)
    end
end

"""
    read_pdos_tot(io)

Read `projwfc.x` total PDOS file.

# Arguments
- `prefix::AbstractString`: The prefix of the PDOS file (e.g., `qe` for `qe.pdos_tot`).

# Returns
- `columns::Matrix{Float64}`: A matrix containing the PDOS data.
- `header::Vector{String}`: A vector containing the first line of the file.
"""
function read_pdos_tot(prefix::AbstractString)
    filename = "$prefix.pdos_tot"
    return read_projwfc_dos(filename)
end

"""
    read_pdos_atm(prefix::AbstractString)

Read all the atom-projected PDOS files.

# Arguments
- `prefix::AbstractString`: The prefix of the PDOS files
    (e.g., `qe` for `qe.pdos_atm#1(Si)_wfc#1(s)`).

# Returns
- PDOS as a vector, which is ordered according to the atom index. Each element
  of the vector is a tuple with the following fields:
  - `atom_index::Int`: Index of the atom (1-based).
  - `atom_label::String`: Label of the atom.
  - `pdos`: A vector of named tuples, each containing:
    - `wfc_index::Int`: Index of the wavefunction.
    - `wfc_label::String`: Label of the wavefunction.
    - `columns::Matrix{Float64}`: The PDOS data matrix.
    - `header::Vector{String}`: The header of the PDOS file.
"""
function read_pdos_atm(prefix::AbstractString)
    dir, name = dirname(prefix), basename(prefix)
    isempty(dir) && (dir = ".")
    files = readdir(dir)
    regex = Regex(name * raw"\.pdos_atm#(\d+)\((.+)\)_wfc#(\d+)\((.+)\)")
    # files = filter(Base.Fix2(startswith, "$name.pdos_atm"), files)
    files = filter(f -> occursin(regex, f), files)
    atom_indices = map(files) do f
        parse(Int, match(regex, f).captures[1])
    end
    natoms = length(unique(atom_indices))
    res = map(1:natoms) do i
        files_i = files[atom_indices .== i]
        atom_label = match(regex, files_i[1]).captures[2]
        wfc_indices = map(f -> parse(Int, match(regex, f).captures[3]), files_i)
        nwfcs_i = length(files_i)
        res = map(1:nwfcs_i) do j
            file = only(files_i[wfc_indices .== j])
            wfc_label = match(regex, file).captures[4]
            # Each pdos file has the same format as pdos_tot file
            columns, header = read_projwfc_dos(joinpath(dir, file))
            return (; wfc_index=j, wfc_label, columns, header)
        end
        return (; atom_index=i, atom_label, pdos=res)
    end
    return res
end
