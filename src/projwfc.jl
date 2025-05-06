using OrderedCollections

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
