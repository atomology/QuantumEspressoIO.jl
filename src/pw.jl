using OrderedCollections

export read_pw_in, write_pw_in

"""
    $(SIGNATURES)

Parse the `atomic_species` card of `pw.x` input.

# Arguments
- `lines::AbstractVector`: The lines of the input file.
- `n_species::Integer`: The number of species in the `atomic_species` card.

# Returns
- A `Pair` of card name to card content. The card option is stored under
    the `:option` key in the card content.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: read_atomic_species!)
lines = [
    "ATOMIC_SPECIES",
    "! a comment line",
    "Si       28.085500  Si.upf",
    "O        15.999000  O.upf",
    "following line",
];
n_species = 2;
card = read_atomic_species!(lines, n_species)
println(card)
println(lines)
# output
:atomic_species => OrderedCollections.OrderedDict{Symbol, Any}(:option => nothing, :species => ["Si", "O"], :masses => [28.0855, 15.999], :pseudos => ["Si.upf", "O.upf"])
["following line"]
```
"""
function read_atomic_species!(lines::AbstractVector, n_species::Integer)
    name = "atomic_species"
    result = parse_card!(lines, name, n_species)
    isnothing(result) && return nothing
    option, content = result

    species = String[]
    masses = Float64[]
    pseudos = String[]
    for line in content
        sp, ma, ps = split(line)
        push!(species, sp)
        push!(masses, parse_float(ma))
        push!(pseudos, ps)
    end

    card = OrderedDict{Symbol, Any}()
    card[:option] = option
    card[:species] = species
    card[:masses] = masses
    card[:pseudos] = pseudos

    return Symbol(name) => card
end

"""
    $(SIGNATURES)

Parse the `atomic_positions` card of `pw.x` input.

# Arguments
- `lines::AbstractVector`: The lines of the input file.
- `n_atoms::Integer`: The number of atoms in the `atomic_positions` card.

# Returns
- A `Pair` of card name to card content. The card option is stored under
    the `:option` key in the card content.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: read_atomic_positions!)
lines = [
    "ATOMIC_POSITIONS crystal",
    "! a comment line",
    "Si        0.0000000000      0.0000000000      0.0000000000",
    "O         0.5000000000      0.5000000000      0.5000000000",
    "following line",
]
n_atoms = 2
card = read_atomic_positions!(lines, n_atoms)
println(card)
println(lines)
# output
:atomic_positions => OrderedCollections.OrderedDict{Symbol, Any}(:option => "crystal", :atoms => ["Si", "O"], :positions => StaticArraysCore.SVector{3, Float64}[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
["following line"]
```
"""
function read_atomic_positions!(lines::AbstractVector, n_atoms::Integer)
    name = "atomic_positions"
    result = parse_card!(lines, name, n_atoms)
    isnothing(result) && return nothing
    option, content = result

    atoms = String[]
    positions = Vec3{Float64}[]
    for line in content
        at, x, y, z = split(line)
        push!(atoms, at)
        push!(positions, Vec3(parse_float.([x, y, z])))
    end

    card = OrderedDict{Symbol, Any}()
    card[:option] = option
    card[:atoms] = atoms
    card[:positions] = positions

    return Symbol(name) => card
end

"""
    $(SIGNATURES)

Parse the `cell_parameters` card of `pw.x` input.

# Arguments
- `lines::AbstractVector`: The lines of the input file.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: read_cell_parameters!)
lines = [
    "CELL_PARAMETERS angstrom",
    "! a comment line",
    "1.0 0.0 0.0",
    "0.0 2.0 0.0",
    "0.0 0.0 3.0",
    "following line",
]
card = read_cell_parameters!(lines)
println(card)
println(lines)
# output
:cell_parameters => OrderedCollections.OrderedDict{Symbol, Any}(:option => "angstrom", :cell => StaticArraysCore.SVector{3, Float64}[[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
["following line"]
```
"""
function read_cell_parameters!(lines::AbstractVector)
    name = "cell_parameters"
    result = parse_card!(lines, name, 3)
    isnothing(result) && return nothing
    option, content = result

    cell = Vec3{Float64}[]
    for line in content
        push!(cell, Vec3(parse_float.(split(line))))
    end

    card = OrderedDict{Symbol, Any}()
    card[:option] = option
    card[:cell] = cell

    return Symbol(name) => card
end

"""
    $(SIGNATURES)

Parse the `k_points` card of `pw.x` input.

# Arguments
- `lines::AbstractVector`: The lines of the input file.

# Returns
- A `Pair` of card name to card content. The card option is stored under
    the `:option` key in the card content.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: read_k_points!)
lines = [
    "K_POINTS crystal",
    "! a comment line",
    "2",
    "0.0 0.0 0.0 1.0",
    "0.5 0.5 0.5 1.0",
    "following line",
]
card = read_k_points!(lines)
println(card)
println(lines)
# output
:k_points => OrderedCollections.OrderedDict{Symbol, Any}(:option => "crystal", :kpoints => StaticArraysCore.SVector{3, Float64}[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], :kweights => [1.0, 1.0])
["following line"]
```
"""
function read_k_points!(lines::AbstractVector)
    name = "k_points"
    # kpoints has variable number of lines,
    # the number of lines is indicated by the card option
    icard = find_card(lines, name)
    # no kpoints card found
    isnothing(icard) && return nothing

    option = parse_card_option(lines[icard])
    loption = lowercase(option)

    card = OrderedDict{Symbol, Any}()
    card[:option] = option

    # Possible options from pw.x input description:
    # tpiba | automatic | crystal | gamma | tpiba_b | crystal_b | tpiba_c | crystal_c
    if loption in ["tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c"]
        # kpoints in reciprocal space
        # 1st parse the next line after card name to get the number of kpoints
        # I use a deepcopy to avoid modifying the original lines
        nkpts = parse_card!(deepcopy(lines[icard:end]), name, 1)[2][1]
        nkpts = parse(Int, nkpts)
        n_lines = nkpts + 1
        content = parse_card!(lines, name, n_lines)[2]
        deleteat!(content, 1)
        kpoints = Vec3{Float64}[]
        kweights = Float64[]
        for line in content
            parts = parse_float.(split(line))
            push!(kpoints, Vec3(parts[1:3]))
            push!(kweights, parts[4])
        end
        card[:kpoints] = kpoints
        card[:kweights] = kweights
    elseif loption == "automatic"
        # automatic kpoints
        n_lines = 1
        content = parse_card!(lines, name, n_lines)[2]
        parts = split(content[1])
        kgrid = parse.(Int, parts[1:3])
        kgrid_shift = parse_float.(parts[4:6])
        card[:kgrid] = kgrid
        card[:kgrid_shift] = kgrid_shift
    elseif loption == "gamma"
        # gamma point, remove the card line
        n_lines = 0
        parse_card!(lines, name, n_lines)
    else
        # unknown option
        error("Unknown kpoints option: $option")
    end

    return Symbol(name) => card
end

"""
    $(SIGNATURES)

Read the `pw.x` input file.

# Arguments
- `io::Union{IO,AbstractString}`: The IO stream or filename to read from.

# Returns
- A dictionary of namelists and cards. The keys are the names of the namelists or cards.

# Examples
```jldoctest
io = IOBuffer(\"""
    &control
        calculation = "scf"
    /
    &system
        ibrav = 0
        nat = 2
        ntyp = 2
    /
    ATOMIC_SPECIES
    Si       28.085500  Si.upf
    O        15.999000  O.upf
    ATOMIC_POSITIONS crystal
    Si        0.0000000000      0.0000000000      0.0000000000
    O         0.5000000000      0.5000000000      0.5000000000
    CELL_PARAMETERS angstrom
    1.0 0.0 0.0
    0.0 2.0 0.0
    0.0 0.0 3.0
    K_POINTS crystal
    2
    0.0000000000      0.0000000000      0.0000000000      1.0000000000
    0.5000000000      0.5000000000      0.5000000000      1.0000000000
\""")
inputs = read_pw_in(io)
println(inputs)
# output
OrderedCollections.OrderedDict{Symbol, Any}(:control => OrderedCollections.OrderedDict{Symbol, Any}(:calculation => "scf"), :system => OrderedCollections.OrderedDict{Symbol, Any}(:ibrav => 0, :nat => 2, :ntyp => 2), :atomic_species => OrderedCollections.OrderedDict{Symbol, Any}(:option => nothing, :species => ["Si", "O"], :masses => [28.0855, 15.999], :pseudos => ["Si.upf", "O.upf"]), :atomic_positions => OrderedCollections.OrderedDict{Symbol, Any}(:option => "crystal", :atoms => ["Si", "O"], :positions => StaticArraysCore.SVector{3, Float64}[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]), :cell_parameters => OrderedCollections.OrderedDict{Symbol, Any}(:option => "angstrom", :cell => StaticArraysCore.SVector{3, Float64}[[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]), :k_points => OrderedCollections.OrderedDict{Symbol, Any}(:option => "crystal", :kpoints => StaticArraysCore.SVector{3, Float64}[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], :kweights => [1.0, 1.0]))
```
"""
function read_pw_in(io::Union{IO,AbstractString})
    # Parse the namelists
    namelists, cards = read_namelists(io; all_lines=true)

    # There are required parameters, also needed for parsing cards
    isnothing(get(namelists, :system, nothing)) && error("Missing namelist: system")
    isnothing(get(namelists[:system], :ntyp, nothing)) && error("Missing parameter: ntyp")
    n_species = namelists[:system][:ntyp]
    isnothing(get(namelists[:system], :nat, nothing)) && error("Missing parameter: nat")
    n_atoms = namelists[:system][:nat]

    # merge the cards into namelists
    params = namelists

    card = read_atomic_species!(cards, n_species)
    isnothing(card) || push!(params, card)

    card = read_atomic_positions!(cards, n_atoms)
    isnothing(card) || push!(params, card)

    card = read_cell_parameters!(cards)
    isnothing(card) || push!(params, card)

    card = read_k_points!(cards)
    isnothing(card) || push!(params, card)

    length(cards) == 0 || @error "Unrecognized cards in the input file: $cards"

    # TODO consider https://github.com/mcmcgrath13/DotMaps.jl for easy access
    return params
end

"""
    $(SIGNATURES)

Write the `atomic_species` card of `pw.x`.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: write_atomic_species)
inputs = Dict(
    :species => ["Si", "O"],
    :masses => [28.0855, 15.999],
    :pseudos => ["Si.upf", "O.upf"],
)
write_atomic_species(stdout, inputs)
# output
ATOMIC_SPECIES
Si       28.085500  Si.upf
O        15.999000  O.upf
```
"""
function write_atomic_species(io::IO, card::AbstractDict)
    println(io, "ATOMIC_SPECIES")
    for (sp, ma, ps) in zip(card[:species], card[:masses], card[:pseudos])
        @printf(io, "%-4s  %12.6f  %s\n", sp, ma, ps)
    end
end

"""
    $(SIGNATURES)

Write the `atomic_positions` card of `pw.x`.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: write_atomic_positions)
inputs = Dict(
    :option => "crystal",
    :atoms => ["Si", "O"],
    :positions => [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
    ],
)
write_atomic_positions(stdout, inputs)
# output
ATOMIC_POSITIONS crystal
Si        0.0000000000      0.0000000000      0.0000000000
O         0.5000000000      0.5000000000      0.5000000000
```
"""
function write_atomic_positions(io::IO, card::AbstractDict)
    option = get(card, :option, "")
    isempty(option) || (option = " $option")
    println(io, "ATOMIC_POSITIONS$option")
    for (sp, pos) in zip(card[:atoms], card[:positions])
        @printf(io, "%-4s  %16.10f  %16.10f  %16.10f\n", sp, pos...)
    end
end

"""
    $(SIGNATURES)

Write the `cell_parameters` card of `pw.x`.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: write_cell_parameters)
inputs = Dict(
    :option => "angstrom",
    :cell => [
        [1.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 3.0],
    ],
)
write_cell_parameters(stdout, inputs)
# output
CELL_PARAMETERS angstrom
    1.0000000000      0.0000000000      0.0000000000
    0.0000000000      2.0000000000      0.0000000000
    0.0000000000      0.0000000000      3.0000000000
```
"""
function write_cell_parameters(io::IO, card::AbstractDict)
    option = get(card, :option, "")
    isempty(option) || (option = " $option")
    println(io, "CELL_PARAMETERS$option")
    for v in card[:cell]
        @printf(io, "%16.10f  %16.10f  %16.10f\n", v...)
    end
end

"""
    $(SIGNATURES)

Write the `k_points` card of `pw.x`.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: write_k_points)
inputs = Dict(
    :option => "crystal",
    :kpoints => [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
    ],
    :kweights => [1.0, 1.0],
)
write_k_points(stdout, inputs)
# output
K_POINTS crystal
2
    0.0000000000      0.0000000000      0.0000000000      1.0000000000
    0.5000000000      0.5000000000      0.5000000000      1.0000000000
```
"""
function write_k_points(io::IO, card::AbstractDict)
    option = get(card, :option, "")
    loption = lowercase(option)

    valid_options_list = ["tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c"]
    valid_options_auto = ["automatic"]
    valid_options_gamma = ["gamma"]
    valid_options = vcat(valid_options_list, valid_options_auto, valid_options_gamma)
    any(option in valid_options) || error("Unknown kpoints option: $option")

    isempty(option) || (option = " $option")
    println(io, "K_POINTS$option")

    if loption in valid_options_list
        # kpoints in reciprocal space
        nkpts = length(card[:kpoints])
        println(io, string(nkpts))
        for (k, w) in zip(card[:kpoints], card[:kweights])
            @printf(io, "%16.10f  %16.10f  %16.10f  %16.10f\n", k..., w)
        end
    elseif loption in valid_options_auto
        # automatic kpoints
        @printf(io, "%d %d %d    %d %d %d\n", card[:kgrid]..., card[:kgrid_shift]...)
    elseif loption in valid_options_gamma
        # gamma point, nothing to do
    else
        # unknown option
        error("Unknown kpoints option: $option")
    end
end

"""
    $(SIGNATURES)

Write the `pw.x` input file.

# Arguments
- `io::IO`: The IO stream to write to.
- `inputs::AbstractDict`: The input data:
    - The keys are the names of the namelists or cards.
    - The values are the corresponding data for the namelists or cards.
    - The namelists are written first, followed by the cards.

# Examples
```jldoctest
# Use OrderedDict to preserve the order of the keys
using OrderedCollections

inputs = OrderedDict(
    :control => OrderedDict(
        :calculation => "scf",
        :prefix => "SiO",
        :outdir => "./out",
        :pseudo_dir => "./pseudo",
    ),
    :system => OrderedDict(
        :ibrav => 0,
        :nat => 2,
        :ntyp => 2,
    ),
    :electrons => OrderedDict(
        :conv_thr => 1e-6,
    ),
    :atomic_species => Dict(
        :species => ["Si", "O"],
        :masses => [28.0855, 15.999],
        :pseudos => ["Si.upf", "O.upf"],
    ),
    :atomic_positions => Dict(
        :option => "crystal",
        :atoms => ["Si", "O"],
        :positions => [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
        ],
    ),
    :cell_parameters => Dict(
        :option => "angstrom",
        :cell => [
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ],
    ),
    :k_points => Dict(
        :option => "crystal",
        :kpoints => [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
        ],
        :kweights => [1.0, 1.0],
    ),
)
write_pw_in(stdout, inputs)
# output
&control
  calculation = 'scf'
  prefix = 'SiO'
  outdir = './out'
  pseudo_dir = './pseudo'
/
&system
  ibrav = 0
  nat = 2
  ntyp = 2
/
&electrons
  conv_thr = 1.0e-6
/
ATOMIC_SPECIES
Si       28.085500  Si.upf
O        15.999000  O.upf
ATOMIC_POSITIONS crystal
Si        0.0000000000      0.0000000000      0.0000000000
O         0.5000000000      0.5000000000      0.5000000000
CELL_PARAMETERS angstrom
    1.0000000000      0.0000000000      0.0000000000
    0.0000000000      2.0000000000      0.0000000000
    0.0000000000      0.0000000000      3.0000000000
K_POINTS crystal
2
    0.0000000000      0.0000000000      0.0000000000      1.0000000000
    0.5000000000      0.5000000000      0.5000000000      1.0000000000
```
"""
function write_pw_in(io::IO, inputs::AbstractDict)
    valid_namelists = [:control, :system, :electrons, :ions, :cell, :fcp, :rism,]
    valid_cards = OrderedDict{Symbol,Function}(
        :atomic_species => write_atomic_species,
        :atomic_positions => write_atomic_positions,
        :cell_parameters => write_cell_parameters,
        :k_points => write_k_points,
        # :additional_k_points,
        # :constraints,
        # :occupations,
        # :atomic_velocities,
        # :atomic_forces,
        # :solvents,
        # :hubbard,
    )

    done_keys = Set{Symbol}()

    for name in valid_namelists
        if haskey(inputs, name)
            write_namelist(io, name, inputs[name])
            push!(done_keys, name)
        end
    end

    for (name, writer) in pairs(valid_cards)
        if haskey(inputs, name)
            writer(io, inputs[name])
            push!(done_keys, name)
        end
    end

    left_keys = setdiff(keys(inputs), done_keys)
    isempty(left_keys) || @error "Unrecognized keys in the input: $left_keys"
    return nothing
end

"""
    $(SIGNATURES)

Write the `pw.x` input file to a file.

# Arguments
- `filename::AbstractString`: The name of the file to write to.
- `inputs::AbstractDict`: See [`write_pw_in(io::IO, inputs::AbstractDict)`](@ref) for details.
"""
function write_pw_in(filename::AbstractString, inputs::AbstractDict)
    return open(filename, "w") do io
        write_pw_in(io, inputs)
    end
end

"""
Read relaxed structures from stdout file of `pw.x`.
"""
function read_structures(filename::AbstractString)
    return open(filename) do io
        natoms = -1
        cell_parameters = []
        atomic_positions = []

        while !eof(io)
            line = readline(io)
            if occursin("number of atoms/cell      =", line)
                natoms = parse(Int, strip(split(line, "=")[end]))
            elseif startswith(line, "CELL_PARAMETERS")
                lines = [line, [readline(io) for _ in 1:3]...]
                push!(cell_parameters, last(read_cell_parameters!(lines)))
            elseif startswith(line, "ATOMIC_POSITIONS")
                (natoms >= 0) || error("Expected `number of atoms/cell` line before ATOMIC_POSITIONS card")
                lines = [line, [readline(io) for _ in 1:natoms]...]
                card = last(read_atomic_positions!(lines, natoms))
                push!(atomic_positions, card)
            end
        end

        return cell_parameters, atomic_positions
    end
end
