export write_pw_in

"""
    $(SIGNATURES)

Write the `atomic_species` card of `pw.x`.

# Examples
```jldoctest
using QuantumEspressoIO: write_atomic_species

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
```jldoctest
using QuantumEspressoIO: write_atomic_positions

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
```jldoctest
using QuantumEspressoIO: write_cell_parameters

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
```jldoctest
using QuantumEspressoIO: write_k_points

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
        @printf(io, "%d %d %d    %f %f %f\n", card[:kgrid]..., card[:kgrid_shift]...)
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
using QuantumEspressoIO
using OrderedCollections

# Use OrderedDict to preserve the order of the keys
inputs = OrderedDict(
    :control => OrderedDict(
        :calculation => "scf",
        :prefix => "SiO2",
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
