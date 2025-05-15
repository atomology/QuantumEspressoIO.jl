using Printf
using OrderedCollections

export read_namelist, write_namelist

"""
    $(SIGNATURES)

Read a fortran namelist file.

# Arguments
- `io::IO`: The IO stream to read from.

# Optional arguments 
- `output_remaining::Bool`: If `true`, return the remaining lines in the file. Default is `false`.

# Returns
- `namelists::OrderedDict`: A dictionary of namelists, each key is a symbol and
    the value is a `OrderedDict` of key-value pairs.
- `others::Vector`: A vector of strings, which are the remaining lines in the file.
    For example, QE `pw.x` input files may contain "cards" for additional parameters
    (e.g. atomic positions, k-points, etc.), whose syntax are customarily defined
    by `pw.x`.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: read_namelist)
io = IOBuffer(\"""
&input
    a = 1
    b = 2.0
    c = 'test'
    d = .true.
/
additional line
\""")
namelists, others = read_namelist(io)
# output
(OrderedCollections.OrderedDict{Symbol, Any}(:input => OrderedCollections.OrderedDict{Symbol, Any}(:a => 1, :b => 2.0, :c => "test", :d => true)), ["additional line"])
```
"""
function read_namelist(io::IO; output_remaining::Bool = false)
    namelists = OrderedDict{Symbol, Any}()
    others = String[]

    current_namelist = nothing

    for rawline in eachline(io)
        line = strip(rawline)
        # note this also remove inline comment
        line = remove_comment(line)
        isempty(line) && continue

        if startswith(line, "&")
            # Parse namelist
            # fortran is case-insensitive, so we convert to lowercase
            name = Symbol(lowercase(line[2:end]))
            namelists[name] = OrderedDict{Symbol, Any}()
            current_namelist = name
        elseif startswith(line, "/")
            # End of namelist
            current_namelist = nothing
            if length(line) > 1
                @warn "Content after / is ignored: $line"
            end
        elseif current_namelist !== nothing
            # Parse key-value pairs
            key_value_pairs = split(line, "="; limit=2)
            # fortran is case-insensitive, so we convert to lowercase
            key = Symbol(lowercase(strip(key_value_pairs[1])))
            # to be safe, we still keep the case of value unchanged
            value = strip(key_value_pairs[2])
            namelists[current_namelist][key] = parse_value(value)
        else
            # Remaining lines, keep raw lines
            push!(others, rawline)
        end
    end

    if output_remaining
        return namelists, others
    else
        return namelists
    end

end

function read_namelist(filename::AbstractString, output_remaining::Bool=false)
    return open(filename) do io
        read_namelist(io, output_remaining=output_remaining)
    end
end

"""
    $(SIGNATURES)

Find a card in the lines. Ignore comments and empty lines.

# Arguments
- `lines::AbstractVector`: The lines to read from.
- `name::AbstractString`: The name of the card to read.

# Returns
- The index of the card in the lines. If not found, return `nothing`.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: find_card)
lines = [
    "  ! This is a comment",
    "  input",
    "  name1 = value1  ! comment 2",
    "  name2 = value2",
]
name = "input"
find_card(lines, name)
# output
2
```
"""
function find_card(lines::AbstractVector, name::AbstractString)
    lname = lowercase(name)
    return findfirst(
        line -> startswith(lowercase(remove_comment(line)), lname),
        lines,
    )
end

"""
    $(SIGNATURES)

Get the option of a card.

# Arguments
- `line::AbstractString`: The line of the card.

# Returns
- The option of the card. If no option is found, return `nothing`.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: parse_card_option)
julia> parse_card_option("POSITIONS angstrom  ! comment")
"angstrom"
```
"""
@inline function parse_card_option(line::AbstractString)
    cardline = remove_comment(line)
    parts = split(cardline; limit=2)
    option =  length(parts) < 2 ? nothing : parts[2]
    return option
end

"""
    $(SIGNATURES)

From lines, get the card with name `name` and `n_lines` following lines.

# Arguments
- `lines`: The lines to read from.
- `name`: The name of the card to read.
- `n_lines`: The number of lines to read after the card name.
    Comment lines are ignored and not counted as `n_lines`.
    If `nothing`, read everything after the card.

# Example
```jldoctest; setup = :(using QuantumEspressoIO: parse_card!)
lines = [
    "  input option1",
    "  ! This is a comment",
    "  name1 = value1  ! comment 2",
    "  name2 = value2",
]
name = "input"
n_lines = 1
card = parse_card!(lines, name, n_lines)
println(card)
println(lines)
# output
("option1", ["name1 = value1"])
["  name2 = value2"]
```
"""
function parse_card!(lines::AbstractVector, name::AbstractString, n_lines::Union{Integer,Nothing}=nothing)
    istart = find_card(lines, name)
    # nothing found
    isnothing(istart) && return nothing

    option = parse_card_option(lines[istart])

    # remove any comments in the card
    content = String[]
    i = istart
    while length(content) < n_lines
        i += 1
        if i > length(lines)
            error("Not enough lines in the card")
        end
        line = remove_comment(lines[i])
        isempty(line) && continue
        push!(content, line)
    end
    iend = i

    deleteat!(lines, istart:iend)

    return option, content
end

"""
    $(SIGNATURES)

Write a Fortran namelist.

# Arguments
- `io::IO`: The IO stream to write to.
- `name::StrOrSym`: The name of the namelist.
- `params::AbstractDict`: The key-value pairs to write, which is a dictionary-like object.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: write_namelist)
using OrderedCollections

name = :input
params = OrderedDict(
    "a" => 1,
    "b" => 2.0,
    "c" => [3, 4, 5],
    "d" => "test",
    "e" => true,
)
write_namelist(stdout, name, params)
# output
&input
  a = 1
  b = 2.0
  c(1) = 3
  c(2) = 4
  c(3) = 5
  d = 'test'
  e = .true.
/
```
"""
function write_namelist(io::IO, name::StrOrSym, params::AbstractDict)
    println(io, "&" * string(name))
    indent = 2
    for (k, v) in params
        # convert to fortran values
        if isnothing(v)
            continue
        elseif isa(v, AbstractVector)
            for (i, vi) in enumerate(v)
                fvi = format_fortran(vi)
                @printf(io, "%s%s(%d) = %s\n", ' '^indent, k, i, fvi)
            end
        else
            fv = format_fortran(v)
            @printf(io, "%s%s = %s\n", ' '^indent, k, fv)
        end
    end
    println(io, "/")
end

"""
    $(SIGNATURES)

Write a Fortran namelist to a file.

# Arguments
- `filename::AbstractString`: The name of the file to write to.
- `name::StrOrSym`: The name of the namelist.
- `params::AbstractDict`: The key-value pairs to write, which is a dictionary-like object.
"""
function write_namelist(filename::AbstractString, name::StrOrSym, params::AbstractDict)
    return open(filename, "w") do io
        write_namelist(io, name, params)
    end
end

function write_namelist(io::Union{IO,AbstractString}, namelist::Pair)
    name = first(namelist)
    params = second(namelist)
    write_namelist(io, name, params)
end

"""
    $(SIGNATURES)

Generic writer for QE input namelist.

# Arguments
- `io::IO`: The IO stream to write to.
- `inputs::AbstractDict`: The input parameters, each key-value pair is treated as a namelist.

# Examples
```jldoctest
using OrderedCollections

# Use OrderedDict to preserve the order of the keys
inputs = OrderedDict(
    :control => OrderedDict("calculation" => "scf", "prefix" => "qe"),
    :system => OrderedDict("ecutwfc" => 30.0, "ecutrho" => 300.0),
    :electrons => OrderedDict("mixing_beta" => 0.7),
)
write_namelist(stdout, inputs)
# output
&control
  calculation = 'scf'
  prefix = 'qe'
/
&system
  ecutwfc = 30.0
  ecutrho = 300.0
/
&electrons
  mixing_beta = 0.7
/
```
"""
function write_namelist(io::IO, inputs::AbstractDict)
    for (name, namelist) in pairs(inputs)
        write_namelist(io, name, namelist)
    end
end

"""
    $(SIGNATURES)

Generic writer for QE input namelist to a file.

# Arguments
- `filename::AbstractString`: The name of the file to write to.
- `inputs::AbstractDict`: The input parameters, see [`write_namelist(io::IO, inputs::AbstractDict)`](@ref) for details.
"""
function write_namelist(filename::AbstractString, inputs::AbstractDict)
    return open(filename, "w") do io
        write_namelist(io, inputs)
    end
end
