using Printf

export write_namelist, write_qe_in

"""
    $(SIGNATURES)

Write a Fortran namelist.

# Arguments
- `io::IO`: The IO stream to write to.
- `name::StrOrSym`: The name of the namelist.
- `namelist::AbstractDict`: The namelist to write, which is a dictionary-like object.

# Examples
```jldoctest
using QuantumEspressoIO

name = :input
namelist = Dict(
    "a" => 1,
    "b" => 2.0,
    "c" => [3, 4, 5],
    "d" => "test",
    "e" => true,
)
write_namelist(stdout, name, namelist)
# output
&input
  c(1) = 3
  c(2) = 4
  c(3) = 5
  e = .true.
  b = 2.0
  a = 1
  d = 'test'
/
```
"""
function write_namelist(io::IO, name::StrOrSym, namelist::AbstractDict)
    println(io, "&" * string(name))
    indent = 2
    for (k, v) in namelist
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
- `namelist::AbstractDict`: The namelist to write, which is a dictionary-like object.
"""
function write_namelist(filename::AbstractString, name::StrOrSym, namelist::AbstractDict)
    return open(filename, "w") do io
        write_namelist(io, name, namelist)
    end
end

"""
    $(SIGNATURES)

Generic writer for QE input namelist.

# Arguments
- `io::IO`: The IO stream to write to.
- `inputs::AbstractDict`: The input parameters, each key-value pair is treated as a namelist.

# Examples
```jldoctest
using QuantumEspressoIO
using OrderedCollections

# Use OrderedDict to preserve the order of the keys
inputs = OrderedDict(
    :control => OrderedDict("calculation" => "scf", "prefix" => "qe"),
    :system => OrderedDict("ecutwfc" => 30.0, "ecutrho" => 300.0),
    :electrons => OrderedDict("mixing_beta" => 0.7),
)
write_qe_in(stdout, inputs)
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
function write_qe_in(io::IO, inputs::AbstractDict)
    for (name, namelist) in pairs(inputs)
        write_namelist(io, name, namelist)
    end
end

"""
    $(SIGNATURES)

Generic writer for QE input namelist to a file.

# Arguments
- `filename::AbstractString`: The name of the file to write to.
- `inputs::AbstractDict`: The input parameters, see [`write_qe_in(io::IO, inputs::AbstractDict)`](@ref) for details.
"""
function write_qe_in(filename::AbstractString, inputs::AbstractDict)
    return open(filename, "w") do io
        write_qe_in(io, inputs)
    end
end
