export write_matdyn_in

"""
    $(SIGNATURES)

Write the `matdyn.x` input file.

# Arguments
- `io::IO`: The IO stream to write to.
- `inputs::AbstractDict`: The input data, containing the following fields:
    - `:input`: An `OrderedDict{Symbol,Any}` containing the namelist parameters.
    - `:qpoints`: A vector of q-points in crystal coordinates, each being a
        vector of three floats.

# Examples
```jldoctest
# Use OrderedDict to preserve the order of the keys
using OrderedCollections

inputs = OrderedDict(
    :input => OrderedDict(
        :flfrc => "q2r.fc",
        :flfrq => "matdyn.freq",
        :q_in_cryst_coord => true,
    ),
    :qpoints => [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
    ],
)
write_matdyn_in(stdout, inputs)
# output
&input
  flfrc = 'q2r.fc'
  flfrq = 'matdyn.freq'
  q_in_cryst_coord = .true.
/
2
    0.0000000000      0.0000000000      0.0000000000
    0.5000000000      0.5000000000      0.5000000000
```
"""
function write_matdyn_in(io::IO, inputs::AbstractDict)
    write_namelist(io, :input, inputs[:input])

    # qpoints in reciprocal space
    nq = length(inputs[:qpoints])
    println(io, string(nq))
    for q in inputs[:qpoints]
        @printf(io, "%16.10f  %16.10f  %16.10f\n", q...)
    end
end

"""
    $(SIGNATURES)

Write the `matdyn.x` input file to a file.

# Arguments
- `filename::AbstractString`: The name of the file to write to.
- `inputs::AbstractDict`: See [`write_matdyn_in(io::IO, inputs::AbstractDict)`](@ref) for details.
"""
function write_matdyn_in(filename::AbstractString, inputs::AbstractDict)
    return open(filename, "w") do io
        write_matdyn_in(io, inputs)
    end
end
