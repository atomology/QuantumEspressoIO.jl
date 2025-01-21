#Inspired from DFControl.jl

const NeedleType = Union{AbstractString, AbstractChar, Regex}

function parse_qe_in(path_to_scf::String)
    scf_parameters = Dict{Symbol, Any}()#:format => "espresso-in", :crystal_coordinates => true
    parse_file(path_to_scf, QE_PW_PARSE_FUNCTIONS, out = scf_parameters)

    return scf_parameters
end

function getfirst(f::Function, A)
    for el in A
        if f(el)
            return el
        end
    end
end

function parse_file(f::AbstractString, parse_funcs::Vector{<:Pair{NeedleType, Any}}; out = Dict{Symbol, Any}())
    lc = 0
    open(f, "r") do file
        while !eof(file)
            line = strip(readline(file))
            lc += 1
            if isempty(line)
                continue
            end
            # Iterate over the pairs in QE_PW_PARSE_FUNCTIONS to find matching lines
            for pf in parse_funcs
                if occursin(pf.first, line)
                    # try
                        entry = pf.second(line, file)
                        #need to append to the diciony out[entry.block]  with entry.flag and entry.result
                        if entry.block in keys(out)
                            out[entry.block][entry.flag] = entry.result
                        else
                            out[entry.block] = Dict{Symbol, Any}(entry.flag => entry.result)
                        end
                        # if entry.flag in keys(out[entry.block])
                        #     out[entry.block][entry.flag] = entry.result
                        # else
                        #     out[entry.block][entry.flag] = entry.result
                        # end
                       
                    # catch e
                        # @warn "File corruption or parsing error detected executing parse function $(pf.second) in file $f at line $lc: \"$line\".\nTrying to continue smoothly. Error: $e"
                    # end
                    break  # Exit the loop once the matching function is found
                end
            end
        end
    end
    return out
end

function parse_file(f::AbstractString, args...; kwargs...)
    open(f, "r") do file
        parse_file(file, args...;kwargs...)
    end
end

function qe_parse_calculation(line, f)
    block = :CONTROL
    flag = :calculation
    result = strip(split(line)[3], ['\''])
    return (;block, flag, result)
end

function qe_parse_verbosity(line, f)
    block = :CONTROL
    flag = :verbosity
    result = strip(split(line)[3], ['\''])
    return (;block, flag, result)
end

function qe_parse_tstress(line, f)
    block = :CONTROL
    flag = :tsress
    result = split(line)[3] == ".true."
    return (;block, flag, result)
end

function qe_parse_tprnfor(line, f)
    block = :CONTROL
    flag = :tprnfor
    result = split(line)[3] == ".true."
    return (;block, flag, result)
end

function qe_parse_outdir(line, f)
    block = :CONTROL
    flag = :outdir
    result = strip(split(line)[3], ['\''])
    return (;block, flag, result)
end

function qe_parse_prefix(line, f)
    block = :CONTROL
    flag = :tsress
    result = strip(split(line)[3], ['\''])
    return (;block, flag, result)
end

function qe_parse_pseudo_dir(line, f)
    block = :CONTROL
    flag = :pseudo_dir
    result = strip(split(line)[3], ['\''])
    return (;block, flag, result)
end

function qe_parse_ibrav(line, f)
    block = :SYSTEM
    flag = :ibrav
    result = parse(Int, split(line)[3])
    return (;block, flag, result)
end

function qe_parse_nbnd(line, f)
    block = :SYSTEM
    flag = :nbnd
    result = parse(Int, split(line)[3])
    return (;block, flag, result)
end

function qe_parse_ecutwfc(line, f)
    block = :SYSTEM
    flag = :ecutwfc
    result = parse(Float64, split(line)[3])
    return (;block, flag, result)
end

function qe_parse_ecutrho(line, f)
    block = :SYSTEM
    flag = :ecutrho
    result = parse(Float64, split(line)[3])
    return (;block, flag, result)
end

function qe_parse_nosym(line, f)
    block = :SYSTEM
    flag = :nosym
    result = split(line)[3] == ".true."
    return (;block, flag, result)
end

function qe_parse_noinv(line, f)
    block = :SYSTEM
    flag = :noinv
    result = split(line)[3] == ".true."
    return (;block, flag, result)
end

function qe_parse_nat(line, f)
    block = :SYSTEM
    flag = :nat
    result = parse(Int, split(line)[3])
    return (;block, flag, result)
end

function qe_parse_diagonalization(line, f)
    block = :ELECTRONS
    flag = :diagonalization
    result = strip(split(line)[3], ['\''])
    return (;block, flag, result)
end

function qe_parse_electrons_maxstep(line, f)
    block = :ELECTRONS
    flag = :electron_maxstep
    result = parse(Int, split(line)[3])
    return (;block, flag, result)
end

function qe_parse_mixing_mode(line, f)
    block = :ELECTRONS
    flag = :mixing_mode
    result = strip(split(line)[3], ['\''])
    return (;block, flag, result)
end

function qe_parse_mixing_beta(line, f)
    block = :ELECTRONS
    flag = :mixing_beta
    result = parse(Float64, split(line)[3])
    return (;block, flag, result)
end

function qe_parse_conv_thr(line, f)
    block = :ELECTRONS
    flag = :conv_thr
    result = parse(Float64, split(line)[3])
    return (;block, flag, result)
end

function qe_parse_kpoints(line, f)
    block = :K_POINTS
    flag = :kpts
    
    line = readline(f)
    values = [parse(Int,val) for val in  split(line)]
    #case of uniform unshifted grid for now
    result = (values[1], values[2], values[3], values[4], values[5], values[6])
    return (;block, flag, result)
end

function qe_parse_species(line, f)
    block = :ATOMIC_SPECIES
    flag = :species

    result = []
    # line = readline(f)
    while !eof(f)
        line = strip(readline(f))
        length(split(line)) != 3 && break
        species = split(line)
        symb = species[1]
        mass = species[2]
        pseudo = species[3]
        push!(result, (;symb, mass, pseudo))
    end

    return (;block, flag, result)
end

function qe_parse_atomic_positions(line, f)
    block = :ATOMIC_POSITIONS
    flag = :positions

    result = []
    # line = readline(f)
    while !eof(f)
        line = strip(readline(f))
        length(split(line)) != 4 && break
        species = split(line)
        symb = species[1]
        x = species[2]
        y = species[3]
        z = species[4]
        push!(result, (;symb, x, y, z))
    end

    return (;block, flag, result)
end

const QE_PW_PARSE_FUNCTIONS::Vector{Pair{NeedleType, Any}}  = [
    #Contol
    "calculation" => qe_parse_calculation,
    "verbosity" => qe_parse_verbosity,
    "tstress" => qe_parse_tstress,
    "tprnfor" => qe_parse_tprnfor,
    "outdir" => qe_parse_outdir,
    "prefix" => qe_parse_prefix,
    "pseudo_dir" => qe_parse_pseudo_dir,
    #System
    "ibrav" => qe_parse_ibrav,
    "nbnd" => qe_parse_nbnd,
    "ecutwfc" => qe_parse_ecutwfc,
    "ecutrho" => qe_parse_ecutrho,
    "nosym" => qe_parse_nosym,
    "noinv" => qe_parse_noinv,
    "nat" => qe_parse_nat,
    #Electrons
    "diagonalization" => qe_parse_diagonalization,
    "electrons_maxstep" => qe_parse_electrons_maxstep,
    "mixing_mode" => qe_parse_mixing_mode,
    "mixing_beta" => qe_parse_mixing_beta,
    "conv_thr" => qe_parse_conv_thr,
    #Atomic_species
    "ATOMIC_SPECIES" => qe_parse_species,
    #K-points
    "K_POINTS" => qe_parse_kpoints,
    # ATOMIC_POSITIONS
    "ATOMIC_POSITIONS" => qe_parse_atomic_positions,
]


# Trying to write scf.in based on input dicitonary

function write_qe_in(path_to_scf::String, scf_parameters::Dict{Symbol, Any})
    order = [:CONTROL, :SYSTEM, :ELECTRONS]

    open(path_to_scf, "w") do file
        for block in order
            write(file, "&$block \n")
            for (key, value) in scf_parameters[block]
                if key in [:format, :crystal_coordinates]
                    continue
                end
                write(file, "  $key = ")
                if isa(value, AbstractString)
                    write(file, "'$value'")
                elseif isa(value, Bool)
                    write(file, value ? ".true." : ".false.")
                elseif isa(value, Tuple)
                    write(file, value)
                else
                    write(file, "$value")
                end
                write(file, "\n")
            end
            write(file, "/\n\n")
        end

        write(file, "ATOMIC_SPECIES\n")
        for (_, values) in scf_parameters[:ATOMIC_SPECIES]
            for value in values
                write(file, "$(value.symb)  $(value.mass)  $(value.pseudo) \n")
            end
        end
        write(file, "\n")

        write(file, "K_POINTS automatic\n")
        write(file, "  ")
        for value in scf_parameters[:K_POINTS][:kpts]
            write(file, "$value ")
        end
        write(file, "\n \n")

        write(file, "ATOMIC_POSITIONS crystal\n")
        for value in scf_parameters[:ATOMIC_POSITIONS][:positions]
            write(file, "$(value.symb)  $(value.x)  $(value.y)  $(value.z) \n")
        end

    end
end
