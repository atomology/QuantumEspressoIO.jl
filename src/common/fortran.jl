"""
Convert a Julia value to a string representation of Fortran value.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: format_fortran)
julia> format_fortran("hello")
"'hello'"
julia> format_fortran(1.0)
"1.0"
julia> format_fortran(true)
".true."
```
"""
function format_fortran(value)
    if isa(value, AbstractString)
        return "'$value'"
    elseif isa(value, Bool)
        return value ? ".true." : ".false."
    elseif isa(value, Number)
        return string(value)
    else
        error("Unsupported type: $(typeof(value))")
    end
end

"""
    $(SIGNATURES)

Parse a string as `Float64`.

The is capable of parsing Fortran outputs, e.g. `1.0D-10`, to the ordinary `1e-10`.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: parse_float)
julia> parse_float("1.0D-10")
1.0e-10
julia> parse_float("1***")
NaN
```
"""
function parse_float(s::AbstractString)
    if occursin("*", s)
        return NaN
    else
        return parse(Float64, replace(lowercase(strip(s)), "d" => "e"))
    end
end

"""
    $(SIGNATURES)

Parse a string as `bool`.

This is capable of parsing Fortran outputs, e.g., `.true.`, `.false.`, `true`, `T`.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: parse_bool)
julia> parse_bool(".true.")
true
julia> parse_bool("false")
false
julia> parse_bool("T")
true
julia> parse_bool("F")
false
julia> parse_bool("1")
true
julia> parse_bool("0")
false
```
"""
function parse_bool(s::AbstractString)
    s = replace(lowercase(strip(s)), "." => "")[1]  # only 1st char
    return s == 't' || s == '1'
end

"""
    $(SIGNATURES)

Parse an integer as `bool`.

- `0`: `false`
- `1` or `-1`: `true`

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: parse_bool)
julia> parse_bool(1)
true
julia> parse_bool(0)
false
```
"""
function parse_bool(i::Integer)
    return i != 0
end


"""
    $(SIGNATURES)

Parse a Fortran value.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: parse_value)
julia> parse_value("'hello'")
"hello"
julia> parse_value(".true.")
true
julia> parse_value("1.0D-10")
1.0e-10
julia> parse_value("1")
1
```
"""
function parse_value(value::AbstractString)
    # string
    str_signs = ["'", '"']
    for s in str_signs
        if startswith(value, s) && endswith(value, s)
            return value[2:end-1]
        end
    end

    # bool
    bool_signs = ["t", "f", "true", "false"]
    for s in bool_signs
        value2 = replace(lowercase(strip(value)), "." => "")
        if value2 == s
            return parse_bool(s)
        end
    end

    # integer, use regex
    if occursin(r"^\s*[-+]?\d+\s*$", value)
        return parse(Int, value)
    end

    # float, use regex
    if occursin(r"^\s*[-+]?\d+(\.\d*)?([DdEe][+-]?\d+)?\s*$", value)
        return parse_float(value)
    end
end

"""
    $(SIGNATURES)

Remove comment and strip empty spaces.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: remove_comment)
julia> remove_comment("  ! This is a comment")
""
julia> remove_comment("  input")
"input"
julia> remove_comment("  name1 = value1  ! comment 2")
"name1 = value1"
julia> remove_comment("  name2 = value2  # comment 3")
"name2 = value2"
```
"""
function remove_comment(line::AbstractString)
    comment_signs = ["!", "#"]

    line = strip(line)
    if any(startswith(line, s) for s in comment_signs)
        return ""
    end

    # This also works for inline comments
    for s in comment_signs
        if occursin(s, line)
            return strip(split(line, s; limit=2)[1])
        end
    end

    return line
end
