"""
Convert a Julia value to a string representation of Fortran value.
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
