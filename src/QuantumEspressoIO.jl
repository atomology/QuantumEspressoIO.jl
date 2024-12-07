module QuantumEspressoIO

# Write your package code here.

include("common/const.jl")
include("common/type.jl")

using FortranFiles

include("scf.jl")
include("xml.jl")
include("bin.jl")

end
