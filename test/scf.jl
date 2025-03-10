using QuantumEspressoIO

using LazyArtifacts
using Artifacts
rootpath = artifact"Si"
path_tst_data = joinpath(rootpath, "scf.in")

#TODO make proper test
params = QuantumEspressoIO.parse_qe_in(path_tst_data)

println("Parsed parameters:")
println(params)

QuantumEspressoIO.write_qe_in("scf_tst.in", params)
