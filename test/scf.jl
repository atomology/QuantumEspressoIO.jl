using QuantumEspressoIO

# params = QuantumEspressoIO.parse_qe_in("scf.in")
params = QuantumEspressoIO.parse_qe_in("/Users/SashaP/Desktop/ge/scf_sc.in")

println("Parsed parameters:")
println(params)


QuantumEspressoIO.write_qe_in("scf_tst.in", params)