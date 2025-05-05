@testitem "Test parsing wavefunctions from binaries of QE" begin
    using LazyArtifacts
    thr = 1e-13
    rootpath = artifact"Si"
    path_tst_data = joinpath(rootpath, "wfc1.dat")

    _, evc_list = QuantumEspressoIO.read_wfc_dat(path_tst_data)
    # calculate_braket(bra, ket) = sum(conj(bra[i]) * ket[i] for i in eachindex(bra))
    function calculate_braket(bra, ket)
        result = sum(conj(bra[i]) * ket[i] for i in eachindex(bra))
        return result
    end

    norm11 = abs(calculate_braket(evc_list[1], evc_list[1]))
    norm12 = abs(calculate_braket(evc_list[1], evc_list[2]))

    @test isapprox(norm11, 1.0; atol=thr)
    @test isapprox(norm12, 0.0; atol=thr)
end
