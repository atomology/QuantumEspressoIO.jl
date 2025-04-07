using QuantumEspressoIO
# using LazyArtifacts
# using Artifacts
# using Test

@testset "parse scf.in" begin
    rootpath = artifact"Si"
    path_tst_data = joinpath(rootpath, "scf.in")
    params = QuantumEspressoIO.parse_qe_in(path_tst_data)

    params_tst = Dict(
                    :SYSTEM => Dict(
                        :nat => 2,
                        :ecutwfc => 100.0,
                        :ntyp => 1,
                        :ibrav => 0
                    ),
                    :ELECTRONS => Dict(
                        :mixing_beta => 0.3,
                        :conv_thr => 1.0e-12,
                        :mixing_mode => "plain",
                        :diagonalization => "david"
                    ),
                    :K_POINTS => Dict(
                        :kpts => [8, 8, 8, 0, 0, 0]
                    ),
                    :CELL_PARAMETERS => Dict(
                        :cell => [-2.71526 0.0 2.71526;
                                0.0 2.71526 2.71526;
                                -2.71526 2.71526 0.0]
                    ),
                    :CONTROL => Dict(
                        :prefix => "silicon",
                        :tstress => true,
                        :tprnfor => true,
                        :verbosity => "high",
                        :calculation => "scf",
                        :outdir => "./tmp/",
                        :pseudo_dir => "pseudo/"
                    ),
                    :ATOMIC_POSITIONS => Dict(
                        :positions => [
                            (symb = "Si", x = 0.0, y = 0.0, z = 0.0),
                            (symb = "Si", x = 0.75, y = 0.75, z = 0.75)
                        ],
                        :format => "crystal"
                    ),
                    :ATOMIC_SPECIES => Dict(
                        :species => [
                            (symb = "Si", mass = "28.085", pseudo = "Si.upf")
                        ]
                    ),
    )

    @test params == params_tst
end

@testset "write scf.in" begin
    rootpath = artifact"Si"
    path_tst_data = joinpath(rootpath, "scf.in")
    params = QuantumEspressoIO.parse_qe_in(path_tst_data)
    QuantumEspressoIO.write_qe_in("scf_tst.in", params)
    params_tst = QuantumEspressoIO.parse_qe_in("scf_tst.in")
    rm("scf_tst.in")

    @test params == params_tst
end
