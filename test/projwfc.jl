@testitem "read_projwfc_up" begin
    using LazyArtifacts

    f = joinpath(artifact"Si2", "si2.projwfc_up")
    params, orbitals, proj = QuantumEspressoIO.read_projwfc_up(f)

    @test params["atm"] == ["Si"]
    nkpts, nbands, nprojs = 511, 16, 8
    @test params["natomwfc"] == nprojs
    @test params["nkstot"] == nkpts
    @test params["nbnd"] == nbands

    @test orbitals == [
        (atom_index=1, atom_label="Si", label="3S", n=1, l=0, m=1)
        (atom_index=1, atom_label="Si", label="3P", n=2, l=1, m=1)
        (atom_index=1, atom_label="Si", label="3P", n=2, l=1, m=2)
        (atom_index=1, atom_label="Si", label="3P", n=2, l=1, m=3)
        (atom_index=2, atom_label="Si", label="3S", n=1, l=0, m=1)
        (atom_index=2, atom_label="Si", label="3P", n=2, l=1, m=1)
        (atom_index=2, atom_label="Si", label="3P", n=2, l=1, m=2)
        (atom_index=2, atom_label="Si", label="3P", n=2, l=1, m=3)
    ]

    @test size(proj) == (nkpts, nbands, nprojs)
    @test proj[2, 1:5, 1] ≈ [0.4981375412, 0.0005278635, 0.0, 0.0, 0.0001201347]
end

@testitem "read_pdos_tot" begin
    using LazyArtifacts

    f = joinpath(artifact"Si2", "si2")
    columns, header = QuantumEspressoIO.read_pdos_tot(f)

    @test size(columns) == (180, 3)
    @test header == ["E (eV)", "dos(E)", "pdos(E)"]
    @test columns[[1, end], :] ≈ [
        -6.234 1.93e-6 1.93e-6
        29.566 -5.32e-5 -7.6e-7
    ]
end

@testitem "read_pdos_atm" begin
    using LazyArtifacts

    f = joinpath(artifact"Si2", "si2")
    pdos = QuantumEspressoIO.read_pdos_atm(f)

    @test length(pdos) == 2
    @test pdos[1].atom_index == 1
    @test pdos[1].atom_label == "Si"
    @test length(pdos[1].pdos) == 2
    @test pdos[1].pdos[1].wfc_label == "s"
    @test pdos[1].pdos[1].columns[[1, end], :] ≈ [
        -6.234   9.58e-7    9.58e-7
        29.566  -9.59e-12  -9.59e-12
    ]
end
