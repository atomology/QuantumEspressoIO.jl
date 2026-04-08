@testitem "read_band_dat" begin
    io = IOBuffer(
        """&plot nbnd=  11, nks=  1 /
        0.0  0.0  0.0
        -1.0  0.0  1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0
        9.0
        """
    )
    res = QuantumEspressoIO.read_band_dat(io)
    @test length(res.kpoints) == length(res.eigenvalues) == 1
    @test res.kpoints[1] == Vec3(0.0, 0.0, 0.0)
    @test res.eigenvalues[1] == [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

    # guess_high_symmetry_kpoints
    kpoints = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.1, 0.1, 0.0], [0.1, 0.1, 0.1]]
end

@testitem "guess_high_symmetry_kpoints" begin
    kpoints = [
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [0.2, 0.0, 0.0],
        [0.2, 0.1, 0.0],
        [0.2, 0.1, 0.1],
    ]
    indices = QuantumEspressoIO.guess_high_symmetry_kpoints(kpoints)
    @test indices == [1, 3, 4, 5]
end
