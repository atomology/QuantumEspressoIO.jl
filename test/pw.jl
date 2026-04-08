@testitem "read_atomic_species!" begin
    using OrderedCollections: OrderedDict
    lines = [
        "ATOMIC_SPECIES",
        "! a comment line",
        "Si       28.085500  Si.upf",
        "O        15.999000  O.upf",
        "following line",
    ]
    lcopy = copy(lines)
    card = QuantumEspressoIO.read_atomic_species!(lcopy, 2)
    @test isa(card, Pair)
    name, content = card
    @test name == "atomic_species"
    expected = OrderedDict{String, Any}(
        "option" => nothing,
        "species" => ["Si", "O"],
        "masses" => [28.0855, 15.999],
        "pseudos" => ["Si.upf", "O.upf"],
    )
    @test content == expected
    @test lcopy == ["following line"]
end

@testitem "read_atomic_positions!" begin
    using OrderedCollections: OrderedDict
    lines = [
        "ATOMIC_POSITIONS crystal",
        "! a comment line",
        "Si        0.0000000000      0.0000000000      0.0000000000",
        "O         0.5000000000      0.5000000000      0.5000000000",
        "following line",
    ]
    lcopy = copy(lines)
    card = QuantumEspressoIO.read_atomic_positions!(lcopy, 2)
    @test isa(card, Pair)
    name, content = card
    @test name == "atomic_positions"
    expected = OrderedDict{String, Any}(
        "option" => "crystal",
        "atoms" => ["Si", "O"],
        "positions" => [Vec3(0.0, 0.0, 0.0), Vec3(0.5, 0.5, 0.5)],
    )
    @test content == expected
    @test lcopy == ["following line"]
end

@testitem "read_cell_parameters!" begin
    using OrderedCollections: OrderedDict
    lines = [
        "CELL_PARAMETERS angstrom",
        "! a comment line",
        "1.0 0.0 0.0",
        "0.0 2.0 0.0",
        "0.0 0.0 3.0",
        "following line",
    ]
    lcopy = copy(lines)
    card = QuantumEspressoIO.read_cell_parameters!(lcopy)
    @test isa(card, Pair)
    name, content = card
    @test name == "cell_parameters"
    expected = OrderedDict{String, Any}(
        "option" => "angstrom",
        "cell" => [Vec3(1.0, 0.0, 0.0), Vec3(0.0, 2.0, 0.0), Vec3(0.0, 0.0, 3.0)],
    )
    @test content == expected
    @test lcopy == ["following line"]
end

@testitem "read_k_points!" begin
    using OrderedCollections: OrderedDict
    lines = [
        "K_POINTS crystal",
        "! a comment line",
        "2",
        "0.0 0.0 0.0 1.0",
        "0.5 0.5 0.5 1.0",
        "following line",
    ]
    lcopy = copy(lines)
    card = QuantumEspressoIO.read_k_points!(lcopy)
    @test isa(card, Pair)
    name, content = card
    @test name == "k_points"
    expected = OrderedDict{String, Any}(
        "option" => "crystal",
        "kpoints" => [Vec3(0.0, 0.0, 0.0), Vec3(0.5, 0.5, 0.5)],
        "kweights" => [1.0, 1.0],
    )
    @test content == expected
    @test lcopy == ["following line"]
end

@testitem "read_pw_in" begin
    using OrderedCollections: OrderedDict
    io = IOBuffer(
        """&control
                calculation = "scf"
            /
            &system
                ibrav = 0
                nat = 2
                ntyp = 2
            /
            ATOMIC_SPECIES
            Si       28.085500  Si.upf
            O        15.999000  O.upf
            ATOMIC_POSITIONS crystal
            Si        0.0000000000      0.0000000000      0.0000000000
            O         0.5000000000      0.5000000000      0.5000000000
            CELL_PARAMETERS angstrom
            1.0 0.0 0.0
            0.0 2.0 0.0
            0.0 0.0 3.0
            K_POINTS crystal
            2
            0.0000000000      0.0000000000      0.0000000000      1.0000000000
            0.5000000000      0.5000000000      0.5000000000      1.0000000000
        """
    )
    params = QuantumEspressoIO.read_pw_in(io)
    @test isa(params, AbstractDict)
    expected = OrderedDict(
        "control" => OrderedDict("calculation" => "scf"),
        "system" => OrderedDict("ibrav" => 0, "nat" => 2, "ntyp" => 2),
        "atomic_species" => OrderedDict("option" => nothing, "species" => ["Si", "O"], "masses" => [28.0855, 15.999], "pseudos" => ["Si.upf", "O.upf"]),
        "atomic_positions" => OrderedDict("option" => "crystal", "atoms" => ["Si", "O"], "positions" => [Vec3(0.0, 0.0, 0.0), Vec3(0.5, 0.5, 0.5)]),
        "cell_parameters" => OrderedDict("option" => "angstrom", "cell" => [Vec3(1.0, 0.0, 0.0), Vec3(0.0, 2.0, 0.0), Vec3(0.0, 0.0, 3.0)]),
        "k_points" => OrderedDict("option" => "crystal", "kpoints" => [Vec3(0.0, 0.0, 0.0), Vec3(0.5, 0.5, 0.5)], "kweights" => [1.0, 1.0]),
    )
    @test params == expected
end

@testitem "write_atomic_species" begin
    using OrderedCollections: OrderedDict
    card = OrderedDict(
        "species" => ["Si", "O"],
        "masses" => [28.0855, 15.999],
        "pseudos" => ["Si.upf", "O.upf"],
    )
    buf = IOBuffer()
    QuantumEspressoIO.write_atomic_species(buf, card)
    out = String(take!(buf))
    expected = """ATOMIC_SPECIES
    Si       28.085500  Si.upf
    O        15.999000  O.upf
    """
    @test out == expected
end

@testitem "write_atomic_positions" begin
    using OrderedCollections: OrderedDict
    card_ap = OrderedDict(
        "option" => "crystal",
        "atoms" => ["Si", "O"],
        "positions" => [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    )
    buf = IOBuffer()
    QuantumEspressoIO.write_atomic_positions(buf, card_ap)
    out = String(take!(buf))
    expected = """ATOMIC_POSITIONS crystal
    Si        0.0000000000      0.0000000000      0.0000000000
    O         0.5000000000      0.5000000000      0.5000000000
    """
    @test out == expected
end

@testitem "write_cell_parameters" begin
    using OrderedCollections: OrderedDict
    card_cp = OrderedDict(
        "option" => "angstrom",
        "cell" => [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
    )
    buf = IOBuffer()
    QuantumEspressoIO.write_cell_parameters(buf, card_cp)
    out = String(take!(buf))
    expected = """CELL_PARAMETERS angstrom
        1.0000000000      0.0000000000      0.0000000000
        0.0000000000      2.0000000000      0.0000000000
        0.0000000000      0.0000000000      3.0000000000
    """
    @test out == expected
end

@testitem "write_k_points (crystal)" begin
    using OrderedCollections: OrderedDict
    card_k = OrderedDict(
        "option" => "crystal",
        "kpoints" => [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        "kweights" => [1.0, 1.0],
    )
    buf = IOBuffer()
    QuantumEspressoIO.write_k_points(buf, card_k)
    out = String(take!(buf))
    expected = """K_POINTS crystal
    2
        0.0000000000      0.0000000000      0.0000000000      1.0000000000
        0.5000000000      0.5000000000      0.5000000000      1.0000000000
    """
    @test out == expected
end

@testitem "write_k_points (automatic)" begin
    using OrderedCollections: OrderedDict
    card_ka = OrderedDict(
        "option" => "automatic",
        "kgrid" => [8, 8, 8],
        "kgrid_shift" => [0, 1, 1],
    )
    buf = IOBuffer()
    QuantumEspressoIO.write_k_points(buf, card_ka)
    out = String(take!(buf))
    expected = """K_POINTS automatic
    8 8 8    0 1 1
    """
    @test out == expected
end

@testitem "write_pw_in (combined)" begin
    using OrderedCollections: OrderedDict
    inputs = OrderedDict(
        "control" => OrderedDict("calculation" => "scf", "prefix" => "SiO", "outdir" => "./out", "pseudo_dir" => "./pseudo"),
        "system" => OrderedDict("ibrav" => 0, "nat" => 2, "ntyp" => 2),
        "electrons" => OrderedDict("conv_thr" => 1.0e-6),
        "atomic_species" => OrderedDict("species" => ["Si", "O"], "masses" => [28.0855, 15.999], "pseudos" => ["Si.upf", "O.upf"]),
        "atomic_positions" => OrderedDict("option" => "crystal", "atoms" => ["Si", "O"], "positions" => [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "cell_parameters" => OrderedDict("option" => "angstrom", "cell" => [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]),
        "k_points" => OrderedDict("option" => "crystal", "kpoints" => [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], "kweights" => [1.0, 1.0]),
    )
    buf = IOBuffer()
    write_pw_in(buf, inputs)
    out = String(take!(buf))
    expected = """&control
      calculation = 'scf'
      prefix = 'SiO'
      outdir = './out'
      pseudo_dir = './pseudo'
    /
    &system
      ibrav = 0
      nat = 2
      ntyp = 2
    /
    &electrons
      conv_thr = 1.0e-6
    /
    ATOMIC_SPECIES
    Si       28.085500  Si.upf
    O        15.999000  O.upf
    ATOMIC_POSITIONS crystal
    Si        0.0000000000      0.0000000000      0.0000000000
    O         0.5000000000      0.5000000000      0.5000000000
    CELL_PARAMETERS angstrom
        1.0000000000      0.0000000000      0.0000000000
        0.0000000000      2.0000000000      0.0000000000
        0.0000000000      0.0000000000      3.0000000000
    K_POINTS crystal
    2
        0.0000000000      0.0000000000      0.0000000000      1.0000000000
        0.5000000000      0.5000000000      0.5000000000      1.0000000000
    """
    @test out == expected
end
