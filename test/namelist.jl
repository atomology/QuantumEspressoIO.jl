@testitem "read_namelists" begin
    using OrderedCollections: OrderedDict
    io = IOBuffer(
        """&input
            a = 1
            b = 2.0
            c = 'test'
            d = .true.
        /
        additional line
        """
    )
    namelists, others = read_namelists(io; all_lines = true)
    @test namelists == OrderedDict("input" => OrderedDict("a" => 1, "b" => 2.0, "c" => "test", "d" => true))
    @test others == ["additional line"]
end

@testitem "read_namelist" begin
    using OrderedCollections: OrderedDict
    io = IOBuffer(
        """&input
            a = 1
            b = 2.0
            c = 'test'
            d = .true.
        /
        additional line
        """
    )
    params, others2 = read_namelist(io, "input"; all_lines = true)
    @test params == OrderedDict("a" => 1, "b" => 2.0, "c" => "test", "d" => true)
    @test others2 == ["additional line"]
end

@testitem "find_card" begin
    using QuantumEspressoIO: find_card
    lines = [
        "  ! This is a comment",
        "  input",
        "  name1 = value1  ! comment 2",
        "  name2 = value2",
    ]
    @test find_card(lines, "input") == 2
end

@testitem "parse_card_option" begin
    using QuantumEspressoIO: parse_card_option
    @test parse_card_option("POSITIONS angstrom  ! comment") == "angstrom"
    @test parse_card_option("POSITIONS {angstrom}  ! comment") == "angstrom"
end

@testitem "parse_card!" begin
    using QuantumEspressoIO: parse_card!
    lines = [
        "  input option1",
        "  ! This is a comment",
        "  name1 = value1  ! comment 2",
        "  name2 = value2",
    ]
    n = copy(lines)
    option, content = QuantumEspressoIO.parse_card!(n, "input", 1)
    @test option == "option1"
    @test content == ["name1 = value1"]
    @test n == ["  name2 = value2"]
end

@testitem "write_namelist" begin
    using OrderedCollections: OrderedDict
    name = "input"
    params_w = OrderedDict(
        "a" => 1,
        "b" => 2.0,
        "c" => [3, 4, 5],
        "d" => "test",
        "e" => true,
    )
    buf = IOBuffer()
    write_namelist(buf, name, params_w)
    out = String(take!(buf))
    expected = """&input
      a = 1
      b = 2.0
      c(1) = 3
      c(2) = 4
      c(3) = 5
      d = 'test'
      e = .true.
    /
    """
    @test out == expected
end

@testitem "write_namelists" begin
    using OrderedCollections: OrderedDict
    inputs = OrderedDict(
        "control" => OrderedDict("calculation" => "scf", "prefix" => "qe"),
        "system" => OrderedDict("ecutwfc" => 30.0, "ecutrho" => 300.0),
        "electrons" => OrderedDict("mixing_beta" => 0.7),
    )
    buf2 = IOBuffer()
    write_namelists(buf2, inputs)
    out2 = String(take!(buf2))
    expected = """&control
      calculation = 'scf'
      prefix = 'qe'
    /
    &system
      ecutwfc = 30.0
      ecutrho = 300.0
    /
    &electrons
      mixing_beta = 0.7
    /
    """
    @test out2 == expected
end
