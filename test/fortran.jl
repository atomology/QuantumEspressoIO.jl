@testitem "format_fortran" begin
    using QuantumEspressoIO: format_fortran
    @test format_fortran("hello") == "'hello'"
    @test format_fortran(1.0) == "1.0"
    @test format_fortran(true) == ".true."
end

@testitem "parse_float" begin
    using QuantumEspressoIO: parse_float
    @test parse_float("1.0D-10") == 1.0e-10
    @test isnan(parse_float("1***"))
end

@testitem "parse_bool" begin
    using QuantumEspressoIO: parse_bool
    @test parse_bool(".true.") == true
    @test parse_bool("false") == false
    @test parse_bool("T") == true
    @test parse_bool("F") == false
    @test parse_bool("1") == true
    @test parse_bool("0") == false
    @test parse_bool(1) == true
    @test parse_bool(0) == false
end

@testitem "parse_value" begin
    using QuantumEspressoIO: parse_value
    @test parse_value("'hello'") == "hello"
    @test parse_value(".true.") == true
    @test parse_value("1.0D-10") == 1.0e-10
    @test parse_value("1") == 1
end

@testitem "remove_comment" begin
    using QuantumEspressoIO: remove_comment
    @test remove_comment("  ! This is a comment") == ""
    @test remove_comment("  input") == "input"
    @test remove_comment("  name1 = value1  ! comment 2") == "name1 = value1"
    @test remove_comment("  name2 = value2  # comment 3") == "name2 = value2"
end
