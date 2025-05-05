@testitem "write_namelist" begin
    tmpfile = tempname(; cleanup=true)

    name = "test"
    namelist = Dict(
        "a" => 1,
    )
    write_namelist(tmpfile, name, namelist)
    ref = """
    &test
      a = 1
    /
    """
    @test readlines(tmpfile) == split(ref)
end
