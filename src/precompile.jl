using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    pw_input = """
    &system
      nat = 1
      ntyp = 1
    /
    ATOMIC_SPECIES
    Si 28.0855 Si.upf
    ATOMIC_POSITIONS crystal
    Si 0.0 0.0 0.0
    CELL_PARAMETERS angstrom
    1.0 0.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0
    K_POINTS crystal
    1
    0.0 0.0 0.0 1.0
    """

    band_output = """
    &plot nbnd=1, nks=1 /
    0.0 0.0 0.0
    -1.0
    """

    @compile_workload begin
        mktempdir() do path
            pwi = read_pw_in(IOBuffer(pw_input))
            tmp = joinpath(path, "pw.in")
            write_pw_in(tmp, pwi)
            read_pw_in(tmp)

            tmp = joinpath(path, "band.dat")
            write(tmp, band_output)
            read_band_dat(tmp)
        end
    end
end
