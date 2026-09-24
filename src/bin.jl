using FortranFiles

"""
Read wavefunction data from a QE's `wfc.dat` file.

# Arguments
- `filename::AbstractString`: The path to the `wfc.dat` file.

# Return
- `miller`: `3 * ngw`, integer matrix of Miller indices for reciprocal lattice vectors.
- `evc`: `igwx × nbnd` matrix of complex wavefunction coefficients, one column
    per band.
"""
function read_wfc_dat(filename::AbstractString)
    f = FortranFile(filename)
    ik, xkx, xky, xkz, ispin = read(f, (Int32,5))
    ngw, igwx, npol, nbnd = read(f, (Int32,4))
    dummy_vector = read(f, (Float64,9))
    miller = reshape(read(f, (Int32,3*igwx)),(3, igwx))

    evc = zeros(ComplexF64, igwx, nbnd)
    for ib in 1:nbnd
        evc[:, ib] = read(f, (ComplexF64,igwx))
    end
    return miller, evc
end
