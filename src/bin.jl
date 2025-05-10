using FortranFiles

"""
Read wavefunction data from a QE's `wfc.dat` file.

# Arguments
- `filename::AbstractString`: The path to the `wfc.dat` file.

# Return
- `miller`: `3 * ngw`, integer matrix of Miller indices for reciprocal lattice vectors.
- `evc_list`: Length-`nbnd` vector, each element is a length-`igwx` vector of
    complex wavefunction coefficients.
"""
function read_wfc_dat(filename::AbstractString)
    f = FortranFile(filename)
    ik, xkx, xky, xkz, ispin = read(f, (Int32,5))
    ngw, igwx, npol, nbnd = read(f, (Int32,4))
    dummy_vector = read(f, (Float64,9))
    miller = reshape(read(f, (Int32,3*igwx)),(3, igwx))

    evc_list = []
    for _ in 1:nbnd
        evc = read(f, (ComplexF64,igwx))
        push!(evc_list,evc)
    end
    return miller, evc_list
end
