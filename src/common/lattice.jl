using LinearAlgebra

"""
    $(SIGNATURES)

Compute reciprocal-space lattice vectors from real-space lattice vectors.

# Arguments
- `lattice`: Can be a vector of lattice vectors, or a [`Mat3`](@ref) matrix.

# Examples
```jldoctest get_recip_lattice; setup = :(using QuantumEspressoIO: get_recip_lattice, mat3)
lattice = [[0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]];
get_recip_lattice(lattice)
# output
3-element Vector{Vector{Float64}}:
 [-2.6927937030769655, 2.243994752564138, 2.019595277307724]
 [1.3463968515384828, -1.121997376282069, 0.5609986881410345]
 [0.4487989505128276, 0.6731984257692414, -0.3365992128846207]
```
```jldoctest get_recip_lattice
lattice = mat3(lattice);
get_recip_lattice(lattice)
# output
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -2.69279   1.3464     0.448799
  2.24399  -1.122      0.673198
  2.0196    0.560999  -0.336599
```
"""
function get_recip_lattice end

function get_recip_lattice(lattice::Mat3)
    # In case of Mat3, return a Mat3 as well
    return 2π * inv(lattice)'
end

function get_recip_lattice(lattice::AbstractVector{<:AbstractVector})
    M = reduce(hcat, lattice)
    recip = 2π * inv(M)'
    # Return a vector of vectors
    return collect.(eachcol(recip))
end

"""
    $(SIGNATURES)

Compute real-space lattice vectors from reciprocal-space lattice vectors.

# Arguments
- `recip_lattice`: Can be a vector of reciprocal lattice vectors, or a [`Mat3`](@ref) matrix.

# Examples
```jldoctest get_lattice; setup = :(using QuantumEspressoIO: get_lattice, mat3, get_recip_lattice)
recip_lattice = [[0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]];
get_lattice(recip_lattice)
# output
3-element Vector{Vector{Float64}}:
 [-2.6927937030769655, 2.243994752564138, 2.019595277307724]
 [1.3463968515384828, -1.121997376282069, 0.5609986881410345]
 [0.4487989505128276, 0.6731984257692414, -0.3365992128846207]
```
```jldoctest get_lattice
recip_lattice = mat3(recip_lattice);
get_lattice(recip_lattice)
# output
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -2.69279   1.3464     0.448799
  2.24399  -1.122      0.673198
  2.0196    0.560999  -0.336599
```
```jldoctest get_lattice
get_recip_lattice(get_lattice(recip_lattice)) ≈ recip_lattice
# output
true
```
"""
const get_lattice = get_recip_lattice

"""
    $(SIGNATURES)

Convert fractional to Cartesian coordinates based on lattice vectors.

# Arguments
- `lattice`: Each element is a lattice vector.
    - For lattice vectors, the unit is usually in angstrom.
    - For reciprocal lattice vectors, the unit is usually in 1/angstrom.
- `vec`: a vector or a list of vectors in fractional coordinates.

# Examples
```jldoctest frac2cart; setup = :(using QuantumEspressoIO: frac2cart)
lattice = [[0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]];
positions = [[0.1, 0.2, 0.3], [1.0, 2.0, 3.0]];
frac2cart(lattice, positions[1])
# output
3-element Vector{Float64}:
 2.1
 1.9
 1.0
```
```jldoctest frac2cart
frac2cart(lattice, positions)
# output
2-element Vector{Vector{Float64}}:
 [2.1, 1.9, 1.0]
 [21.0, 19.0, 10.0]
```
"""
function frac2cart end

function frac2cart(lattice::AbstractVector{<:AbstractVector}, vecs::AbstractVector{<:AbstractVector})
    mat = reduce(hcat, lattice)
    carts = Ref(mat) .* vecs
    return carts
end

function frac2cart(lattice::AbstractVector{<:AbstractVector}, vec::AbstractVector{<:Real})
    return frac2cart(lattice, [vec])[1]
end

"""
    $(SIGNATURES)

Convert Cartesian to fractional coordinates based on lattice vectors.

# Arguments
- `lattice`: Each element is a lattice vector.
    - For lattice vectors, the unit is usually in angstrom.
    - For reciprocal lattice vectors, the unit is usually in 1/angstrom.
- `vec`: a vector or a list of vectors in Cartesian coordinates.

# Examples
```jldoctest cart2frac; setup = :(using QuantumEspressoIO: cart2frac, frac2cart)
lattice = [[0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]];
positions = [[2.1, 1.9, 1.0], [21.0, 19.0, 10.0]];
frac2cart(lattice, cart2frac(lattice, positions[1])) ≈ positions[1]
# output
true
```
```jldoctest cart2frac
frac2cart(lattice, cart2frac(lattice, positions)) ≈ positions
# output
true
```
"""
function cart2frac end

function cart2frac(lattice::AbstractVector{<:AbstractVector}, vecs::AbstractVector{<:AbstractVector})
    mat = reduce(hcat, lattice)
    inv_mat = inv(mat)
    frac = Ref(inv_mat) .* vecs
    return frac
end

function cart2frac(lattice::AbstractVector{<:AbstractVector}, vec::AbstractVector{<:Real})
    return cart2frac(lattice, [vec])[1]
end
