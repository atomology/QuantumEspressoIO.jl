using LinearAlgebra

"""
    $(SIGNATURES)

Convert fractional to Cartesian coordinates based on lattice vectors.

# Arguments
- `lattice`: Each element is a lattice vector.
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
