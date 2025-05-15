using LinearAlgebra

"""
Convert crystall coordinates to Cartesian based on lattice matrix.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: crystal2cart)
julia> lattice = [[0.0, 2.715265, 2.715265], [2.715265, 0.0, 2.715265], [2.715265, 2.715265, 0.0]]
julia> crystal2cart(lattice, [0.25, 0.25, 0.25])
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 1.3576325
 1.3576325
 1.3576325
```
"""
function crystal2cart(lattice::AbstractVector{<:AbstractVector}, crystal_coord::AbstractVector{<:Real})
    cart_coord =  Vec3(map(row -> dot(row, crystal_coord), lattice)) # in Angstrom if was read from read_pw_xml 
    return  cart_coord
end

"""
Convert crystall coordinates of vector of atoms to Cartesian based on lattice matrix.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: crystal2cart)
julia> lattice = [[0.0, 2.715265, 2.715265], [2.715265, 0.0, 2.715265], [2.715265, 2.715265, 0.0]]
julia> positions = [[0,0,0],[0.25,0.25,0.25]]
julia> crystal2cart(lattice, positions)
2-element Vector{SVector{3, Float64}}:
 [0.0, 0.0, 0.0]
 [1.3576325, 1.3576325, 1.3576325]
```
"""
function crystal2cart(lattice::AbstractVector{<:AbstractVector}, crystal_coords::AbstractVector{<:AbstractVector})
    cart_coords = [crystal2cart(lattice, crystal_coord) for crystal_coord in crystal_coords]
    return  cart_coords
end

"""
Convert Cartesian coordinates to crystall based on lattice matrix.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: cart2crystall)
julia> lattice = [[0.0, 2.715265, 2.715265], [2.715265, 0.0, 2.715265], [2.715265, 2.715265, 0.0]]
julia> cart2crystall(lattice, [1.3576325, 1.3576325, 1.3576325])
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 0.25
 0.25
 0.25
```
"""
function cart2crystall(lattice::AbstractVector{<:AbstractVector}, cart_coord::AbstractVector{<:Real})
    lattice_mat = hcat(lattice...)'
    inv_lattice = inv(lattice_mat)
    return  Vec3(inv_lattice * cart_coord)
end

"""
Convert Cartesian coordinates of vector of atoms to crystall based on lattice matrix.

# Examples
```jldoctest; setup = :(using QuantumEspressoIO: cart2crystall)
julia> lattice = [[0.0, 2.715265, 2.715265], [2.715265, 0.0, 2.715265], [2.715265, 2.715265, 0.0]]
julia> positions = [[0,0,0],[1.3576325, 1.3576325, 1.3576325]]
julia> crystal2cart(lattice, positions)
2-element Vector{SVector{3, Float64}}:
 [0.0, 0.0, 0.0]
 [0.25, 0.25, 0.25]
```
"""
function cart2crystall(lattice::AbstractVector{<:AbstractVector}, cart_coords::AbstractVector{<:AbstractVector})
    crystal_coords =  [cart2crystall(lattice, cart_coord) for cart_coord in cart_coords]
    return  crystal_coords
end
