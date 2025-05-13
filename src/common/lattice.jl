using LinearAlgebra

# Function to transform a vector from crystal to Cartesian coordinates
function crystal2cart(lattice::Vector{SVector{3, Float64}}, crystal_coord::SVector{3, Float64})
    cart_coord =  SA_F64[map(row -> dot(row, crystal_coord), lattice)...] # in Angstrom if was read from read_qe_xml 
    return  cart_coord
end

function crystal2cart(lattice::Vector{SVector{3, Float64}}, crystal_coords::Vector{SVector{3, Float64}})
    cart_coords = [crystal2cart(lattice, crystal_coord) for crystal_coord in crystal_coords]
    return  cart_coords
end

# Function to transform a vector from Cartesian (Angstrom) to crystal coordinates
function cart2crystall(lattice::Vector{SVector{3, Float64}}, cart_coord::SVector{3, Float64})
    lattice_mat = hcat(lattice...)'
    inv_lattice = inv(lattice_mat)
    return  inv_lattice * cart_coord
end

function cart2crystall(lattice::Vector{SVector{3, Float64}}, cart_coords::Vector{SVector{3, Float64}})
    crystal_coords =  [cart2crystall(lattice, cart_coord) for cart_coord in cart_coords]
    return  crystal_coords
end
