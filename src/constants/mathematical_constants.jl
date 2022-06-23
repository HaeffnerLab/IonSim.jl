export eye3, c_rank1, c_rank2, x̂, ŷ, ẑ, ndot

"""`eye3`: Three-dimensional identity matrix"""
const eye3 = [1 0 0; 0 1 0; 0 0 1]

"""`c_rank1`: Matrix of spherical basis vectors (defined in e.g. Quantum dynamics of cold trapped ions with application to quantum computation, Appl. Phys. B 66, 181-190 (1998).

Useful for converting between coordinates of rank-1 spherical tensors and Cartesian coordinates"""
const c_rank1 = [
    1 / sqrt(2) * [1 im 0] #q=-1
    [0 0 1] #q=0
    -1 / sqrt(2) * [1 -im 0]
] #q=1

"""`c_rank2`: Matrix of spherical basis rank-2 tensors (defined in e.g. Quantum dynamics of cold trapped ions with application to quantum computation, Appl. Phys. B 66, 181-190 (1998).

Useful for converting between coordinates of rank-2 spherical tensors and Cartesian coordinates"""
const c_rank2 = cat(
    1 / sqrt(6) * [[1, im, 0] [im, -1, 0] [0, 0, 0]], #q=-2
    1 / sqrt(6) * [[0, 0, 1] [0, 0, im] [1, im, 0]], #q=-1
    1 / 3 * [[-1, 0, 0] [0, -1, 0] [0, 0, 2]], #q=0
    1 / sqrt(6) * [[0, 0, -1] [0, 0, im] [-1, im, 0]], #q=1
    1 / sqrt(6) * [[1, -im, 0] [-im, -1, 0] [0, 0, 0]];
    dims = 3
) #q=2


#############################################################################################
# Cartesian unit vectors
#############################################################################################

const x̂ = (x = 1, y = 0, z = 0)
const ŷ = (x = 0, y = 1, z = 0)
const ẑ = (x = 0, y = 0, z = 1)

function _print_axis(a::NamedTuple{(:x, :y, :z)})
    if a == x̂
        return "x̂"
    elseif a == ŷ
        return "ŷ"
    elseif a == ẑ
        return "ẑ"
    else
        return string(a)
    end
end

ndot(a::NamedTuple{(:x, :y, :z)}, b::NamedTuple{(:x, :y, :z)}) =
    a.x * b.x + a.y * b.y + a.z * b.z
function Base.:+(a::NamedTuple{(:x, :y, :z)}, b::NamedTuple{(:x, :y, :z)})
    return (x = a.x + b.x, y = a.y + b.y, z = a.z + b.z)
end
function Base.:-(a::NamedTuple{(:x, :y, :z)}, b::NamedTuple{(:x, :y, :z)})
    return (x = a.x - b.x, y = a.y - b.y, z = a.z - b.z)
end
Base.:/(a::NamedTuple{(:x, :y, :z)}, b::Number) = (x = a.x / b, y = a.y / b, z = a.z / b)
Base.:*(a::NamedTuple{(:x, :y, :z)}, b::Number) = (x = a.x * b, y = a.y * b, z = a.z * b)
Base.:*(b::Number, a::NamedTuple{(:x, :y, :z)}) = (x = a.x * b, y = a.y * b, z = a.z * b)
