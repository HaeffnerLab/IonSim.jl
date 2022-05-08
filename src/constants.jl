module PhysicalConstants
using Unitful

import Base.sqrt

export μB, ħ, c, e, ϵ₀, α, eye3, c_rank1, c_rank2, INVERSE_TIME, MAGNETIC

#############################################################################################
# Physical constants (everything in SI units)
#############################################################################################

"""`μB` = 9.27400994e-24 J⋅T⁻¹ <br> (Bohr Magneton)"""
const μB = u"μB"
"""`ħ` = 1.0545718e-34 m²kg/s <br> (Planck's constant / 2π)"""
const ħ = u"ħ"
"""`c` = 2.99792458e8 m/s <br> (speed of light in vacuum)"""
const c =  u"c0"
"""`e` = 1.60217662e-19 C <br> (charge of electron)"""
const e = u"q"
"""`ϵ₀` = 8.85418782e-12 ``(s^4A^2) / (m^3 kg)``"""
const ϵ₀ = u"ϵ0"
"""`α` = e²/4πϵ₀ħc``"""
const α = u"q"*u"q"/4/pi/u"ϵ0"/u"ħ"/u"c0" |> NoUnits

const INVERSE_TIME = Union{typeof(1u"Hz"), typeof(1.0u"Hz")}
const MAGNETIC = Union{typeof(1u"T"), typeof(1.0u"T")}

#############################################################################################
# 3D real-space tensors
#############################################################################################

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

end  # module

export x̂, ŷ, ẑ, ndot

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
