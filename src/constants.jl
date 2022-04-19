module PhysicalConstants

import Base.sqrt

export μB, ħ, c, e, ϵ₀, α, kB, eye3, c_rank1, c_rank2

"""
    PhysicalConstant(x::Real)
Useful physical constants. Values are in SI units.
"""
struct PhysicalConstant <: Real
    x::Real
    units::String
end

#############################################################################################
# Physical constants (everything in SI units)
#############################################################################################

"""`μB` = 9.27400994e-24 J⋅T⁻¹ <br> (Bohr Magneton)"""
const μB = PhysicalConstant(9.27400994e-24, "J⋅T⁻¹")
"""`ħ` = 1.0545718e-34 m²kg/s <br> (Planck's constant / 2π)"""
const ħ = PhysicalConstant(1.0545718e-34, "m²kg/s")
"""`c` = 2.99792458e8 m/s <br> (speed of light in vacuum)"""
const c = PhysicalConstant(2.99792458e8, "m/s")
"""`e` = 1.60217662e-19 C <br> (charge of electron)"""
const e = PhysicalConstant(1.60217662e-19, "C")
"""`ϵ₀` = 8.85418782e-12 ``(s^4A^2) / (m^3 kg)``"""
const ϵ₀ = PhysicalConstant(8.85418782e-12, "(s^4A^2) / (m^3 kg)")
"""`α` = e²/4πϵ₀ħc``"""
const α = PhysicalConstant(0.007297352557920479, "")
"""`kB` = 1.38064852e-23 ``m^2kg/(s^2K)``"""
const kB = PhysicalConstant(1.38064852e-23, "m^2kg/(s^2K)")


#############################################################################################
# 3D real-space tensors
#############################################################################################

"""`eye3`: Three-dimensional identity matrix"""
const eye3 = [1 0 0; 0 1 0; 0 0 1]

"""`c_rank1`: Matrix of spherical basis vectors (defined in e.g. Quantum dynamics of cold trapped ions with application to quantum computation, Appl. Phys. B 66, 181-190 (1998).

Useful for converting between coordinates of rank-1 spherical tensors and Cartesian coordinates"""
const c_rank1 = [1  /  sqrt(2)  *  [1 im 0]; #q=-1
                [0 0 1]; #q=0
                -1  /  sqrt(2)  *  [1 -im 0]] #q=1

"""`c_rank2`: Matrix of spherical basis rank-2 tensors (defined in e.g. Quantum dynamics of cold trapped ions with application to quantum computation, Appl. Phys. B 66, 181-190 (1998).

Useful for converting between coordinates of rank-2 spherical tensors and Cartesian coordinates"""
                const c_rank2 = cat(1/sqrt(6)  *  [[1, im, 0] [im, -1, 0] [0, 0, 0]], #q=-2
                1/sqrt(6)  *  [[0, 0, 1] [0, 0, im] [1, im, 0]], #q=-1
                    1/3  *  [[-1, 0, 0] [0, -1, 0] [0, 0, 2]], #q=0
                1/sqrt(6)  *  [[0, 0, -1] [0, 0, im] [-1, im, 0]], #q=1
                1/sqrt(6)  *  [[1, -im, 0] [-im, -1, 0] [0, 0, 0]]; dims = 3) #q=2



Base.print(pc::PhysicalConstant) = print("$(pc.x) [$(pc.units)]")
Base.show(io::IO, pc::PhysicalConstant) = print(io, "$(pc.x) [$(pc.units)]")

Base.convert(::Type{<:Number}, x::PhysicalConstant) = map(x->x.x, x)
Base.promote_rule(::Number, ::PhysicalConstant) = Number
Base.Fix2(f, x::PhysicalConstant) = Base.Fix2(f, x.x)

Base.:*(x1::PhysicalConstant, x2::PhysicalConstant) = x1.x * x2.x
Base.:/(x1::PhysicalConstant, x2::PhysicalConstant) = x1.x / x2.x
Base.:+(x1::PhysicalConstant, x2::PhysicalConstant) = x1.x + x2.x
Base.:-(x1::PhysicalConstant, x2::PhysicalConstant) = x1.x - x2.x
Base.:^(x1::PhysicalConstant, x2::PhysicalConstant) = x1.x^x2.x
Base.:(>)(x1::PhysicalConstant, x2::PhysicalConstant) = Base.:(>)(x1.x, x2.x)
Base.:(<)(x1::PhysicalConstant, x2::PhysicalConstant) = Base.:(<)(x1.x, x2.x)
Base.:(≥)(x1::PhysicalConstant, x2::PhysicalConstant) = Base.:(≥)(x1.x, x2.x)
Base.:(≤)(x1::PhysicalConstant, x2::PhysicalConstant) = Base.:(≤)(x1.x, x2.x)

sqrt(x1::PhysicalConstant) = sqrt(x1.x)

end  # module


export x̂, ŷ, ẑ, ndot

const x̂ = (x=1, y=0, z=0)
const ŷ = (x=0, y=1, z=0)
const ẑ = (x=0, y=0, z=1)

function _print_axis(a::NamedTuple{(:x,:y,:z)})
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

ndot(a::NamedTuple{(:x,:y,:z)}, b::NamedTuple{(:x,:y,:z)}) = a.x * b.x + a.y * b.y + a.z * b.z
function Base.:+(a::NamedTuple{(:x,:y,:z)}, b::NamedTuple{(:x,:y,:z)})
    (x=a.x + b.x, y=a.y + b.y, z=a.z + b.z)
end
function Base.:-(a::NamedTuple{(:x,:y,:z)}, b::NamedTuple{(:x,:y,:z)})
    (x=a.x - b.x, y=a.y - b.y, z=a.z - b.z)
end
Base.:/(a::NamedTuple{(:x,:y,:z)}, b::Number) = (x=a.x/b, y=a.y/b, z=a.z/b)
Base.:*(a::NamedTuple{(:x,:y,:z)}, b::Number) = (x=a.x*b, y=a.y*b, z=a.z*b)
Base.:*(b::Number, a::NamedTuple{(:x,:y,:z)}) = (x=a.x*b, y=a.y*b, z=a.z*b)