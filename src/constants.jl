module PhysicalConstants

import Base.sqrt

export μB, ħ, m_ca40, c, e, α, ϵ₀, kB, ca40_qubit_transition_frequency

"""
    PhysicalConstant(x::Real)
Useful physical constants. Values are in SI units.
"""
struct PhysicalConstant <: Real
    x::Real
    units::String
end

# some useful constants, everything in SI units
"""`μB` = 9.27400994e-24 J⋅T⁻¹ <br> (Bohr Magneton)"""
const μB = PhysicalConstant(9.27400994e-24, "J⋅T⁻¹")
"""`ħ` = 1.0545718e-34 m²kg/s <br> (Planck's constant / 2π)"""
const ħ = PhysicalConstant(1.0545718e-34, "m²kg/s")
"""`m_ca40` = 6.635943757345042e-26 kg <br> (mass of 40Ca)"""
const m_ca40 = PhysicalConstant(6.635943757345042e-26, "kg")
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
"""`ca40_qubit_transition_frequency` = c / 729.147e-9 ``Hz`` """
const ca40_qubit_transition_frequency = PhysicalConstant(2.99792458e8 / 729.147e-9, "Hz")

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
