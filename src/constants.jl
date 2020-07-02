module PhysicalConstants

import Base.sqrt

export μB, ħ, m_ca40, m_ca43, m_be9, m_yb171, m_ba138, m_sr88, m_mg25, m_hg198, m_hg199, c, e, α, ϵ₀, kB, ca40_qubit_transition_frequency

"""
    PhysicalConstant(x::Real)
Useful physical constants. Values are in SI units.
"""
struct PhysicalConstant <: Real
    x::Real
    units::String
end

# some useful constants, everything in SI units
""" ## `m_ca40` = 6.6359443331e-26 kg <br> (mass of 40Ca)"""
const m_ca40 = PhysicalConstant(6.6359443331e-26, "kg")
""" ## `m_ca43` = 7.133470993e-26 <br> (mass of 43Ca)"""
const m_ca43 = PhysicalConstant(7.133470993e-26, "kg")
""" ## `m_be9` = 1.496508205e-26 <br> (mass of 9Be)"""
const m_be9 = PhysicalConstant(1.496508205e-26, "kg")
""" ## `m_yb171` = 2.838464542e-25 <br> (mass of 171Yb)"""
const m_yb171 = PhysicalConstant(2.838464542e-25, "kg")
""" ## `m_ba138` = 2.2899705013e-25 <br> (mass of 138Ba)"""
const m_ba138 = PhysicalConstant(2.2899705013e-25, "kg")
""" ## `m_sr88` = 1.459707037e-25 <br> (mass of 88Sr)"""
const m_sr88 = PhysicalConstant(1.459707037e-25, "kg")
""" ## `m_mg25` = 4.1489958410e-26 <br> (mass of 25Mg)"""
const m_mg25 = PhysicalConstant(4.1489958410e-26, "kg")
""" ## `m_hg198` = 3.2873155315e-25 <br> (mass of 198Hg)"""
const m_hg198 = PhysicalConstant(3.2873155315e-25, "kg")
""" ## `m_hg199` = 3.3039460302e-25 <br> (mass of 199Hg)"""
const m_hg199 = PhysicalConstant(3.3039460302e-25, "kg")

""" ## `μB` = 9.27400994e-24 J⋅T⁻¹ <br> (Bohr Magneton)"""
const μB = PhysicalConstant(9.27400994e-24, "J⋅T⁻¹")
""" ## `ħ` = 1.0545718e-34 m²kg/s <br> (Planck's constant / 2π)"""
const ħ = PhysicalConstant(1.0545718e-34, "m²kg/s")
""" ## `c` = 2.99792458e8 m/s <br> (speed of light in vacuum)"""
const c = PhysicalConstant(2.99792458e8, "m/s")
""" ## `e` = 1.60217662e-19 C <br> (charge of electron)"""
const e = PhysicalConstant(1.60217662e-19, "C")
""" ## `ϵ₀` = 8.85418782e-12 ``(s^4A^2) / (m^3 kg)``"""
const ϵ₀ = PhysicalConstant(8.85418782e-12, "(s^4A^2) / (m^3 kg)")
""" ## `α` = e²/4πϵ₀ħc``"""
const α = PhysicalConstant(0.007297352557920479, "")
""" ## `kB` = 1.38064852e-23 ``m^2kg/(s^2K)``"""
const kB = PhysicalConstant(1.38064852e-23, "m^2kg/(s^2K)")
""" ## `ca40_qubit_transition_frequency` = c / 729.147e-9 ``Hz`` """
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