module PhysicalConstants

import Base.sqrt
export μB, ħ, c, e, ϵ₀, α, kB

"""
    PhysicalConstant(x::Real)
Useful physical constants. Values are in SI units.
"""
struct PhysicalConstant <: Real
    x::Real
    units::String
end

Base.print(pc::PhysicalConstant) = print("$(pc.x) [$(pc.units)]")
Base.show(io::IO, pc::PhysicalConstant) = print(io, "$(pc.x) [$(pc.units)]")

Base.convert(::Type{<:Number}, x::PhysicalConstant) = map(x -> x.x, x)
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

end  # module
