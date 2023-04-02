"""
    MagneticField
Contains information about the static magnetic field (except for fluctuations) in the
vicinity of the ions in a Chamber.
**fields**
* `B::Real` Mean magnitude of the B-field [Tesla]
* `Bhat::NamedTuple{(:x,:y,:z)}`: Direction of the B-field
* `δB::Function`: Time-dependence of the magnitude B-field [Tesla]
"""
abstract type MagneticField end


"""
    LinearMagneticField(B::Real;
                δB::Union{Real,Function}=0,
                Bhat::NamedTuple{(:x,:y,:z)=ẑ,
                ∇B::Real=0
    )
MageticField for linear ion chains. Features a static gradient of the magnitude in the
z-direction only.
**user-defined fields**
* `B`: A real value describing the mean magnitude of the B-field [Tesla].
* `Bhat::NamedTuple{(:x,:y,:z)}`: Describes the direction of the B-field (defaults to ẑ).
* `δB::Function`: Time-dependence of the B-field [Tesla]
* `∇B`: Magnitude of the B-field gradient. We assume that the gradient always points along
        the z-direction. [Tesla / meter]
**derived fields**
* `_cnst_δB::Bool`: A Boolean flag signifying whether or not `δB` is a constant function.
"""
mutable struct LinearMagneticField <: MagneticField
    B::Real
    δB::Function
    Bhat::NamedTuple{(:x, :y, :z)}
    ∇B::Real
    _cnst_δB::Bool
    function MageticField(B;
        Bhat=ẑ,
        ∇B=0,
        δB::TδB=0,
    ) where {TδB}
        @assert isapprox(norm(Bhat), 1, rtol=1e-6) "!(|$Bhat| = 1)"
        if TδB <: Number
            _cnst_δB = true
            δBt(t) = δB
        else
            _cnst_δB = false
            δBt = δB
        end
        return new(B, Bhat, ∇B, δBt, _cnst_δB)
    end
end

magnitude(magfield::MagneticField) = magfield.B
magnitude_fluctuation(magfield::MagneticField) = magfield.δB
direction(magfield::MagneticField) = magfield.Bhat
gradient(magfield::LinearMagneticField) = magfield.∇B

function magnitude!(magfield::MagneticField, B::Real)
    magfield.B = B
end

function direction!(magfield::MagneticField, Bhat::NamedTuple{(:x, :y, :z)})
    @assert isapprox(norm(Bhat), 1, rtol=1e-6) "!(|$Bhat| = 1)"
    magfield.Bhat = Bhat
end

function magnitude_fluctuation!(magfield::MagneticField, δB::Function)
    magfield.δB = δB
    magfield._cnst_δB = false
end

function magnitude_fluctuation!(magfield::Chamber, δB::Real)
    magfield.δB = (t -> δB)
    magfield._cnst_δB = true
end

function gradient!(magfield::LinearMagneticField, ∇B::Real)
    magfield.∇B = ∇B
end