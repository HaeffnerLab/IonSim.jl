export Bfield, set_gradient!

"""
    Bfield(T::Trap, ion::Ion)
Retuns the value of the magnetic field in `T` at the location of `ion`, including both the trap's overall B-field and its B-field gradient.
"""
function Bfield(T::Trap, ion::Ion)
    @assert ionintrap(T, ion) "trap does not contain ion"
    return T.B + T.∇B * ionposition(ion)
end

#TODO: change this function name to set_B_field_gradient!, or something similarily clear, and update unit testing
"""
    set_gradient!(
            T::Trap, ion_indxs::Tuple{Int,Int}, transition::Tuple, f::Real
        )
Sets the Bfield gradient in place to achieve a detuning `f` between the `transition` of two
ions, which are assumed to be of the same species. `ion_indxs` refer to the
ordering of the ions in the chain.
"""
function set_gradient!(T::Trap, ion_indxs::Tuple{Int, Int}, transition::Tuple, f::Real)
    ionA = T.configuration.ions[ion_indxs[1]]
    ionB = T.configuration.ions[ion_indxs[2]]
    separation = abs(ionposition(ionA) - ionposition(ionB))

    (SL1, SL2) = transition
    L1 = sublevel2level(ionA, SL1)
    L2 = sublevel2level(ionA, SL2)
    g1 = landegf(ionA, L1)
    g2 = landegf(ionA, L2)
    m1 = quantumnumbers(ionA, SL1).m
    m2 = quantumnumbers(ionA, SL2).m
    # Calculate Zeeman shifts with a unit B-field using a method of zeeman_shift that ensures a nonlinear term is not used
    E1 = zeeman_shift(1.0, g1, m1)
    E2 = zeeman_shift(1.0, g2, m2)
    return T.∇B = f / (abs(E2 - E1) * separation)
end
