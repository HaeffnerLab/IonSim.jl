import QuantumOptics.timeevolution
import QuantumOptics.stochastic
using QuantumOptics: Basis, Operator, CompositeBasis, Ket, StateVector
using QuantumOptics.timeevolution: DecayRates


#############################################################################################
# Wrap QuantumOptics solvers in order to implement our form of basis check
#############################################################################################

# function timeevolution.master_dynamic(tspan, rho0::T, f::Function;
#             rates::DecayRates=nothing, fout::Union{Function,Nothing}=nothing, kwargs...
#         ) where {B<:Basis,T<:Operator{B,B}}
#     check_bases(f(0., 0.).basis_l, rho0.basis_l)
#     timeevolution.master_dynamic(tspan, rho0, f; rates=rates, fout=fout, kwargs...)
# end

# function timeevolution.mcwf_dynamic(tspan, psi0::T, f::Function;
#     seed=rand(UInt), rates::DecayRates=nothing,
#     fout=nothing, display_beforeevent=false, display_afterevent=false,
#     kwargs...) where {T<:Ket}
#     check_bases(f(0., 0.).basis_l, psi0.basis)
#     timeevolution.mcwf_dynamic(tspan, psi0, f; seed=seed, rates=rates, fout=fout, 
#             display_beforeevent=display_beforeevent, display_afterevent=display_afterevent,
#             kwargs...
#         )
# end

# function stochastic.schroedinger_dynamic(tspan, psi0::T, fdeterm::Function, fstoch::Function;
#             fout::Union{Function,Nothing}=nothing, noise_processes::Int=0,
#             normalize_state::Bool=false, kwargs...
#         ) where T<:Ket
#     check_bases(fdeterm(0.0, 0.0), fstoch)
#     stochastic.schroedinger_dynamic(tspan, psi0, fdeterm, fstoch;
#             fout=fout, noise_processes=noise_processes, normalize_state=normalize_state, 
#             kwargs...
#         )
# end

# function stochastic.master_dynamic(tspan, rho0::T, fdeterm::Function, fstoch::Function;
#             rates::DecayRates=nothing, fout::Union{Function,Nothing}=nothing,
#             noise_processes::Int=0, kwargs...
#         ) where {B<:Basis,T<:Operator{B,B}}
#         check_bases()
#         stochastic.master_dynamic(tspan, rho0, fdeterm, fstoch; rates=rates, fout=fout,
#                 noise_processes=noise_processes, kwargs...
#             )
# end

# function timeevolution.schroedinger_dynamic(tspan, psi0::T, f::Function;
#             fout::Union{Function,Nothing}=nothing, kwargs...
#         ) where T<:StateVector
#     check_bases(f(0., 0.).basis_l, psi0.basis)
#     timeevolution.schroedinger_dynamic(tspan, psi0, f; fout=fout, kwargs...)
# end

# when solving for unitary evolution, it is convenient to switch back and forth between an 
# initial state that is pure and one that is mixed
function timeevolution.schroedinger_dynamic(tspan, rho0::T, f::Function;
            fout::Union{Function,Nothing}=nothing, kwargs...
        ) where {B<:Basis,T<:Operator{B,B}}
    check_bases(f(0., 0.).basis_l, rho0.basis_l)
    Jvec = []
    timeevolution.master_dynamic(tspan, rho0, (t,rho) -> (f(t, rho), Jvec, Jvec);
        fout=fout, kwargs...)
end

check_bases(b1::CompositeBasis, b2::CompositeBasis) = b1 == b2 || throw(IncompatibleBases())

mutable struct IncompatibleBases <: Exception end