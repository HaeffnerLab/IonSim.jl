using IonSim: _alaguerre
export two_ion_ms, rabi_flop


"""
    two_ion_ms(tspan, Ω::Real, ν::Real, δ::Real, η::Real, n̄::Real)
[ref](https://doi.org/10.1103/PhysRevA.62.022311) <br>
Assumes vibrational mode starts in a thermal state with: ``\\langle a^\\dagger a\\rangle = n̄`` 
and ions start in doubly ground state. Returns `(ρgg, ρee)`, the population in the doubly 
ground and doubly excited state, respectively. ``[Ω], [ν], [δ] = Hz``
"""
function two_ion_ms(tspan, Ω::Real, ν::Real, δ::Real, η::Real, n̄::Real)
    ρgg = Float64[]; ρee = Float64[] 
    n̄ *= 1.0; Ω *= 2π; ν *= 2π; δ *= 2π
    for t in tspan
        F = -√2 * η * Ω * sin((ν - δ) * t) / (ν - δ)
        G =  -√2 * η * Ω * (1 - cos((ν - δ) * t)) / (ν - δ)
        A = -((η * Ω)^2 / (ν - δ)) * (t - (1/(2 * (ν - δ))) * sin(2 * (ν - δ) * t))
        term1 = (1/2) * exp(-(F^2 + G^2)/4)
        term2 = (1/8) * exp(-(F^2 + G^2))
        ρee_t = 0.0
        ρgg_t  = 0.0
        for n in 0:200
            pn = _Pn(n̄, n)
            isnan(pn) && continue
            t1_1 = _laguerre((F^2 + G^2) / 2, n) * cos(A + (1/2) * F * G) * term1
            t2_2 = _laguerre(2 * (F^2 + G^2), n) * term2
            ρgg_t += pn * (t1_1 + t2_2)
            ρee_t += pn * (-t1_1 + t2_2)
        end
        push!(ρgg, ρgg_t)
        push!(ρee, ρee_t)
    end
    3/8 .+ ρgg, 3/8 .+ ρee
end

"""
    rabi_flop(tspan, Ω::Real, η::Real, n̄::Real; s::Int=0) <br>
Single ion rabi flop. Returns:
``\\sum_{n=0}^∞ p_n sin^2(\\Omega_n t)`` <br> with
``\\Omega_n = Ωe^{-η^2/2}η^s\\sqrt{\\frac{n!}{(n+s)!}}L_{n}^{s}(η^2)`` <br>
where ``s`` is the order of the (blue) sideband that we are driving and ``L_{n}^{s}`` is the
associated Laguerre polynomial. [ref](https://doi.org/10.1103/RevModPhys.75.281)

"""
function rabi_flop(tspan, Ω::Real, η::Real, n̄::Real; s::Int=0)
    n̄ *= 1.0; Ω *= 2π
    p = Vector{Float64}(undef, 0)
    for t in tspan
        pi = 0.0
        for n in 0:300
            pi += _Pn(n̄, n) * sin(Ω/2 * exp(-η^2/2) * η^s * sqrt(1/prod(n+1:n+s)) * _alaguerre(η^2, n, s) * t)^2
        end
        push!(p, pi)
    end
    p
end

function _laguerre(x, n)
    L = 1.0, -x + 1
    if n < 2
        return L[n+1]
    end
    for i in 2:n
        L = L[2], ((2i - 1 - x) * L[2] - (i - 1) * L[1]) / i
    end
    L[2]
end

_Pn(n̄::Real, n::Int) = (n̄ / (n̄ + 1))^n / (n̄ + 1)

