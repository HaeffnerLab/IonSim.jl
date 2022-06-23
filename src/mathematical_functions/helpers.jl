# internal functions

# computes iⁿ(-i)ᵐ * (s! / ((s+1) * √(m!n!)))
function _pf(s::Int, n::Int, m::Int)
    n -= 1
    m -= 1
    s -= 1
    @assert n <= s && m <= s
    val = 1.0 / (s + 1)
    for i in 0:(s - 2)
        if (m - i > 0) && (n - i > 0)
            val *= (s - i) / (√((m - i) * (n - i)))
        elseif m - i > 0
            val *= (s - i) / (√(m - i))
        elseif n - i > 0
            val *= (s - i) / (√(n - i))
        else
            val *= (s - i)
        end
    end
    return (-1im)^n * 1im^m * val
end

# computes the coefficients for the 'probabilist's' Hermite polynomial of order n
function _He(n::Int)
    a = zeros(Float64, n + 2, n + 2)
    a[1, 1] = 1
    a[2, 1] = 0
    a[2, 2] = 1
    for i in 2:(n + 1), j in 1:(n + 1)
        if j ≡ 1
            a[i + 1, j] = -(i - 1) * a[i - 1, j]
        else
            a[i + 1, j] = a[i, j - 1] - (i - 1) * a[i - 1, j]
        end
    end
    return [a[n + 1, k + 1] for k in 0:n]
end

# computes He_n(x) (nth order Hermite polynomial)
function _fHe(x::Real, n::Int)
    n -= 1
    He = 1.0, x
    if n < 2
        return He[n + 1]
    end
    for i in 2:n
        He = He[2], x * He[2] - (i - 1) * He[1]
    end
    return He[2]
end

# computes the matrix elements ⟨m|Dˢ(α)|n⟩ for the truncated displacement operator Dˢ(α)
# which exists in a Hilbert space of dimension s
function _Dtrunc(Ω, Δ, η, ν, rs, s, n, prefactor, timescale, L, t)
    d = complex(1, 0)
    for i in 1:L
        val = 0.0
        Δn = n[1][i] - n[2][i]
        for r in rs[i]
            val +=
                exp(im * r * abs(η[i])) * _fHe(r, n[2][i]) * _fHe(r, n[1][i]) /
                _fHe(r, s[i])^2
        end
        d *= (
            exp(im * Δn * (2π * ν[i] * timescale * t + π / 2 + π * (sign(η[i] < 0)))) *
            val *
            prefactor[i]
        )
    end
    g = Ω * exp(-1im * t * Δ)
    return g * d, g * conj(d)
end

function _laguerre(x, n)
    L = 1.0, -x + 1
    if n < 2
        return L[n + 1]
    end
    for i in 2:n
        L = L[2], ((2i - 1 - x) * L[2] - (i - 1) * L[1]) / i
    end
    return L[2]
end

_Pn(n̄::Real, n::Int) = (n̄ / (n̄ + 1))^n / (n̄ + 1)

# associated Laguerre polynomial
function _alaguerre(x::Real, n::Int, k::Int)
    L = 1.0, -x + k + 1
    if n < 2
        return L[n + 1]
    end
    for i in 2:n
        L = L[2], ((k + 2i - 1 - x) * L[2] - (k + i - 1) * L[1]) / i
    end
    return L[2]
end

# matrix elements of the displacement operator in the Fock Basis, assuming an
# infinite-dimensional Hilbert space. https://doi.org/10.1103/PhysRev.177.1857
function _Dnm(ξ::Number, n::Int, m::Int)
    if n < m
        return (-1)^isodd(abs(n - m)) * conj(_Dnm(ξ, m, n))
    end
    n -= 1
    m -= 1
    s = 1.0
    for i in (m + 1):n
        s *= i
    end
    ret = sqrt(1 / s) * ξ^(n - m) * exp(-abs2(ξ) / 2.0) * _alaguerre(abs2(ξ), m, n - m)
    if isnan(ret)
        return 1.0 * (n == m)
    end
    return ret
end
