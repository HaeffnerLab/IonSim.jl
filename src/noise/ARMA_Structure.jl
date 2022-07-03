
# == ARMA Power Spectrum Estimate == #
@inline function arma_psd(sigma, phi, theta, freq, T; phi0 = -1.0, theta0 = phi0)
    Zlist = @. exp(-1im * 2pi * freq * T)
    # phi0 = theta0 = -1 in notation
    Θ = map(z -> -theta0 - dot(theta, z .^ (1:length(theta))), Zlist)
    Φ = map(z -> -phi0 - dot(phi, z .^ (1:length(phi))), Zlist)
    return @. sigma^2 * abs(Θ / Φ)^2
end
ar_psd(sigma, phi, freq, T) = arma_psd(sigma, phi, [0.0], freq, T)
ma_psd(sigma, theta, freq, T) = arma_psd(sigma, [0.0], theta, freq, T)
wold_psd(sigma, psi, freq, T) =
    arma_psd(sigma, [0.0], -psi[2:end], freq, T; theta0 = -psi[1])
arma_acf(sigma, psi) =
    sigma^2 .* [dot(psi[1:(end - i)], psi[(1 + i):end]) for i in 0:(length(psi) - 1)]

# == IIR Filter Solver (MA/Wold) == #

# Compute the impulse response function associated with ARMA process arma
function generate_IIR_coefficients(phi, theta, impulse_length = max(length(phi), 1024))
    p = length(phi)

    if isa(impulse_length, Nothing)
        impulse_length = max(p, 1024)
    else
        err_msg = "Impulse length must be greater than number of AR coefficients"
        @assert impulse_length >= p err_msg
    end

    # == Pad theta with zeros at the end == #
    psi_zero = 1.0
    psi = zeros(impulse_length)
    theta = [theta; zeros(impulse_length - length(theta))]
    for j in 1:(impulse_length)
        psi[j] -= theta[j]
        for i in 1:min(j, p)
            psi[j] += phi[i] * (j - i > 0 ? psi[j - i] : psi_zero)
        end
    end
    return [psi_zero; psi[1:(end - 1)]]
end

function extend_IIR_coefficients(psi, phi, impulse_length = max(2 * length(psi), 1024))
    p = length(phi)
    q = length(psi)

    if isa(impulse_length, Nothing)
        impulse_length = max(2 * q, 1024)
    else
        err_msg = "Impulse length must be greater than number of AR coefficients"
        @assert impulse_length >= p err_msg
    end

    # == Pad theta with zeros at the end == #
    psi = [psi; zeros(impulse_length - q)]
    for j in (q + 1):impulse_length
        for i in 1:min(j, p)
            psi[j] += phi[i] * (j - i > 0 ? psi[j - i] : psi_zero)
        end
    end
    return psi
end

function WoldtoARMA(psi, phi; order = length(phi) - 1)
    theta = zeros(order)
    for j in 1:(order)
        theta[j] = -psi[j + 1]
        for i in 1:min(j, order + 1)
            theta[j] += phi[i] * (j - i > 0 ? psi[j - i + 1] : 1.0)
        end
    end
    return theta
end

function WoldDecompACF(psi; λ = 10)
    acf = map(i -> dot(psi[1:(end - i)], psi[(1 + i):end]), 0:(length(psi) - 1))
    return (acf ./ acf[1]) .+ λ * (1.0 - psi[1])
end

@inline function WoldDecompACF_Jac(t, psi; λ = 10)
    N = length(psi)
    jac = Array{eltype(psi)}(undef, length(t), N) # J_jk
    ρ0 = dot(psi[1:end], psi[1:end])
    @inbounds @simd for i in CartesianIndices(jac) # psi_0 starts at index 1
        jac[i] = (i[2] >= i[1]) ? psi[i[2] - i[1] + 1] : 0
        jac[i] += (i[2] + i[1] < N) ? psi[i[2] + i[1] + 1] : 0
        jac[i] -= dot(psi[1:(end - i[1])], psi[(1 + i[1]):end]) * (2 * psi[i[2]])
        jac[i] /= ρ0
        if i[1] == 1
            jac[i] -= λ
        end
    end
    return jac
end

function estimate_IIR_coefficients(acf, psi_guess, coeff_length = length(psi_guess))
    # == Fitting Algorithm for ψ values (useful for CARMAProcess) == #s
    ts = 0:(coeff_length - 1)
    ys = acf[1:coeff_length] ./ acf[1]
    p0 = psi_guess[1:coeff_length]
    psi, _ = so.curve_fit(
        (x, d...) -> WoldDecompACF(collect(d)),
        ts,
        ys,
        p0,
        jac = (x, d...) -> WoldDecompACF_Jac(x, collect(d))
    )
    psi ./= psi[1]

    return psi
end

# == FIR Filter Solver (AR) == #

@inline function YuleWalker!(phi, ρ; ρtrue = ρ)
    n = length(ρ) - 1
    Rmat = zeros(n, length(phi))
    for j in 1:n
        Rmat[j, :] .= ρ[abs.(j .- (1:length(phi))) .+ 1]
    end
    return lsqr!(
        phi,
        Rmat,
        ρtrue[2:end],
        atol = 1e-6,
        btol = 1e-6,
        maxiter = max(maximum(size(Rmat)), 16)
    ) #IterativeSolvers
end

function AR_CharacteristicEqn(phi)
    # == Find Roots & Check Stationary Condition == #
    charpoly = [-reverse(phi); 1.0]
    roots = np.roots(charpoly) # roots of |z|
    if any(abs.(roots) .<= 1.0) # |z| > 1 s.t. |w| < 1 to be stationary
        @warn "AR(p) process is not stationary, must satisfy characteristic equation with |z| > 1. Consider trying different coefficient length."
    end
    return roots
end

@inline function generate_AR_coefficients(acf::Vector, coeff_length; rtol = 1e-3)

    # == Check coeff_length == #
    if !(coeff_length < length(acf))
        coeff_length = length(acf) - 1
    elseif !(coeff_length >= 3)
        coeff_length = 3 # lsqr must be vector; ensures CARMA can be calculated
    end

    # == PACF cutoff == #
    if !isa(rtol, Nothing)
        pacf = zeros(1)
        phi = zeros(1)
        @inbounds for k in 1:coeff_length
            YuleWalker!(phi, acf[1:(k + 1)] / acf[1])
            pacf[k] = phi[end]
            if (k > 1) && isapprox(pacf[k], pacf[k - 1], atol = 100eps(), rtol = rtol)
                # double check that within error bounds
                if abs(pacf[k]) < (1 / sqrt(length(acf)))
                    break
                end
            end
            push!(phi, 0.0)
            push!(pacf, 0.0)
        end

        #remove last push (if needed)
        (length(phi) > coeff_length) || (phi[end] == 0.0) ? pop!(phi) : nothing

        # == No cutoff == #
    else
        phi = zeros(coeff_length)
        phi[1] = 1.0
        YuleWalker!(phi, acf[1:coeff_length] / acf[1])
    end

    # == Find roots == #
    roots = AR_CharacteristicEqn(phi)

    return phi, roots
end

function WoldtoAR(psi; order = length(psi) - 1)
    phi = zeros(order)
    phi[1] = psi[2]
    for i in 2:(order)
        j = i + 1 # psi[1] = 1.0
        phi[i] = psi[j] - dot(phi[1:(i - 1)], psi[(j - 1):-1:(j - length(1:(i - 1)))])
    end
    return phi, AR_CharacteristicEqn(phi)
end

function generate_ARMA_coefficients(acf::Vector, coeff_length, impulse_length; rtol = 1e-3)
    # == AR Coefficients (FIR Filter) == #
    phi, ar_roots = generate_AR_coefficients(acf, coeff_length; rtol = rtol)

    # == MA Coefficients (IIR Filter) == #
    psi = generate_IIR_coefficients(phi, [0.0], impulse_length) # initial psi guess
    # @assert false length(phi)
    psi = estimate_IIR_coefficients(acf, psi, max(128, min(2length(phi), coeff_length)))
    psi = extend_IIR_coefficients(psi, phi, length(acf))
    # sigma_est = sqrt(acf[1] - dot(phi,acf[2:length(phi)+1]))
    # acf_est = arma_acf(sigma_est, psi)
    # psi[2:end] .*= sqrt(acf[1] / acf_est[1])
    theta = WoldtoARMA(psi, phi)

    return ((phi, ar_roots), theta, psi)
end

# == ARMA Structure == #

mutable struct ARMA
    # Paramters
    freq::Union{Vector, Frequencies} # Frequency values
    psd::Vector      # (Discrete) Power Spectral Density
    ar_psd::Vector   # AR Power Spectral Denisty Estimate
    arma_psd::Vector   # ARMA Power Spectral Denisty Estimate
    acf::Vector      # Ensemble Autocorrelation Function
    ar_acf::Vector   # AR ACF Estimate (using Wold's IIR ACF)
    arma_acf::Vector   # ARMA ACF Estimate (using Wold's IIR ACF)
    phi::Vector      # AR parameters phi_1, ..., phi_p
    theta::Vector      # MA parameters theta_1, ..., theta_p
    psi::Vector      # Wold parameters psi_i for i in coeff_length
    ar_roots::Vector # AR roots for each phi_i

    # Extra
    mean::Real       # Mean of noise process
    sigma::Real      # Process standard deviation (√variance)
    ma_sigma::Real
    dt::Real         # Correlation timestep
    timescale::Real  # Timescale
    normalize_sigma::Bool # process variance normalized
    coeff_length::Int # AR coefficient length
    impulse_length::Int # MA Coefficient length
    window_length::Int # Fitting Window
    freq_bin::Real   # Frequency bin size (cycle / (window_length * dt))
    arma_rtol::Real   # AR Coefficient Stopping Criteria
end

@inline function ARMA(
    psd::AbstractVector,
    freq::Union{AbstractVector, Frequencies},
    dt::Real = 1.0,
    timescale::Real = 1.0;
    mean::Real = 0.0,
    normalize_sigma::Bool = true,
    coeff_length::Union{Integer, Nothing} = nothing,
    impulse_length::Union{Integer, Nothing} = nothing,
    rtol = 1e-6
)

    # == Autocorrelation Function == #
    acf = real.(ifft(psd))
    acf = acf[1:(length(psd) ÷ 2 + 1)]

    # == Find ARMA Coefficients == #
    (coeff_length === nothing) ? coeff_length = length(acf) : nothing
    (impulse_length === nothing) ? impulse_length = length(acf) : nothing
    (phi, ar_roots), theta, psi =
        generate_ARMA_coefficients(acf, coeff_length, impulse_length; rtol = rtol)

    # == Power Spectrum & Noise estimation == #
    ar_ps = ar_psd(1.0, phi, freq, dt * timescale)
    ma_ps = arma_psd(1.0, phi, theta, freq, dt * timescale)

    ar_sigma = sqrt(sqrt(sum((psd ./ ar_ps) .^ 2) / length(psd)))
    # ar_sigma = sqrt(acf[1] - dot(phi,acf[2:length(phi)+1])) # est from Yule-Walker
    ar_ps .*= ar_sigma^2

    ma_sigma = sqrt(sqrt(sum((psd ./ ma_ps) .^ 2) / length(psd)))
    # ma_sigma = sqrt(acf[1]) / sqrt(dot(psi,psi)) # est from acf fitting
    ma_ps .*= ma_sigma^2

    # == AR/ARMA Fitted PSD & ACF Est. == #
    # psi_tmp = generate_IIR_coefficients(phi, theta, length(acf))
    ma_acf = arma_acf(ma_sigma, psi)

    psi_tmp = generate_IIR_coefficients(phi, [0.0], length(acf))
    ar_acf = arma_acf(ar_sigma, psi_tmp)

    # == Normalize Noise variance == #
    if normalize_sigma # process with unit variance
        ar_sigma /= sqrt(ar_acf[1])
        ma_sigma /= sqrt(ma_acf[1])
        ar_ps ./= ar_acf[1]
        ma_ps ./= ma_acf[1]
        ar_acf ./= ar_acf[1]
        ma_acf ./= ma_acf[1]
        psd ./= acf[1]
        acf ./= acf[1]
    end

    return ARMA(
        freq,
        psd,
        ar_ps,
        ma_ps,
        acf,
        ar_acf,
        ma_acf,
        phi,
        theta,
        psi,
        ar_roots,
        mean,
        ar_sigma,
        ma_sigma,
        dt,
        timescale,
        normalize_sigma,
        coeff_length,
        impulse_length,
        length(freq),
        abs(freq[2] - freq[1]),
        rtol
    )
end

# # https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1171974&tag=1
# # solve for residual power spectrum
# q = length(arma.phi) - 1

# res_psd = (arma.psd ./ arma.ar_psd)
# freq = arma.freq
# Δt = arma.dt * arma.timescale
# Zmat = freq .* transpose(-q:q) # ck is symmetric
# Zmat = @. exp(1im*2pi* Zmat *Δt)

# # == Solve for θ values == #
# ys = res_psd
# p0 = ones(q+1)
# ck_tmp = Array{Float64}(undef, q+1)

# function ck_values(ts, theta, ck_tmp; q = length(theta)-1) # q = length(theta) minus theta0 term
#     for (idx,k) in enumerate(0:q)
#         ck_tmp[idx] = dot(theta[1:(q-k)+1],theta[k+1:q+1])
#     end
#     # return real.(Zmat * (ck_tmp .+ 100(-1.0 - theta[1])))
#     return real.(Zmat * [reverse(ck_tmp[2:end]); ck_tmp])
#     # return real.(Zmat * ck_tmp)
# end

# θ, _ = so.curve_fit((x,theta...)->ck_values(x, collect(theta), ck_tmp), ts, ys, p0)
# scale = -θ[1]
# θ./=scale

# plot((ys .- scale^2 * ck_values(ts, θ, ck_tmp))./ys) |> display

# arma_psd = GaussianNoiseProcess.arma_psd(arma.sigma * scale, arma.phi, θ[2:end], arma.freq, arma.dt*arma.timescale)
# # arma_psd .*= maximum(arma.psd) / maximum(arma_psd)

# plot(arma.psd[1:length(arma.psd)÷2+1], label="Model")
# plot!(arma_psd[1:length(arma.psd)÷2+1], linewidth=2, label="ARMA Fit")
# plot!(arma.ar_psd[1:length(arma.psd)÷2+1], label="AR Fit")
# plot!(xscale=:log, yscale=:log)
