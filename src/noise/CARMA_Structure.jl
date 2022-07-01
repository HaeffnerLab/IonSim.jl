# == CARMA poles == #
function ARMAtoCARMAPoles(ar_roots, arma_τ)
  ar_ρ = abs.(ar_roots)
  ar_θ = @. atan(imag(ar_roots) / real(ar_roots))
  return @. -(log(ar_ρ) + 1im * ar_θ) / arma_τ
end

# == Continuous ACF Fitting Functions == #
mutable struct CARMA_ACFSolver{T1,T2,T3,T4,T5,T6}
  acf::T1
  fit_length::T2
  Λ::T3
  ck::T4
  λ::T4
  exp_λτ::T5
  jac::T3
  carma_sigma::T6
  function CARMA_ACFSolver(acf, λr, sigma, dt)
    fit_length = min(length(λr),length(acf))
    denom_inv = 1.0 ./ (λr .+ transpose(λr))
    t = (0:length(acf)-1) .* dt
    expλt = -exp.(t * transpose(λr))
    expλt .*= sigma^2
    new{typeof(acf), typeof(fit_length), typeof(denom_inv), typeof(λr), typeof(expλt), typeof(sigma)}(acf,fit_length,denom_inv,similar(λr),λr,expλt,similar(denom_inv),sigma)
  end
end

function ck_ReIm_values(X,α)
  α = α[1:length(α)÷2] .+ 1im.*α[length(α)÷2+1:end]
  X.ck .= transpose(transpose(α) * X.Λ)
  X.ck .*= α
  # X.ck .+= 2*(X.acf[1] - X.carma_sigma^2 * sum(-X.ck))
  return [real(X.ck); imag(X.ck)]
end

function ck_values!(X,α)
  X.ck .= transpose(transpose(α) * X.Λ)
  X.ck .*= α
end

# function jac_ck_ReIm_values(X,α)
#   α = α[1:length(α)÷2] .+ 1im.*α[length(α)÷2+1:end]
#   for i in CartesianIndices(X.jac)
#     X.jac[i] = α[i[1]] / (X.λ[i[2]] + X.λ[i[1]])
#     if i[1] == i[2]
#       X.jac[i] += sum(α .* (1.0 ./ (X.λ .+ X.λ[i[2]])))
#     end
#   end
#   @assert false size(α)
#   return hcat(real(X.jac), imag(X.jac))
# end


@inline function (X::CARMA_ACFSolver)(psi, arma_sigma; p0=nothing, fit_args=nothing)
  # == Initialize guess == #
  if isa(p0, Nothing)
    p0 = ones(2*length(ck))
    X.ck .= zeros(X.fit_length)
  elseif eltype(p0) <: Complex
    ck_values!(X,p0)
    p0 = [real(p0);imag(p0)]
  end

  # == Solve for c_k values == #
  lsqr!(X.ck, X.exp_λτ[1:X.fit_length,:], X.acf[1:X.fit_length],
        atol=1e-9,btol=1e-9, maxiter=maximum(size(X.exp_λτ)))

  # == Solve for α values == #
  if isa(fit_args,Nothing)
    carma_α, _ = so.curve_fit((x,d...)->ck_ReIm_values(X, collect(d)),
    [1:X.fit_length], [real(X.ck);imag(X.ck)], p0) #, jac=(x,d...)->jac_ck_ReIm_values(X,collect(d)))
  else
    carma_α, _ = so.curve_fit((x,d...)->ck_ReIm_values(X, collect(d)),
    [1:X.fit_length], [real(X.ck);imag(X.ck)], p0, fit_args...)
  end

  return carma_α[1:length(X.ck)] .+ 1im .* carma_α[length(X.ck)+1:end]
end

# == ARMA to CARMA Solver == #
@inline function ARMAtoCARMASolver(arma_τ, acf, ar_phi, psi, carma_poles, arma_sigma; fit_args=nothing)
  # == Setup solvers == #
  carma_sigma = arma_sigma / sqrt(arma_τ)
  Solver = CARMA_ACFSolver(acf, carma_poles, carma_sigma, arma_τ)

  # == Impulse invariance (weaker statement) -> faster; less accurate for smaller # of parameters == #
  μmatrix = -Solver.exp_λτ[1:min(length(acf),length(psi)), :] ./ carma_sigma^2
  μTmatrix = μmatrix'

  carma_α = lsqr(μTmatrix * μmatrix, μTmatrix * psi,
            atol=1e-9, btol=1e-9, maxiter=maximum(size(Solver.exp_λτ))) #iterative solvers

  # == ACF invariance (stronger statement; slower numerically; less prone to numerical error) == #
  if (length(ar_phi) <= 256)
    carma_α = Solver(psi, arma_sigma; p0=carma_α, fit_args=fit_args)
  end

  # == Best Fit Parameters == #
  ck_values!(Solver, carma_α)
  carma_acf = Solver.exp_λτ * Solver.ck
  carma_α .*= sqrt(acf[1] / abs(carma_acf[1]))
  carma_acf .*= acf[1] / abs(carma_acf[1])
  carma_psd = real(fft([carma_acf;reverse(carma_acf[2:end-1])]))

  return carma_α, carma_sigma, carma_psd, real(carma_acf)
end

# mutable struct CARMA_ACFSolver{T1,T2,T3}
#   sqmat::T1
#   denom_inv::T1
#   res1::T2
#   res2::T2
#   exp_λτ::T2
#   σ::T3
#   function CARMA_ACFSolver(λr, sigma, t, dt)
#     denom_inv = 1.0 ./ (λr .+ transpose(λr))
#     expλt = exp.((t .* dt) * transpose(λr))
#     new{typeof(denom_inv), typeof(expλt), typeof(sigma)}(similar(denom_inv), denom_inv, similar(expλt),similar(expλt), expλt, sigma)
#   end
# end

# @inline function cont_autocorr(X,α)
#   α = α[1:length(α)÷2] .+ 1im.*α[length(α)÷2+1:end]
#   X.sqmat .= (α * transpose(α)) # d[i] * d[j] (row = i, column = j)
#   X.sqmat .*= X.denom_inv # κ[i] + κ[j]
#   X.res1 .= X.exp_λτ * X.sqmat
#   return -X.σ^2 .* sum(X.res1, dims=2)[:]
# end

# @inline function jac_autocorr(X,α)
#   α = α[1:length(α)÷2] .+ 1im.*α[length(α)÷2+1:end]

#   X.sqmat .= α .* X.denom_inv # Γik

#   dγdc1 = X.res1
#   dγdc1 .= X.exp_λτ .* sum(X.sqmat,dims=1)
#   dγdc1 .+= X.exp_λτ * X.sqmat
#   # dγdc2 = 1im .* dγdc1

#   return -X.σ^2 .* hcat([real.(dγdc1); imag.(dγdc1)],[-imag.(dγdc1); real.(dγdc1)])
# end


    # # == Setup == #
    # ts = 0:1:min(2*length(ar_phi)-1,length(acf)-1)
    # ys = acf[1:length(ts)]
    # p0 = ones(Float64, length(ts))
    # Solver = CARMA_ACFSolver(carma_poles, arma_sigma / sqrt(arma_τ), ts, arma_τ)

    # # == curve fitting via scipy curve_fit == #
    # if isa(fit_args,Nothing)
    #     carma_α, _ = so.curve_fit((x,d...)->Solver(x,collect(d)), [ts;ts], [ys;zero(ys)],
    #                 p0, jac=(x,d...)->jac_autocorr(Solver, collect(d)))
    # else
    #     carma_α, _ = so.curve_fit((x,d...)->Solver(x,collect(d)), [ts;ts], [ys;zero(ys)],
    #                 p0, jac=(x,d...)->jac_autocorr(Solver, collect(d)), fit_args...)
    # end
    # carma_α = carma_α[1:length(ar_phi)] .+ 1im .* carma_α[length(ar_phi)+1:end]

        # # == Setup == #
        # ts = 0:1:(length(psi)-1)
        # ts *= arma_τ

        # # == Impulse invariance via lsqr solver == #
        # μmatrix = Array{ComplexF64}(undef, length(ts), length(carma_poles))
        # for i in CartesianIndices(μmatrix)
        #   μmatrix[i] = exp(carma_poles[i[2]] * ts[i[1]])
        # end
        # μTmatrix = μmatrix'
        # carma_α = lsqr(μTmatrix * μmatrix, μTmatrix * psi, atol=1e-9, btol=1e-9) #iterative solvers
        # # carma_α = lsqr(μmatrix, psi, atol=1e-9, btol=1e-9) #iterative solvers

        # ts = 0:1:(length(acf)-1)
        # carma_sigma = arma_sigma/sqrt(arma_τ) # estimate discrete -> continuous
        # Solver = CARMA_ACFSolver(carma_poles, carma_sigma, ts, arma_τ)
        # carma_acf = Solver([ts;ts], [real.(carma_α);imag.(carma_α)]) |> sol -> sol[ts .+ 1] .+ 1im * sol[(ts[end]+2):end]
        # # carma_sigma *= sqrt(acf[1] / real(carma_acf[1]))
        # # carma_acf .*= acf[1] / real(carma_acf[1])

        # carma_psd = real.(fft([carma_acf; reverse(carma_acf[2:end-1])]))
        # # carma_sigma *= sqrt(np.mean(abs.(psd ./ carma_psd))) # approximate scale factor
