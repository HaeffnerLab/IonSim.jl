import QuantumOptics
const qo = QuantumOptics
using Test, IonSim
using Suppressor

@suppress_err begin

    # setup system
    C = Ca40([("S1/2", -1 // 2), ("D5/2", -1 // 2)])
    chain = LinearChain(
        ions = [C, C],
        com_frequencies = (x = 2u"Hz", y = 2u"Hz", z = 1u"Hz"),
        vibrational_modes = (x = [1], y = [], z = [1])
    )
    T = Trap(configuration = chain)
    modes = get_vibrational_modes(chain)

    @testset "operators -- VibrationalMode operators" begin
        # test creation of VibrationalMode functions by comparison with equiv. QO functions
        fb = qo.FockBasis(10)
        @test create(modes[1]).data == qo.create(fb).data
        @test destroy(modes[1]).data == qo.destroy(fb).data
        @test number(modes[1]).data == qo.number(fb).data
        n = rand(0:10)
        @test fockstate(modes[1], n).data == qo.fockstate(fb, n).data

        # displacement operator
        α = rand(0:1e-3:7) + im * rand(0:1e-3:7)
        @test displace(modes[1], α).data ≈ qo.displace(fb, α).data
        fb2 = qo.FockBasis(200)
        modes[1].N = 200
        i = rand(1:5)
        j = rand(1:5)
        @test displace(modes[1], α, method = "analytic").data[i, j] ≈
              qo.displace(fb2, α).data[i, j] atol = 1e-6

        # test that mean excitation of thermalstate is as expected
        modes[1].N = 500
        n̄ = abs(2randn())
        @test expect(number(modes[1]), thermalstate(modes[1], n̄)) ≈ n̄
        @test expect(number(modes[1]), thermalstate(modes[1], n̄, method = "analytic")) ≈ n̄

        # test coherentstate matches QO results
        α = 10 * (randn() + im * randn())
        coherentstate(modes[1], α).data == qo.coherentstate(fb, α).data

        # test coherenthermalstate
        N = 500
        modes[1].N = N
        n̄ = rand(0:1e-6:10)
        @test coherentthermalstate(modes[1], n̄, 0, method = "analytic").data ≈
              thermalstate(modes[1], n̄).data
        @test coherentthermalstate(modes[1], 0, α, method = "analytic").data ≈
              dm(coherentstate(modes[1], α)).data rtol = 1e-3 * N^2
        @test coherentthermalstate(modes[1], n̄, 0).data ≈ thermalstate(modes[1], n̄).data
        @test coherentthermalstate(modes[1], 0, α).data ≈
              dm(coherentstate(modes[1], α)).data rtol = 1e-3 * N^2

        # shouldn't be able to have a mean phonon occupation greater than Hilbert space dimension
        @test_throws AssertionError coherentthermalstate(modes[1], N + 1, 0)
        @test_throws AssertionError coherentthermalstate(modes[1], 0, N + 1)
        @test_throws AssertionError coherentstate(modes[1], N + 1)
        @test_throws AssertionError thermalstate(modes[1], N + 1)
        @test_throws AssertionError displace(modes[1], N + 1)
    end

    @testset "operators -- Ion operators" begin
        # test that ionstate constructs the appropriate state for a single ion
        @test ionstate(C, ("S1/2", -1 // 2)).data == ionstate(C, 1).data == ComplexF64[1; 0]
        @test ionstate(C, ("D5/2", -1 // 2)).data == ionstate(C, 2).data == ComplexF64[0; 1]

        # test ionstate for an IonConfiguration input
        @test ionstate(chain, ("S1/2", -1 // 2), ("D5/2", -1 // 2)).data ==
              kron(ComplexF64[0; 1], ComplexF64[1; 0])
        @test ionstate(chain, 1, 2).data == kron(ComplexF64[0; 1], ComplexF64[1; 0])

        # test ionstate for an Trap input
        @test ionstate(T, ("S1/2", -1 // 2), ("D5/2", -1 // 2)).data ==
              kron(ComplexF64[0; 1], ComplexF64[1; 0])

        # test sigma(ion::Ion, ψ1::T, ψ2::T) where {T<:Union{String,Int}}
        @test sigma(C, ("S1/2", -1 // 2), ("D5/2", -1 // 2)).data ==
              sigma(C, 1, 2).data ==
              ComplexF64[0 1; 0 0]
        # test sigma(ion::Ion, ψ1::T<:Union{String,Int})
        @test sigma(C, ("S1/2", -1 // 2)).data == sigma(C, 1).data == ComplexF64[1 0; 0 0]

        # test ionprojector for IonConfiguration input
        ψ = ionprojector(chain, ("S1/2", -1 // 2), ("D5/2", -1 // 2), only_ions = true)
        @test ψ.data == kron(
            ComplexF64[0; 1] * ComplexF64[0; 1]',
            ComplexF64[1; 0] * ComplexF64[1; 0]'
        )
        @test ionprojector(chain, ("S1/2", -1 // 2), ("D5/2", -1 // 2)) ==
              ψ ⊗ one(modes[1]) ⊗ one(modes[2])
        @test ψ == ionprojector(T, ("S1/2", -1 // 2), ("D5/2", -1 // 2), only_ions = true)
    end

    @testset "operators -- internal functions" begin
        # test _pf(s, n, m)
        s = rand(1:12)
        n = rand(1:s)
        m = rand(1:s)
        v1 = IonSim._pf(s, n, m)
        v2 = 1im^n * (-1im)^m * factorial(s) / ((s + 1) * √(factorial(m) * factorial(n)))
        isapprox(v1, v2)

        # test He(n) gives the correct Hermite polynomial for order 10
        IonSim._He(10) == [-945, 0, 4725, 0, -3150, 0, 630, 0, -45, 0, 1]

        # test fHe(x, He)
        IonSim._fHe(1, 10) == sum(IonSim._He(10))

        # _alaguerre
        L32 =
            (1 / 6) * (
                -2^3 + 3 * (2 + 3) * 2^2 - 3 * (2 + 2) * (2 + 3) * 2 +
                (2 + 1) * (2 + 2) * (2 + 3)
            )
        @test IonSim._alaguerre(2, 3, 2) ≈ L32

        # _Dnm
        ξ = im * exp(2π * im)
        d = displace(modes[1], ξ).data
        diff = 0.0
        for i in 1:100, j in 1:100
            diff += abs(d[i, j] - IonSim._Dnm(ξ, i, j))
        end
        @test diff < 100  # <1% difference of L1 norm
        # Note: displace() is an approximation, whereas _Dnm should not be
    end
end  # end suppress
