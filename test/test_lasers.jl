using Test, IonSim
using Suppressor

@testset "lasers -- Laser" begin
    # should not be able to set unnormalized k-vectors or polarizations
    @test_throws AssertionError Laser(ϵ = (x = 2, y = 0, z = 0))
    @test_throws AssertionError Laser(k = (x = 2, y = 0, z = 0))
    # The below test is commented out because we have temporarily removed the non-orthogonal warning for convenience (see lasers.jl)
    # should not be able to initialize with non-orthogonal ϵ̂, k̂
    # @test_throws AssertionError Laser(ϵ=(x=1, y=0, z=0), k=(x=1, y=0, z=0))
    # pointing cannot be overspecified 
    @test_throws AssertionError Laser(pointing = [(1, 0.5), (1, 0.5)])
    # pointing magnitude must be within [0,1]
    @test_throws AssertionError Laser(pointing = [(1, 2.0)])

    # should be able to set L.E, L.ϕ to a constant Real value 
    t = 0:1:100
    L = Laser(E = 1, ϕ = 1)
    ones = [1 for _ in t]
    # should be able to set L.E, L.ϕ to a function of time
    @test L.E.(t) == ones
    @test L.ϕ.(t) == ones
    L = Laser(E = sin, ϕ = sin)
    @test L.E.(t) == sin.(t)
    @test L.ϕ.(t) == sin.(t)

    # test that normalization is enforced for altered L.ϵ/L.k
    @test_throws AssertionError polarization!(L, (x = 2, y = 0, z = 0))
    @test_throws AssertionError wavevector!(L, (x = 2, y = 0, z = 0))
    # test that pointing constraints are enforced for altered values
    @test_throws AssertionError pointing!(L, [(1, 0.5), (1, 0.5)])
    @test_throws AssertionError pointing!(L, [(1, 2.0)])
    L.λ = 1
    @test L.λ == 1

    # just make sure print(L) isn't broken
    @suppress print(L)

    # Commeneted out until non-orthogonal warnings are updated
    # # test non-orthogonal warnings
    # @test_logs (:warn, "!(ϵ ⟂ k)") L.ϵ = (x̂+ŷ+ẑ)/√3
    # @test_logs (:warn, "!(ϵ ⟂ k)") L.k = (x̂+ŷ+ẑ)/√3

    # test comparison
    L2 = copy(L)
    @test L == L2
    efield!(L, 7)
    @test L ≠ L2
end
