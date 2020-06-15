import QuantumOptics
const qo = QuantumOptics
using Test, IonSim


# setup system
C = Ca40(["S-1/2", "D-1/2"])
chain = LinearChain(
    ions=[C, C], com_frequencies=(x=2,y=2,z=1), selected_modes=(x=[1], y=[], z=[1])
)
T = Trap(configuration=chain)
modes = get_vibrational_modes(chain) 


@testset "VibrationalMode operators" begin
    # test creation of VibrationalMode functions by comparison with equiv. QO functions
    fb = qo.FockBasis(10)
    @test create(modes[1]).data == qo.create(fb).data
    @test destroy(modes[1]).data == qo.destroy(fb).data
    @test number(modes[1]).data == qo.number(fb).data
    n = rand(0:10)
    @test fockstate(modes[1], n).data == qo.fockstate(fb, n).data
    # @test displace

    # test that mean excitation of thermalstate is as expected
    modes[1].N = 200
    n̄ = abs(2randn())
    @test expect(number(modes[1]), thermalstate(modes[1], n̄)) ≈ n̄

    # test coherentstate matches QO results
    α = 10*(randn() + im*randn()) 
    coherentstate(modes[1], α).data == qo.coherentstate(fb, α).data

    # test coherenthermalstate

end


@testset "Ion operators" begin
    # test that ionstate constructs the appropriate state for a single ion
    @test ionstate(C, "S-1/2").data == ionstate(C, 1).data == ComplexF64[1; 0]
    @test ionstate(C, "D-1/2").data == ionstate(C, 2).data == ComplexF64[0; 1]

    # test ionstate for an IonConfiguration input
    @test ionstate(chain, "S-1/2", "D-1/2").data == kron(ComplexF64[0; 1], ComplexF64[1; 0])
    @test ionstate(chain, 1, 2).data == kron(ComplexF64[0; 1], ComplexF64[1; 0])
    
    # test ionstate for an Trap input
    @test ionstate(T, "S-1/2", "D-1/2").data == kron(ComplexF64[0; 1], ComplexF64[1; 0])

    # test sigma(ion::Ion, ψ1::T, ψ2::T) where {T<:Union{String,Int}}
    @test sigma(C, "S-1/2", "D-1/2").data == sigma(C, 1, 2).data == ComplexF64[0 1; 0 0]

    # test ionprojector for IonConfiguration input
    ψ = ionprojector(chain, "S-1/2", "D-1/2", only_ions=true)
    @test  ψ.data == kron(ComplexF64[0; 1] * ComplexF64[0; 1]', ComplexF64[1; 0] * ComplexF64[1; 0]')
    @test ionprojector(chain, "S-1/2", "D-1/2") == ψ ⊗ one(modes[1]) ⊗ one(modes[2])
    @test ψ == ionprojector(T, "S-1/2", "D-1/2", only_ions=true)
end