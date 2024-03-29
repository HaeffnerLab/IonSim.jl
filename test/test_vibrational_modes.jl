using Test, IonSim
using Suppressor

@testset "vibrational_modes -- VibrationalMode" begin
    # setup system
    vm = VibrationalMode(1, [1, 1], δν=1)

    # make sure δν can appropriately handle constant or function
    t = 0:1:100
    ones = [1 for _ in t]
    @test vm.δν.(t) == ones
    @test vm._cnst_δν
    frequency_fluctuation!(vm, sin)
    @test vm.δν.(t) == sin.(t)
    @test !vm._cnst_δν
    vm = VibrationalMode(1, [1, 1], δν=t -> t)
    @test !vm._cnst_δν

    # test ==
    vm1 = VibrationalMode(1, [1, 1], δν=1)
    @test vm == vm1

    # test that updating cutoff also appropriately updates shape field
    modecutoff!(vm, 20)
    @test vm.shape == [21]

    # # make sure user can't directly update shape field
    # vm.shape = 100
    # @test vm.shape == [21]
    # @test vm.N == 20

    # make sure user can't set field N to a float or negative value
    @test_throws MethodError modecutoff!(vm, 1.0)
    @test_throws AssertionError modecutoff!(vm, -1)

    # just make sure no errors on print
    @suppress print(vm)
    @suppress show(vm)

    # # make sure user can't directly alter :modestructure or :axis after initialization
    # ms = vm.modestructure
    # vm.modestructure = "something"
    # @test vm.modestructure == ms
    # axis = vm.axis
    # vm.axis = "something"
    # @test vm.axis == axis

    # test indexing of VibrationalMode, which should return a fockstate
    modecutoff!(vm, 1)
    @test vm[0].data == ComplexF64[1; 0]
    @test vm[1].data == ComplexF64[0; 1]
end
