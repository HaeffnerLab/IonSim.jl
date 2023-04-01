using Test, IonSim

@testset "termsymbols" begin
    # check string macros
    @test LSTerm(n=4, l=2, s=1//2, j=3//2) == ls"4²D_3/2"
    @test J₁KTerm(k=3//2, s=1, j=1//2)  == jk"³[3/2]_1/2"

    # check format checks
    # try 
    #     ls"4²D"
    # end

    # check print not broken
    print(ls"4²D_3/2")
    println(ls"4²D_3/2")
end