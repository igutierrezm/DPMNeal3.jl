using DPMNeal3
using Random
using Test

@testset "GenericBlock.jl" begin
    rng = MersenneTwister(1)
    @test_throws AssertionError GenericBlock(rng, 10; K0 = 20)
    @test_throws AssertionError GenericBlock(rng, 10; K0 = -1)
    @test_throws AssertionError GenericBlock(rng, 10; a0 = -1)
    @test_throws AssertionError GenericBlock(rng, 10; b0 = -1)
    @test_throws TypeError      GenericBlock(rng, 10; K0 = .5)
    @test_throws MethodError    GenericBlock(rng, .5)
    N  = 4
    K0 = 2
    gb = GenericBlock(rng, N; K0 = K0)
    @test unique(gb.d) == collect(1:K0)
    @test length(gb.d) == N
    @test gb.A == Set(1:K0)
    @test gb.P == [K0 + 1]
    @test gb.K[] == K0
end
