using Parameters
using DPMNeal3
using Random
using Test
include("utils.jl")

struct MyData
    x::Vector{Int}
    y::Vector{Float64}
end

@testset "GenericBlock" begin
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

@testset "update_α!" begin
    N  = 4
    K0 = 2
    rng = MersenneTwister(1)
    gb = GenericBlock(rng, N; K0 = K0)
    DPMNeal3.update_α!(rng, gb)
    @test gb.α[] < Inf
    @test gb.α[] > 0
end

@testset "SpecificBlock" begin
    G = 4
    sb = SpecificBlock(G; v0 = 2, r0 = 3, u0 = 3.0, s0 = 9.0)
    @test sb.G  == 4
    @test sb.v0 == 2
    @test sb.r0 == 3
    @test sb.u0 == 3.0
    @test sb.s0 == 9.0
    @test sb.v1 == [2 * ones(G)]
    @test sb.r1 == [3 * ones(G)]
    @test sb.u1 == [3 * ones(G)]
    @test sb.s1 == [9 * ones(G)]
    @test sb.γ  == zeros(Int, G)
end

@testset "resize!" begin
    G = 4
    n = 2
    sb = SpecificBlock(G; v0 = 2, r0 = 3, u0 = 3.0, s0 = 9.0)
    resize!(sb, n)
    @test length(sb.v1) == n
    @test length(sb.r1) == n
    @test length(sb.u1) == n
    @test length(sb.s1) == n
end

@testset "update_sb! (1)" begin
    rng = MersenneTwister(1)
    data = MyData([1], [1.0])
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    sb = SpecificBlock(1; v0, r0, u0, s0)
    gb = GenericBlock(rng, 1)
    update_sb!(sb, gb, data)
    @test sb.v1[1][1] ≈ 2.0
    @test sb.v1[2][1] ≈ 1.0
    @test sb.r1[1][1] ≈ 2.0
    @test sb.r1[2][1] ≈ 1.0
    @test sb.u1[1][1] ≈ 0.5
    @test sb.u1[2][1] ≈ 0.0
    @test sb.s1[1][1] ≈ 1.5
    @test sb.s1[2][1] ≈ 1.0
end

@testset "update_sb! (2)" begin
    rng = MersenneTwister(1)
    data = MyData([1], [1.0])
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    sb = SpecificBlock(1; v0, r0, u0, s0)
    gb = GenericBlock(rng, 1)
    update_sb!(sb, gb, data)
    update_sb!(sb, gb, data, 1, 1, 2)
    @test sb.v1[1][1] ≈ 1.0
    @test sb.v1[2][1] ≈ 2.0
    @test sb.r1[1][1] ≈ 1.0
    @test sb.r1[2][1] ≈ 2.0
    @test sb.u1[1][1] ≈ 0.0
    @test sb.u1[2][1] ≈ 0.5
    @test sb.s1[1][1] ≈ 1.0
    @test sb.s1[2][1] ≈ 1.5
end

@testset "logh" begin
    rng = MersenneTwister(1)
    data = MyData([1], [1.0])
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    sb = SpecificBlock(1; v0, r0, u0, s0)
    gb = GenericBlock(rng, 1)
    update_sb!(sb, gb, data)
    @test logh(sb, gb, data, 1) ≈ (
        0.5 * 1.0 * log(1.0) -
        0.5 * 2.0 * log(1.5) +
        loggamma(2.0 / 2) -
        loggamma(1.0 / 2) +
        0.5 * log(1.0 / 2.0) -
        0.5 * log(π)
    )
end

@testset "logq" begin
    rng = MersenneTwister(1)
    data = MyData([1, 1], [1.0, 0.0])
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    sb = SpecificBlock(1; v0, r0, u0, s0)
    gb = GenericBlock(rng, 2)
    update_sb!(sb, gb, data)
    @test sb.v1[1][1] ≈ 3.0
    @test sb.r1[1][1] ≈ 3.0
    @test sb.u1[1][1] ≈ 1/3
    @test sb.s1[1][1] ≈ 5/3

    @test logq(sb, gb, data, 2, 1) ≈ (
        0.5 * 2 * log(1.5) -
        0.5 * 3 * log(5/3) +
        loggamma(3/2) -
        loggamma(2/2) +
        0.5 * log(2/3) -
        0.5 * log(π)
    )
end

# x = StatsBase.denserank(x)
