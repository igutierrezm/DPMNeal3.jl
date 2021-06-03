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
    @test sb.v1 == [zeros(G)]
    @test sb.r1 == [zeros(G)]
    @test sb.u1 == [zeros(G)]
    @test sb.s1 == [zeros(G)]
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
    N = 10 # sample size
    F = 2 # number of factors
    J = [3, 2] # levels per factor
    x = [zeros(Int, F) for i in 1:N]
    for i = 1:N, f = 1:F
        x[i][f] = rand(1:J[f])
    end
    x = StatsBase.denserank(x)
    y = randn(N)
    data = MyData(x, y)
    G = length(unique(x))
    sb = SpecificBlock(G)
    gb = GenericBlock(rng, N)
    update_sb!(sb, gb, data)
end

@testset "update_sb! (2)" begin
    rng = MersenneTwister(1)
    N = 10 # sample size
    F = 2 # number of factors
    J = [3, 2] # levels per factor
    x = [zeros(Int, F) for i in 1:N]
    for i = 1:N, f = 1:F
        x[i][f] = rand(1:J[f])
    end
    x = StatsBase.denserank(x)
    y = randn(N)
    data = MyData(x, y)
    G = length(unique(x))
    sb = SpecificBlock(G)
    gb = GenericBlock(rng, N)
    update_sb!(sb, gb, data)
    update_sb!(sb, gb, data, 1, 1, 2)
end

@testset "log_pl" begin
    rng = MersenneTwister(1)
    N = 10 # sample size
    F = 2 # number of factors
    J = [3, 2] # levels per factor
    x = [zeros(Int, F) for i in 1:N]
    for i = 1:N, f = 1:F
        x[i][f] = rand(1:J[f])
    end
    x = StatsBase.denserank(x)
    y = randn(N)
    data = MyData(x, y)
    G = length(unique(x))
    sb = SpecificBlock(G)
    gb = GenericBlock(rng, N)
    update_sb!(sb, gb, data)
    update_sb!(sb, gb, data, 1, 1, 2)
    log_pl(sb, gb, data, 1, 1)
end