using Parameters
using DPMNeal3
using Random
using Test
using StatsBase
using Statistics
using SpecialFunctions
using Lazy: @forward

# import DPMNeal3: update_suffstats!, logpredlik
# include("utils.jl")

struct Data
    x::Vector{Int}
    y::Vector{Float64}
end

@testset "DPM" begin
    rng = MersenneTwister(1)
    @test_throws AssertionError DPM(rng, 10; K0 = 20)
    @test_throws AssertionError DPM(rng, 10; K0 = -1)
    @test_throws AssertionError DPM(rng, 10; a0 = -1)
    @test_throws AssertionError DPM(rng, 10; b0 = -1)
    @test_throws TypeError      DPM(rng, 10; K0 = .5)
    @test_throws MethodError    DPM(rng, .5)
    N  = 4
    K0 = 2
    gb = DPM(rng, N; K0 = K0)
    @test unique(gb.d) == collect(1:K0)
    @test length(gb.d) == N
    @test gb.A == Set(1:K0)
    @test gb.P == Set(K0 + 1)
    @test gb.K[1] == K0
end

@testset "NormalDPM (1)" begin
    rng = MersenneTwister(1)
    m = NormalDPM(rng, 5)
    @test m.v0 == 2
    @test m.r0 == 1
    @test m.u0 == 0.0
    @test m.s0 == 1.0
    @test m.v1 == [2]
    @test m.r1 == [1]
    @test m.u1 == [0]
    @test m.s1 == [1]
end

@testset "NormalDPM (2)" begin
    rng = MersenneTwister(1)
    m = NormalDPM(rng, 5; v0 = 2.0, r0 = 3.0, u0 = 3.0, s0 = 9.0)
    @test m.v0 == 2
    @test m.r0 == 3
    @test m.u0 == 3.0
    @test m.s0 == 9.0
    @test m.v1 == [2]
    @test m.r1 == [3]
    @test m.u1 == [3]
    @test m.s1 == [9]
end

@testset "update_α!" begin
    rng = MersenneTwister(1)
    m = NormalDPM(rng, 5)
    DPMNeal3.update_α!(rng, m)
    @test dp_mass(m) < Inf
    @test dp_mass(m) > 0
end

@testset "add_cluster!" begin
    N, K0 = 5, 2
    rng = MersenneTwister(1)
    m = NormalDPM(rng, N; K0)
    DPMNeal3.add_cluster!(m)
    @test length(m.v1) == 2
    @test length(m.r1) == 2
    @test length(m.u1) == 2
    @test length(m.s1) == 2
end

@testset "update_suffstats! (1)" begin
    N = 1
    y = ones(N)
    rng = MersenneTwister(1)
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    m = NormalDPM(rng, N; v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m, y)
    @test m.v1[1] ≈ 2.0
    @test m.v1[2] ≈ 1.0
    @test m.r1[1] ≈ 2.0
    @test m.r1[2] ≈ 1.0
    @test m.u1[1] ≈ 0.5
    @test m.u1[2] ≈ 0.0
    @test m.s1[1] ≈ 1.5
    @test m.s1[2] ≈ 1.0
end

@testset "update_sb! (2)" begin
    N = 1
    y = ones(N)
    rng = MersenneTwister(1)
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    m = NormalDPM(rng, N; v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m, y)
    DPMNeal3.update_suffstats!(m, y, 1, 1, 2)
    @test m.v1[1] ≈ 1.0
    @test m.v1[2] ≈ 2.0
    @test m.r1[1] ≈ 1.0
    @test m.r1[2] ≈ 2.0
    @test m.u1[1] ≈ 0.0
    @test m.u1[2] ≈ 0.5
    @test m.s1[1] ≈ 1.0
    @test m.s1[2] ≈ 1.5
end

@testset "logpredlik (empty clusters)" begin
    N = 1
    y = ones(N)
    rng = MersenneTwister(1)
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    m = NormalDPM(rng, N; v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m, y)
    @test DPMNeal3.logpredlik(m, y, 1, first(passive_clusters(m))) ≈ (
        0.5 * 1.0 * log(1.0) -
        0.5 * 2.0 * log(1.5) +
        loggamma(2.0 / 2) -
        loggamma(1.0 / 2) +
        0.5 * log(1.0 / 2.0) -
        0.5 * log(π)
    )
end

@testset "logpredlik (non-empty clusters)" begin
    N = 2
    y = [1.0, 0.0]
    rng = MersenneTwister(1)
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    m = NormalDPM(rng, N; v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m, y)
    @test m.v1[1] ≈ 3.0
    @test m.r1[1] ≈ 3.0
    @test m.u1[1] ≈ 1/3
    @test m.s1[1] ≈ 5/3
    @test DPMNeal3.logpredlik(m, y, 2, 1) ≈ (
        0.5 * 2 * log(1.5) -
        0.5 * 3 * log(5/3) +
        loggamma(3/2) -
        loggamma(2/2) +
        0.5 * log(2/3) -
        0.5 * log(π)
    )
end

@testset "update! (1)" begin
    N = 2
    y = [1.0, 0.0]
    rng = MersenneTwister(1)
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    m = NormalDPM(rng, N; v0, r0, u0, s0)
    update!(rng, m, y)
end

@testset "update! (2)" begin
    N = 1000
    rng = MersenneTwister(1)
    m = NormalDPM(rng, N)
    y = randn(rng, N)
    update!(rng, m, y)
end

# @testset "update_γ!" begin
#     rng = MersenneTwister(1)
#     N, F = 1000, 1
#     y = randn(rng, N)
#     x = [rand(rng, 1:3, F) for _ in 1:N]
#     x = StatsBase.denserank(x)
#     G = length(unique(x))
#     data = Data(x, y)
#     sb = SpecificBlock(G)
#     gb = DPMGB(rng, N)
#     update!(rng, sb, gb, data)
#     update_γ!(rng, sb, gb, data)
# end

# @testset "final_example" begin
#     rng = MersenneTwister(1)
#     N, F = 1000, 1
#     y = randn(rng, N)
#     x = [rand(rng, 1:3, F) for _ in 1:N]
#     x = StatsBase.denserank(x)
#     for i = 1:N
#         if x[i] == 3
#             y[i] += 10.0
#         end
#     end
#     y .= (y .- mean(y)) ./ √var(y)
#     G = length(unique(x))
#     data = Data(x, y)
#     sb = SpecificBlock(G)
#     gb = DPMGB(rng, N; K0 = 1)
#     for t in 1:10
#         update!(rng, sb, gb, data)
#         update_γ!(rng, sb, gb, data)
#         println(sb.γ[:])
#     end
# end
