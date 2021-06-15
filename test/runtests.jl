using Parameters
using DPMNeal3
using Random
using Test
# using StatsBase
# using Statistics
using SpecialFunctions: loggamma
import DPMNeal3: parent_dpm, logpredlik, update_hyperpars!, update_suffstats!
include("normaldpm.jl")

struct Data
    x::Vector{Int}
    data::Vector{Float64}
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

@testset "DPM Accessors" begin
    rng = MersenneTwister(1)
    m = DPM(rng, 5; K0 = 5)
    @test n_clusters(m) == 5
    @test cluster_capacity(m) == 6
    @test cluster_labels(m) == collect(1:5)
    @test cluster_sizes(m) == [ones(Int, 5); 0]
    @test active_clusters(m) == Set(1:5)
    @test passive_clusters(m) == Set(6)
    @test dp_mass(m) < Inf
    @test dp_mass(m) > 0.0
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
    add_cluster!(m)
    @test length(m.v1) == 2
    @test length(m.r1) == 2
    @test length(m.u1) == 2
    @test length(m.s1) == 2
end

@testset "update_suffstats! (1)" begin
    N = 1
    data = ones(N)
    rng = MersenneTwister(1)
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    m = NormalDPM(rng, N; v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m, data)
    @test m.v1[1] ≈ 2.0
    @test m.v1[2] ≈ 1.0
    @test m.r1[1] ≈ 2.0
    @test m.r1[2] ≈ 1.0
    @test m.u1[1] ≈ 0.5
    @test m.u1[2] ≈ 0.0
    @test m.s1[1] ≈ 1.5
    @test m.s1[2] ≈ 1.0
end

@testset "update_suffstats! (2)" begin
    N = 1
    data = ones(N)
    rng = MersenneTwister(1)
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    m = NormalDPM(rng, N; v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m, data)
    DPMNeal3.update_suffstats!(m, data, 1, 1, 2)
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
    data = ones(N)
    rng = MersenneTwister(1)
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    m = NormalDPM(rng, N; v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m, data)
    @test DPMNeal3.logpredlik(m, data, 1, first(passive_clusters(m))) ≈ (
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
    data = [1.0, 0.0]
    rng = MersenneTwister(1)
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    m = NormalDPM(rng, N; v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m, data)
    @test m.v1[1] ≈ 3.0
    @test m.r1[1] ≈ 3.0
    @test m.u1[1] ≈ 1/3
    @test m.s1[1] ≈ 5/3
    @test DPMNeal3.logpredlik(m, data, 2, 1) ≈ (
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
    data = [1.0, 0.0]
    rng = MersenneTwister(1)
    v0, r0, u0, s0 = 1.0, 1.0, 0.0, 1.0
    m = NormalDPM(rng, N; v0, r0, u0, s0)
    update!(rng, m, data)
end

@testset "update! (2)" begin
    N = 1000
    rng = MersenneTwister(1)
    m = NormalDPM(rng, N)
    data = randn(rng, N)
    update!(rng, m, data)
end
