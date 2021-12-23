# using Revise
using DPMNeal3
using Test
using SpecialFunctions: loggamma

# import DPMNeal3: parent_dpm, logpredlik, update_hyperpars!, update_suffstats!
# include("normaldpm.jl")

@testset "Skeleton" begin
    N, M, a_α0, b_α0 = 10, 5, 1.0, 2.0
    skl = DPMNeal3.Skeleton(N = N, M = M)
    skl = DPMNeal3.Skeleton(; N, M, a_α0, b_α0)
    @test unique(DPMNeal3.cluster_labels(skl)) == [1]
    @test unique(DPMNeal3.cluster_sizes(skl)) == [N]
    @test length(DPMNeal3.predictive_pdf(skl)) == M
    @test length(DPMNeal3.cluster_labels(skl)) == N
    @test length(DPMNeal3.cluster_sizes(skl)) == 1
    @test DPMNeal3.n_clusters(skl) == 1
    @test DPMNeal3.active_clusters(skl) == Set(1)
    @test DPMNeal3.passive_clusters(skl) == Set(2)
end

@testset "DPMNormal - add_cluster!" begin
    m = DPMNeal3.DPMNormal(y = randn(10))
    DPMNeal3.add_cluster!(m)
    @test length(m.v1) == 2
    @test length(m.r1) == 2
    @test length(m.u1) == 2
    @test length(m.s1) == 2
end

@testset "DPMNormal - update_suffstats! (1)" begin
    y, v0, r0, u0, s0 = ones(1), 1.0, 1.0, 0.0, 1.0
    m = DPMNeal3.DPMNormal(; y, v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m)
    @test m.v1[1] ≈ 2.0
    @test m.v1[2] ≈ 1.0
    @test m.r1[1] ≈ 2.0
    @test m.r1[2] ≈ 1.0
    @test m.u1[1] ≈ 0.5
    @test m.u1[2] ≈ 0.0
    @test m.s1[1] ≈ 1.5
    @test m.s1[2] ≈ 1.0
end

@testset "DPMNormal - update_suffstats! (2)" begin
    y, v0, r0, u0, s0 = ones(1), 1.0, 1.0, 0.0, 1.0
    m = DPMNeal3.DPMNormal(; y, v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m)
    DPMNeal3.update_suffstats!(m, 1, 1, 2)
    @test m.v1[1] ≈ 1.0
    @test m.v1[2] ≈ 2.0
    @test m.r1[1] ≈ 1.0
    @test m.r1[2] ≈ 2.0
    @test m.u1[1] ≈ 0.0
    @test m.u1[2] ≈ 0.5
    @test m.s1[1] ≈ 1.0
    @test m.s1[2] ≈ 1.5
end

@testset "DPMNormal - in_sample_logpredlik (empty clusters)" begin
    y, v0, r0, u0, s0 = ones(1), 1.0, 1.0, 0.0, 1.0
    m = DPMNeal3.DPMNormal(; y, v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m)
    @test DPMNeal3.in_sample_logpredlik(m, 1, 2) ≈ (
        0.5 * 1.0 * log(1.0) -
        0.5 * 2.0 * log(1.5) +
        loggamma(2.0 / 2) -
        loggamma(1.0 / 2) +
        0.5 * log(1.0 / 2.0) -
        0.5 * log(π)
    )
end

@testset "DPMNormal - logpredlik (non-empty clusters)" begin
    y, v0, r0, u0, s0 = [1.0, 0.0], 1.0, 1.0, 0.0, 1.0
    m = DPMNeal3.DPMNormal(; y, v0, r0, u0, s0)
    DPMNeal3.update_suffstats!(m)
    @test m.v1[1] ≈ 3.0
    @test m.r1[1] ≈ 3.0
    @test m.u1[1] ≈ 1/3
    @test m.s1[1] ≈ 5/3
    @test DPMNeal3.in_sample_logpredlik(m, 2, 1) ≈ (
        0.5 * 2 * log(1.5) -
        0.5 * 3 * log(5/3) +
        loggamma(3/2) -
        loggamma(2/2) +
        0.5 * log(2/3) -
        0.5 * log(π)
    )
end

@testset "DPMNormal - update! (1)" begin
    y, v0, r0, u0, s0 = [1.0, 0.0], 1.0, 1.0, 0.0, 1.0
    m = DPMNeal3.DPMNormal(; y, v0, r0, u0, s0)
    DPMNeal3.update!(m)
end

@testset "DPMNormal - update! (2)" begin
    y, v0, r0, u0, s0 = randn(100), 1.0, 1.0, 0.0, 1.0
    m = DPMNeal3.DPMNormal(; y, v0, r0, u0, s0)
    DPMNeal3.update!(m)
end

@testset "DPMNormal - Check that the interface must be implemented" begin
    struct MyDPM <: DPMNeal3.AbstractModel end
    m = MyDPM()
    @test_throws ErrorException DPMNeal3.skeleton(m)
    @test_throws ErrorException DPMNeal3.update_suffstats!(m)
    @test_throws ErrorException DPMNeal3.update_suffstats!(m, 1, 1, 2)
    @test_throws ErrorException DPMNeal3.in_sample_logpredlik(m, 1, 1)
    @test_throws ErrorException DPMNeal3.out_of_sample_logpredlik(m, 1, 1)
end