@testset "Skeleton" begin
    N, a0, b0 = 10, 1.0, 2.0
    d = τ = collect(1:N)
    α = Ref(1.0)
    skl = DPMNeal3.Skeleton(; N)
    skl = DPMNeal3.Skeleton(; N, a0, b0, d, τ, α)
    @test_throws DimensionMismatch DPMNeal3.Skeleton(; N = 1, d = zeros(Int, 0))
    @test_throws DimensionMismatch DPMNeal3.Skeleton(; N = 1, τ = zeros(Int, 0))
    @test_throws ArgumentError DPMNeal3.Skeleton(; N = 0, d = zeros(Int, 0))
    @test_throws DomainError DPMNeal3.Skeleton(; N = 1, d = [-1])
    @test_throws DomainError DPMNeal3.Skeleton(; N = 1, d = [+2])
    @test_throws DomainError DPMNeal3.Skeleton(; N = 2, d = [2, 3])
    @test_throws DomainError DPMNeal3.Skeleton(; N = 2, d = [1, 3])
    @test_throws DomainError DPMNeal3.Skeleton(; N = 2, τ = [1, 3])
    @test_throws DomainError DPMNeal3.Skeleton(; N = 2, τ = [2, 3])
    @test_throws DomainError DPMNeal3.Skeleton(; N = 2, τ = [2, -1])
end

@testset "n_clusters" begin
    N = 100
    d = Distributions.sample(1:5, N)
    skl = DPMNeal3.Skeleton(; N, d)
    @test DPMNeal3.n_clusters(skl) == skl.K[]
    @test DPMNeal3.n_clusters(skl) == length(unique(d))
end

@testset "cluster_sizes" begin
    N = 100
    d = Distributions.sample(1:5, N)
    n = [sum(d .== i) for i in 1:maximum(d)]
    skl = DPMNeal3.Skeleton(; N, d)
    @test DPMNeal3.cluster_sizes(skl) == n
    @test DPMNeal3.cluster_sizes(skl) !== skl.n
    @test DPMNeal3.cluster_sizes(skl; make_copy = false) === skl.n
end

@testset "active_clusters" begin
    N = 100
    d = Distributions.sample(1:5, N)
    A = Set(unique(d))
    skl = DPMNeal3.Skeleton(; N, d)
    @test DPMNeal3.active_clusters(skl) == A
    @test DPMNeal3.active_clusters(skl) !== skl.A
    @test DPMNeal3.active_clusters(skl; make_copy = false) === skl.A
end

@testset "passive_clusters" begin
    N = 100
    d = Distributions.sample(1:5, N)
    P = Set(maximum(d) + 1)
    skl = DPMNeal3.Skeleton(; N, d)
    @test DPMNeal3.passive_clusters(skl) == P
    @test DPMNeal3.passive_clusters(skl) !== skl.P
    @test DPMNeal3.passive_clusters(skl; make_copy = false) === skl.P
end

@testset "cluster_labels" begin
    N = 100
    d = Distributions.sample(1:5, N)
    skl = DPMNeal3.Skeleton(; N, d)
    @test DPMNeal3.cluster_labels(skl) == d
    @test DPMNeal3.cluster_labels(skl) !== d
    @test DPMNeal3.cluster_labels(skl; make_copy = false) === d
end
