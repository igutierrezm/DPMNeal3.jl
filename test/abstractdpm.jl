struct DummyDPM <: DPMNeal3.AbstractDPM
    skl::DPMNeal3.Skeleton
end

DPMNeal3.skeleton(m::DummyDPM) = m.skl
DPMNeal3.update_suffstats!(m::DummyDPM) = nothing
DPMNeal3.update_suffstats!(m::DummyDPM, i::Int, k0::Int, k1::Int) = nothing
DPMNeal3.logpredlik(m::DummyDPM, i::Int, k::Int) = 1.0

@testset "update_mass_parameter!" begin
    N = 10
    skl = DPMNeal3.Skeleton(; N)
    dpm = DummyDPM(skl)
    DPMNeal3.update_mass_parameter!(dpm)
    @test dpm.skl.α[] > 0
end

@testset "update_cluster_indicators!" begin
    N = 10
    skl = DPMNeal3.Skeleton(; N)
    dpm = DummyDPM(skl)
    DPMNeal3.update_cluster_indicators!(dpm)
    @test Set(unique(skl.d)) == skl.A
    @test !isempty(skl.P)
    @test minimum(skl.d) > 0
    @test isempty(intersect(skl.A, skl.P))
end
