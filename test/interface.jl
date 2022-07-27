struct DummyDPM2 <: DPMNeal3.AbstractDPM
end

@testset "skeleton" begin
    m = DummyDPM2()
    @test_throws ErrorException DPMNeal3.skeleton(m)
end

@testset "update_suffstats!" begin
    m = DummyDPM2()
    @test_throws ErrorException DPMNeal3.update_suffstats!(m)
    @test_throws ErrorException DPMNeal3.update_suffstats!(m, 1, 2, 3)
end

@testset "logpredlik" begin
    m = DummyDPM2()
    @test_throws ErrorException DPMNeal3.logpredlik(m, 1, 2)
end
