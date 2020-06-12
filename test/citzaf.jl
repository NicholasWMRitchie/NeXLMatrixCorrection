using Test

@testset "CitZAF" begin
    # Mostly the testing is performed in testagainstheinrich.jmd
    # These tests just ensure stablity.
    θ = deg2rad(40.0)
    algs = zafcorrection(CitZAF, ReedFluorescence, NullCoating, mat"NaAlSi3O8", mat"NaF", n"Na K", 15.0e3)
    @test isapprox(A(algs..., n"Na K-L3",θ,θ), 1.0404, atol=0.0001)
    @test isapprox(Z(algs...), 0.96463, atol=0.0001)
    algs = zafcorrection(CitZAF, ReedFluorescence, NullCoating, mat"NaAlSi3O8", mat"Al2O3", n"Al K", 15.0e3)
    @test isapprox(A(algs..., n"Al K-L3",θ,θ), 0.89712, atol=0.0001)
    @test isapprox(Z(algs...), 0.9896, atol=0.0001)
    algs = zafcorrection(CitZAF, ReedFluorescence, NullCoating, mat"NaAlSi3O8", mat"Al2O3", n"O K", 15.0e3)
    @test isapprox(A(algs..., n"O K-L3",θ,θ), 0.94331, atol=0.0001)
    @test isapprox(Z(algs...), 0.98962, atol=0.0001)
    algs = zafcorrection(CitZAF, ReedFluorescence, NullCoating, mat"NaAlSi3O8", mat"SiO2", n"Si K", 15.0e3)
    @test isapprox(A(algs..., n"Si K-L3",θ,θ), 0.84532, atol=0.0001)
    @test isapprox(Z(algs...), 1.0065, atol=0.0001)
end
