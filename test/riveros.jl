using Test

@testset "Riveros" begin
    # Mostly the testing is performed in testagainstheinrich.jmd
    # These tests just ensure stablity.
    θ = deg2rad(40.0)
    algs = zafcorrection(Riveros1993, ReedFluorescence, NullCoating, mat"NaAlSi3O8", mat"NaF", n"Na K", 15.0e3)
    @test isapprox(A(algs..., n"Na K-L3",θ,θ), 1.0401, atol=0.0001)
    @test isapprox(Z(algs...), 0.96363, atol=0.0001)
    algs = zafcorrection(Riveros1993, ReedFluorescence, NullCoating, mat"NaAlSi3O8", mat"Al2O3", n"Al K", 15.0e3)
    @test isapprox(A(algs..., n"Al K-L3",θ,θ), 0.92632, atol=0.0001)
    @test isapprox(Z(algs...), 0.9900, atol=0.0001)
    algs = zafcorrection(Riveros1993, ReedFluorescence, NullCoating, mat"NaAlSi3O8", mat"Al2O3", n"O K", 15.0e3)
    @test isapprox(A(algs..., n"O K-L3",θ,θ), 0.94999, atol=0.0001)
    @test isapprox(Z(algs...), 0.9901, atol=0.0001)
    algs = zafcorrection(Riveros1993, ReedFluorescence, NullCoating, mat"NaAlSi3O8", mat"SiO2", n"Si K", 15.0e3)
    @test isapprox(A(algs..., n"Si K-L3",θ,θ), 0.89044, atol=0.0001)
    @test isapprox(Z(algs...), 1.0067, atol=0.0001)
end
