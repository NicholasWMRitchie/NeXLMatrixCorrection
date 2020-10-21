using Test
using NeXLMatrixCorrection

#@testset "XPhi" begin
    @testset "Mg in Al at 25 keV" begin
        # See Figure 1 of Merlet 1994
        m, cxr, e0, toa = mat"Al", n"Mg K-L3", 25.0e3, deg2rad(40.0)
        xp = matrixcorrection(XPhi, m, inner(cxr), e0)

        @test isapprox(NeXLMatrixCorrection.ϕ0(xp), 1.4, atol=0.1)
        @test isapprox(NeXLMatrixCorrection.ϕ0(xp), ϕ(xp, 0.0), rtol=1.0e-5)
        @test isapprox(max(xp), 2.65, atol=0.1)
        @test isapprox(ϕ(xp, 1.0e-3), 0.8, atol=0.1)
        @test isapprox(ϕ(xp, 1.6e-3), 0.1, atol=0.1)
    end
    @testset "Cd in Al at 25 keV" begin
        m, cxr, e0, toa = mat"Al", n"Cd L3-M5", 25.0e3, deg2rad(40.0)
        xp = matrixcorrection(XPhi, m, inner(cxr), e0)

        @test isapprox(ϕ(xp, 0.0), 1.5, atol=0.05)
        @test isapprox(ϕ(xp, xp.ρzm), 2.45, atol=0.1)
        @test isapprox(ϕ(xp, 1.0e-3), 0.75, atol=0.1)
        @test isapprox(ϕ(xp, 1.8e-3), 0.2, atol=0.1)
    end
    @testset "Cd in Au at 25 keV" begin
        m, cxr, e0, toa = mat"Au", n"Cd L3-M5", 25.0e3, deg2rad(40.0)
        xp = matrixcorrection(XPhi, m, inner(cxr), e0)

        @test isapprox(ϕ(xp, 0.0), 1.5, atol=0.05)
        @test isapprox(ϕ(xp, xp.ρzm), 3.0, atol=0.1)
        @test isapprox(ϕ(xp, 1.0e-3), 0.2, atol=0.1)
    end
    @testset "Al in C at 15 keV" begin
        m, cxr, e0, toa = mat"C", n"Al K-L3", 15.0e3, deg2rad(40.0)
        xp = matrixcorrection(XPhi, m, inner(cxr), e0)

        @test isapprox(ϕ(xp, 0.0), 1.1, atol=0.05)
        @test isapprox(ϕ(xp, xp.ρzm), 1.70, atol=0.1)
        @test isapprox(ϕ(xp, 0.3e-3), 0.75, atol=0.1)
        @test isapprox(ϕ(xp, 0.6e-3), 0.1, atol=0.1)
    end
#end
