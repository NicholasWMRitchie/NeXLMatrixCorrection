using Test
using NeXLMatrixCorrection

@testset "Standard" begin
        fe1 = KRatio(
            characteristic(n"Fe", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.6517, 0.02),
        ) 
        fe2 = KRatio(
            characteristic(n"Fe", ltransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.3217, 0.02),
        )
        ca1 = KRatio(
            characteristic(n"Ca", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"CaF2",
            uv(0.7280, 0.02),
        )
        ca2 = KRatio(
            characteristic(n"Ca", kbeta),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"CaF2",
            uv(0.7326, 0.02),
        )
        si1 = KRatio(
            characteristic(n"Si", ktransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Si",
            uv(0.3801, 0.02),
        )
        festd = Standard(mat"Fe2O3", [ fe1, fe2 ])
        castd = Standard(mat"Ca5(PO4)3F", [ ca1, ca2 ])
        sistd = Standard(mat"SiO2", [ si1 ])
        @test NeXLCore.findmatch(KRatio(
            characteristic(n"Fe", ltransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.6517, 0.02),
            ), festd) == fe2
        @test isnothing(NeXLCore.findmatch(KRatio(
            characteristic(n"Fe", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.45, 0.02),
        ), castd))
        sfek = standardize(KRatio(
            characteristic(n"Fe", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.45, 0.02),
        ),[festd, castd, sistd])
        @test value(sfek) ==0.45/0.6517
        @test isequal(sfek.standard, mat"Fe2O3")
        scak=standardize(KRatio(
            characteristic(n"Ca", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"CaF2",
            uv(0.32, 0.02),
        ),[festd, castd, sistd])
        @test value(scak)==0.32/0.7280
        @test isequal(scak.standard, mat"Ca5(PO4)3F")
        nos = standardize(KRatio(
            characteristic(n"F", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"LiF",
            uv(0.333, 0.02),
        ), [festd, castd, sistd])
        @test isequal(nos.standard, mat"LiF")
        @test value(nos)==0.333
    end